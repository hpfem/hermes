// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/

#include "../../hermes_common/matrix.h"
#include "discrete_problem.h"
#include "weakform.h"
#include "space.h"

#include "../../hermes_common/error.h"
#include "../../hermes_common/callstack.h"

//  Solvers
#include "../../hermes_common/solver/solver.h"


static int _precalculated = 0;

DiscreteProblem::DiscreteProblem(WeakForm* wf, Space* space, bool is_linear) : wf(wf), 
space(space), is_linear(is_linear)
{
  if(space->get_n_eq() != wf->get_neq())
        error("WeakForm does not have as many equations as Space in DiscreteProblem::DiscreteProblem()");
  if (_precalculated == 0) {
    // precalculating values and derivatives 
    // of all polynomials at all possible 
    // integration points
    info("Precalculating Legendre polynomials...");
    precalculate_legendre_1d();
    precalculate_legendre_1d_left();
    precalculate_legendre_1d_right();

    info("Precalculating Lobatto shape functions...");
    precalculate_lobatto_1d();
    precalculate_lobatto_1d_left();
    precalculate_lobatto_1d_right();
    _precalculated = 1;
  }
}

// process volumetric weak forms
void DiscreteProblem::process_vol_forms(SparseMatrix *mat, Vector *rhs, bool rhsonly) {
  int n_eq = space->get_n_eq();
  Element *elems = space->get_base_elems();
  int n_elem = space->get_n_base_elem();
  Iterator *I = new Iterator(space);

  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    //printf("Processing elem %d\n", m);
    int    pts_num;                                     // num of quad points
    double phys_pts[MAX_QUAD_PTS_NUM];                  // quad points
    double phys_weights[MAX_QUAD_PTS_NUM];              // quad weights
    double phys_u[MAX_QUAD_PTS_NUM];                    // basis function 
    double phys_dudx[MAX_QUAD_PTS_NUM];                 // basis function x-derivative
    double phys_v[MAX_QUAD_PTS_NUM];                    // test function
    double phys_dvdx[MAX_QUAD_PTS_NUM];                 // test function x-derivative
    if (n_eq > MAX_EQN_NUM) error("number of equations exceeded in process_vol_forms().");
    // all previous solutions (all components)
    double phys_u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM];     
    // x-derivatives of all previous solutions (all components)
    double phys_du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM];  
    // decide quadrature order and set up 
    // quadrature weights and points in element m
    // CAUTION: This is heuristic
    int order = 4*e->p;

    // prepare quadrature points and weights in element 'e'
    create_phys_element_quadrature(e->x1, e->x2,  
                               order, phys_pts, phys_weights, &pts_num); 

    // evaluate previous solution and its derivative 
    // at all quadrature points in the element, 
    // for every solution component
    // 0... in the entire element
    for(int sln=0; sln < e->n_sln; sln++) {
      e->get_solution_quad(0, order, phys_u_prev[sln], phys_du_prevdx[sln], sln); 
    }

    // volumetric bilinear forms
    for (unsigned int ww = 0; ww < this->wf->matrix_forms_vol.size(); ww++) {
	    WeakForm::MatrixFormVol *mfv = &this->wf->matrix_forms_vol[ww];
	    if (e->marker == mfv->marker ||  mfv->marker == ANY) {
  	    int c_i = mfv->i;  
	      int c_j = mfv->j;  

	      // loop over test functions (rows)
	      for(int i=0; i<e->p + 1; i++) {
	        // if i-th test function is active
	        int pos_i = e->dof[c_i][i]; // row in matrix
          //printf("elem (%g, %g): pos_i = %d\n", e->x1, e->x2, pos_i);
	        if(pos_i != -1) {
	          // transform i-th test function to element 'm'
            //printf("Elem (%g, %g): i = %d, order = %d\n", e->x1, e->x2, i, order);
	          element_shapefn(e->x1, e->x2, i, order, phys_v, phys_dvdx); 
	          // if we are constructing the matrix
            // loop over basis functions (columns)
            for(int j=0; j < e->p + 1; j++) {
	            int pos_j = e->dof[c_j][j]; // matrix column
              //printf("elem (%g, %g): pos_j = %d\n", e->x1, e->x2, pos_j);
	            // if j-th basis function is active
	              // transform j-th basis function to element 'm'
              element_shapefn(e->x1, e->x2,  
		          j, order, phys_u, phys_dudx); 
              // evaluate the bilinear form
              double val_ij = mfv->fn(pts_num, phys_pts,
	            phys_weights, phys_u, phys_dudx, phys_v, phys_dvdx,
              phys_u_prev, phys_du_prevdx, mfv->space); 
              //truncating
              if (fabs(val_ij) < 1e-12) val_ij = 0.0; 
              // add the result to the matrix
              if (val_ij != 0) {
                if(pos_j != -1) {
                  if(!rhsonly) {
                    mat->add(pos_i, pos_j, val_ij);
                    if (DEBUG_MATRIX) {
                      info("Adding to matrix pos %d, %d value %g (comp %d, %d)", 
                      pos_i, pos_j, val_ij, c_i, c_j);
                    }
                  }
                }
                else
                  if(this->is_linear && rhs != NULL)
                    rhs->add(pos_i, -val_ij * e->coeffs[c_j][c_j][j]);
            }
            }
          }
        }
      }
    }

    // volumetric part of residual
    for (unsigned int ww = 0; ww < this->wf->vector_forms_vol.size(); ww++) {
      WeakForm::VectorFormVol *vfv = &this->wf->vector_forms_vol[ww];
      int c_i = vfv->i;  
      // loop over test functions (rows)
      for(int i=0; i<e->p + 1; i++) {
        // if i-th test function is active
        int pos_i = e->dof[c_i][i]; // row in residual vector
        if(pos_i != -1) {
          // transform i-th test function to element 'm'
          element_shapefn(e->x1, e->x2,  
		      i, order, phys_v, phys_dvdx); 
          // contribute to residual vector
          double val_i = vfv->fn(pts_num, phys_pts, phys_weights, 
			       phys_u_prev, phys_du_prevdx, phys_v,
             phys_dvdx, vfv->space);
          // truncating
          if(fabs(val_i) < 1e-12) val_i = 0.0; 
          // add the contribution to the residual vector
          if (val_i != 0 && rhs != NULL) rhs->add(pos_i, val_i);
          if (DEBUG_MATRIX)
	          if (val_i != 0) {
	            info("Adding to residual pos %d value %g (comp %d)", 
                  pos_i, val_i, c_i);
              }
        }
      }
    }
  } // end while
  delete I;
}

// process boundary weak forms
void DiscreteProblem::process_surf_forms(SparseMatrix *mat, Vector *rhs, int bdy_index, bool rhsonly) {
  Iterator *I = new Iterator(space);
  Element *e; 

  // evaluate previous solution and its derivative at the end point
  // FIXME: maximum number of equations limited by [MAX_EQN_NUM][MAX_QUAD_PTS_NUM]
  double phys_u_prev[MAX_SLN_NUM][MAX_EQN_NUM], 
         phys_du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM]; // at the end point

  // decide whether we are on the left-most or right-most one
  double x_ref, x_phys; 
  if(bdy_index == BOUNDARY_LEFT) {
    e = I->first_active_element(); 
    x_ref = -1; // left end of reference element
    x_phys = space->get_left_endpoint();
  }
  else {
    e = I->last_active_element(); 
    x_ref = 1;  // right end of reference element
    x_phys = space->get_right_endpoint();
  }

  // get solution value and derivative at the boundary point
  for(int sln=0; sln < e->n_sln; sln++) {
    e->get_solution_point(x_phys, phys_u_prev[sln], phys_du_prevdx[sln], sln); 
  }

  // surface bilinear forms
  for (unsigned int ww = 0; ww < this->wf->matrix_forms_surf.size(); ww++) {
    WeakForm::MatrixFormSurf *mfs = &this->wf->matrix_forms_surf[ww];
    if (mfs->bdy_index != bdy_index) continue;
    int c_i = mfs->i;  
    int c_j = mfs->j;  

    // loop over test functions on the boundary element
    for(int i=0; i<e->p + 1; i++) {
      double phys_v, phys_dvdx; 
      int pos_i = e->dof[c_i][i]; // matrix row
      if(pos_i != -1) {
        // transform j-th basis function to the boundary element
        element_shapefn_point(x_ref, e->x1, e->x2, i, phys_v, 
                              phys_dvdx); 
        // loop over basis functions on the boundary element
        for(int j=0; j < e->p + 1; j++) {
          double phys_u, phys_dudx;
          int pos_j = e->dof[c_j][j]; // matrix column
          // if j-th basis function is active
          // transform j-th basis function to the boundary element
          element_shapefn_point(x_ref, e->x1, e->x2, j, phys_u, 
                                phys_dudx); 
          // evaluate the surface bilinear form
          double val_ij_surf = mfs->fn(x_phys,
                           phys_u, phys_dudx, phys_v, 
                           phys_dvdx, phys_u_prev, phys_du_prevdx, 
                           NULL); 
          // truncating
          if(fabs(val_ij_surf) < 1e-12) val_ij_surf = 0.0; 
          // add the result to the matrix
          if (val_ij_surf != 0) {
            if(pos_j != -1) {
              if(!rhsonly) {
                mat->add(pos_i, pos_j, val_ij_surf);
                if (DEBUG_MATRIX) {
                    info("Adding to matrix pos %d, %d value %g (comp %d, %d)", 
                    pos_i, pos_j, val_ij_surf, c_i, c_j);
                }
              }
            }
            else
              if(this->is_linear && rhs != NULL)
                  rhs->add(pos_i, -val_ij_surf);
          }
        }
      }
    }
  }

  // surface part of residual
  for (unsigned int ww = 0; ww < this->wf->vector_forms_surf.size(); ww++)
  {
    WeakForm::VectorFormSurf *vfs = &this->wf->vector_forms_surf[ww];
    if (vfs->bdy_index != bdy_index) continue;
    int c_i = vfs->i;  

    // loop over test functions on the boundary element
    for(int i=0; i<e->p + 1; i++) {
      double phys_v, phys_dvdx; 
      int pos_i = e->dof[c_i][i]; // matrix row
      if(pos_i != -1) {
        // transform j-th basis function to the boundary element
        element_shapefn_point(x_ref, e->x1, e->x2, i, phys_v, 
                              phys_dvdx); 
        // evaluate the surface bilinear form
        double val_i_surf = vfs->fn(x_phys,
                        phys_u_prev, phys_du_prevdx, phys_v, phys_dvdx, 
                        NULL);
        // truncating
        if(fabs(val_i_surf) < 1e-12) val_i_surf = 0.0; 
        // add the result to the matrix
        if (val_i_surf != 0 && rhs != NULL) rhs->add(pos_i, val_i_surf);
      }
    }
  }
  delete I;
}

// construct Jacobi matrix or residual vector
void DiscreteProblem::assemble(scalar *coeff_vec, SparseMatrix *mat, Vector *rhs, bool rhsonly,
                               bool force_diagonal_blocks, bool add_dir_lift, Table* block_weights) {
  // number of equations in the system
  int n_eq = space->get_n_eq();

  // total number of unknowns
  int ndof = Space::get_num_dofs(space);

  // Copy coefficients from vector coeff_vec to elements.
  if (coeff_vec != NULL) set_coeff_vector(coeff_vec, space);

  // Reallocate the matrix and residual vector.
  if (rhs != NULL) rhs->alloc(ndof);
  // Zero the vector, which should be done by the appropriate implementation anyway.
  if (rhs != NULL) rhs->zero();
  if (mat != NULL)
  {
    mat->free();
    mat->prealloc(ndof);
    for(int i = 0; i < ndof; i++)
      for(int j = 0; j< ndof; j++)
        mat->pre_add_ij(i, j);
    mat->alloc();
    // Zero the matrix, which should be done by the appropriate implementation anyway.
    mat->zero();
  }

  // process volumetric weak forms via an element loop
  process_vol_forms(mat, rhs, rhsonly);

  // process surface weak forms for the left boundary
  process_surf_forms(mat, rhs, BOUNDARY_LEFT, rhsonly);

  // process surface weak forms for the right boundary
  process_surf_forms(mat, rhs, BOUNDARY_RIGHT, rhsonly);

  // DEBUG: print Jacobi matrix
  if(DEBUG_MATRIX) {
    info("Jacobi matrix:");
    for(int i=0; i<ndof; i++) {
      for(int j=0; j<ndof; j++) {
        info("%g ", mat->get(i, j));
      }
    }
  }

  // DEBUG: print residual vector
  if(DEBUG_MATRIX && rhs != NULL) {
    info("Residual:");
    for(int i=0; i<ndof; i++) {
      info("%g ", rhs->get(i));
    }
  }
} 

void J_dot_vec_jfnk(DiscreteProblem *dp, Space *space, Vector* vec,
                    scalar* y_orig, Vector* f_orig, 
                    Vector* J_dot_vec,
                    double jfnk_epsilon, int ndof, MatrixSolverType matrix_solver) 
{
  scalar* y_perturbed = new scalar[ndof];
  Vector* f_perturbed = create_vector(matrix_solver);
  for (int i = 0; i < ndof; i++)
    y_perturbed[i] = y_orig[i] + jfnk_epsilon*vec->get(i);
  // NULL stands for that we are not interested in the matrix, just the vector.
  dp->assemble(y_perturbed, NULL, f_perturbed, true); 
  set_coeff_vector(y_orig, space);
  for (int i = 0; i < ndof; i++) {
    J_dot_vec->set(i,(f_perturbed->get(i) - f_orig->get(i))/jfnk_epsilon);
  }
  delete [] y_perturbed;
  delete [] f_perturbed;
}

// CG method adjusted for JFNK
// NOTE: 
void jfnk_cg(DiscreteProblem *dp, Space *space, 
             double matrix_solver_tol, int matrix_solver_maxiter, 
	     double jfnk_epsilon, double tol_jfnk, int jfnk_maxiter, 
             MatrixSolverType matrix_solver, bool verbose)
{
  int ndof = Space::get_num_dofs(space);
  // vectors for JFNK
  Vector* f_orig = create_vector(matrix_solver);
  scalar* y_orig = new scalar[ndof];
  Vector* vec = create_vector(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);

  // vectors for the CG method
  Vector* r = create_vector(matrix_solver);
  Vector* p = create_vector(matrix_solver);
  Vector* J_dot_vec = create_vector(matrix_solver);

  // fill vector y_orig using dof and coeffs arrays in elements
  get_coeff_vector(space, y_orig);

  // JFNK loop
  int jfnk_iter_num = 1;
  while (1) {
    // construct residual vector f_orig corresponding to y_orig
    // (f_orig stays unchanged through the entire CG loop)
    // NULL stands for that we are not interested in the matrix, just the vector.
    dp->assemble(y_orig, NULL, f_orig, true); 

    // calculate L2 norm of f_orig
    double res_norm_squared = 0;
    for(int i = 0; i < ndof; i++) res_norm_squared += f_orig->get(i)*f_orig->get(i);

    // If residual norm less than 'tol_jfnk', break
    if (verbose) {
        info("Residual norm: %.15f", sqrt(res_norm_squared));
    }
    if(res_norm_squared < tol_jfnk*tol_jfnk) break;

    if (verbose) {
        info("JFNK iteration: %d", jfnk_iter_num);
    }

    // right-hand side is negative residual
    // (rhs stays unchanged through the CG loop)
    for(int i = 0; i < ndof; i++) rhs->set(i, -f_orig->get(i));

    // beginning CG method
    // r = rhs - A*vec0 (where the initial vector vec0 = 0)
    for (int i=0; i < ndof; i++) r->set(i, rhs->get(i));
    // p = r
    for (int i=0; i < ndof; i++) p->set(i, r->get(i));

    // CG loop
    int iter_current = 0;
    double tol_current_squared;
    // initializing the solution vector with zero
    for(int i = 0; i < ndof; i++) vec->set(i, 0);
    while (1) {
      J_dot_vec_jfnk(dp, space, p, y_orig, f_orig,
                     J_dot_vec, jfnk_epsilon, ndof, matrix_solver);
      double r_times_r = vec_dot(r, r, ndof);
      double alpha = r_times_r / vec_dot(p, J_dot_vec, ndof); 
      for (int i=0; i < ndof; i++) {
        vec->set(i, vec->get(i) + alpha*p->get(i));
        r->set(i, r->get(i) - alpha*J_dot_vec->get(i));
      }
      /*
      // debug - output of solution vector
      printf("   vec = ");
      for (int i=0; i < ndof; i++) {
        printf("%g ", vec[i]);
      }
      printf("\n");
      */

      double r_times_r_new = vec_dot(r, r, ndof);
      iter_current++;
      tol_current_squared = r_times_r_new;
      if (tol_current_squared < matrix_solver_tol*matrix_solver_tol 
          || iter_current >= matrix_solver_maxiter) break;
      double beta = r_times_r_new/r_times_r;
      for (int i=0; i < ndof; i++) p->set(i, r->get(i) + beta*p->get(i));
    }
    // check whether CG converged
    if (verbose) {
        info("CG (JFNK) made %d iteration(s) (tol = %g)",
           iter_current, sqrt(tol_current_squared));
    }
    if(tol_current_squared > matrix_solver_tol*matrix_solver_tol) {
      error("CG (JFNK) did not converge.");
    }

    // updating vector y_orig by new solution which is in x
    for(int i=0; i<ndof; i++) y_orig[i] = y_orig[i] + vec->get(i);

    jfnk_iter_num++;
    if (jfnk_iter_num >= jfnk_maxiter) {
      error("JFNK did not converge.");
    }
  }

  delete f_orig;
  // copy updated vector y_orig to space
  set_coeff_vector(y_orig, space);
}

Space* construct_refined_space(Space* space, int order_increase)
{
  Space* ref_space = space->replicate();
  Iterator *I = new Iterator(ref_space);
  Element *e;
  while ((e = I->next_active_element()) != NULL) {
    int3 cand = {1, e->p + order_increase, e->p + order_increase};
    e->refine(cand);
    if (cand[0] == 1) 
      ref_space->n_active_elem++; // if hp-refinement
  }
  ref_space->assign_dofs();
  return ref_space;
}

double get_l2_norm(Vector* vec) 
{
  _F_
  double val = 0;
  for (unsigned int i = 0; i < vec->length(); i++)
    val = val + vec->get(i)*vec->get(i);
  return sqrt(std::abs(val));
}

void set_zero(scalar *vec, int length) 
{
  memset(vec, 0, length * sizeof(scalar));
}

int DiscreteProblem::get_num_dofs() { return space->get_num_dofs(); };
  
bool DiscreteProblem::is_matrix_free() { return false; };

// Signature of this function is identical in H1D, H2D, H3D, but it is currently unused in H1D.
void DiscreteProblem::create_sparse_structure(SparseMatrix* matrix, Vector* rhs, bool rhsonly,
                                              bool force_diagonal_blocks, Table* block_weights) { return; };