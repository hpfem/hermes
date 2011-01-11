#define HERMES_REPORT_INFO

// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "hermes2d.h"

void create_stage_wf(double current_time, double time_step, ButcherTable* bt, 
                     DiscreteProblem* dp, WeakForm* stage_wf_right, WeakForm* stage_wf_left) 
{
  // Number of stages.
  int num_stages = bt->get_size();

  // Original weak formulation.
  WeakForm* wf = dp->get_weak_formulation();

  // Extract mesh from (the first space of) the discrete problem.
  Mesh* mesh = dp->get_space(0)->get_mesh();

  // Create a constant Solution to represent the stage time 
  // stage_time = current_time + c_i*time_step.
  // (Temporary workaround. these should be passed as numbers.) 
  Solution** stage_time_sol = new Solution*[num_stages];
  for (int i = 0; i < num_stages; i++) {
    stage_time_sol[i] = new Solution(mesh);
    stage_time_sol[i]->set_const(mesh, current_time + bt->get_C(i)*time_step);   
  }

  // Extracting volume and surface matrix and vector forms from the 
  // original weak formulation.
  if (wf->get_neq() != 1) error("wf->neq != 1 not implemented yet.");
  Hermes::Tuple<WeakForm::MatrixFormVol> mfvol_base = wf->get_mfvol();
  Hermes::Tuple<WeakForm::MatrixFormSurf> mfsurf_base = wf->get_mfsurf();
  Hermes::Tuple<WeakForm::VectorFormVol> vfvol_base = wf->get_vfvol();
  Hermes::Tuple<WeakForm::VectorFormSurf> vfsurf_base = wf->get_vfsurf();

  // Duplicate matrix volume forms, scale them according 
  // to the Butcher's table, enhance them with additional 
  // external solutions, and anter them as blocks to the 
  // new stage Jacobian.
  for (unsigned int m = 0; m < mfvol_base.size(); m++) {
    WeakForm::MatrixFormVol mfv_base = mfvol_base[m];
    for (int i = 0; i < num_stages; i++) {
      for (int j = 0; j < num_stages; j++) {
        WeakForm::MatrixFormVol mfv_ij;
        mfv_ij.i = i;
        mfv_ij.j = j;
        mfv_ij.sym = mfv_base.sym;
        mfv_ij.area = mfv_base.area;
        mfv_ij.fn = mfv_base.fn;
        mfv_ij.ord = mfv_base.ord;
        mfv_ij.ext.copy(mfv_base.ext);
        mfv_ij.scaling_factor = -time_step * bt->get_A(i, j);

        // Add stage_time_sol[i] as an external function to the form.
        mfv_ij.ext.push_back(stage_time_sol[i]);

        // Add the matrix form to the corresponding block of the 
        // stage Jacobian matrix.
        stage_wf_right->add_matrix_form(&mfv_ij);
      }
    }
  }

  // Add mass volumetric forms to diagonal blocks
  for (int i = 0; i < num_stages; i++) {
    WeakForm::MatrixFormVol mfv_ii;
    mfv_ii.i = i;
    mfv_ii.j = i;
    mfv_ii.sym = HERMES_SYM;
    mfv_ii.area = HERMES_ANY;
    mfv_ii.fn = l2_form<double, scalar>;
    mfv_ii.ord = l2_form<Ord, Ord>;
    mfv_ii.ext = Hermes::Tuple<MeshFunction*> ();
    mfv_ii.scaling_factor = 1.0;
    // Add the matrix form to the diagonal block.
    stage_wf_left->add_matrix_form(&mfv_ii);
  }

  // Duplicate matrix surface forms, enhance them with 
  // additional external solutions, and anter them as
  // blocks of the stage Jacobian.
  for (unsigned int m = 0; m < mfsurf_base.size(); m++) {
    WeakForm::MatrixFormSurf mfs_base = mfsurf_base[m];
    for (int i = 0; i < num_stages; i++) {
      for (int j = 0; j < num_stages; j++) {
        WeakForm::MatrixFormSurf mfs_ij;
        mfs_ij.i = i;
        mfs_ij.j = j;
        mfs_ij.area = mfs_base.area;
        mfs_ij.fn = mfs_base.fn;
        mfs_ij.ord = mfs_base.ord;
        mfs_ij.ext.copy(mfs_base.ext);
        mfs_ij.scaling_factor = -time_step * bt->get_A(i, j);

        // Add stage_time_sol[i] as an external function to the form.
        mfs_ij.ext.push_back(stage_time_sol[i]);

        // Add the matrix form to the corresponding block of the 
        // stage Jacobian matrix.
        stage_wf_right->add_matrix_form_surf(&mfs_ij);
      }
    }
  }

  // Duplicate vector volume forms, enhance them with 
  // additional external solutions, and anter them as
  // blocks of the stage residual.
  for (unsigned int m = 0; m < vfvol_base.size(); m++) {
    WeakForm::VectorFormVol vfv_base = vfvol_base[m];
    for (int i = 0; i < num_stages; i++) {
      WeakForm::VectorFormVol vfv_i;
      vfv_i.i = i;
      vfv_i.area = vfv_base.area;
      vfv_i.fn = vfv_base.fn;
      vfv_i.ord = vfv_base.ord;
      vfv_i.ext.copy(vfv_base.ext);
      vfv_i.scaling_factor = -1.0;

      // Add stage_time_sol[i] as an external function to the form.
      vfv_i.ext.push_back(stage_time_sol[i]);

      // Add the matrix form to the corresponding block of the 
      // stage Jacobian matrix.
      stage_wf_right->add_vector_form(&vfv_i);
    }
  }

  // Add mass volumetric vector forms to diagonal blocks.
  for (int i = 0; i < num_stages; i++) {
    WeakForm::VectorFormVol vfv_i;
    vfv_i.i = i;
    vfv_i.area = HERMES_ANY;
    vfv_i.fn = l2_residual_form<double, scalar>;
    vfv_i.ord = l2_residual_form<Ord, Ord>;
    vfv_i.ext = Hermes::Tuple<MeshFunction*> ();
    vfv_i.scaling_factor = 1.0;
    // Add the matrix form to the diagonal block.
    stage_wf_left->add_vector_form(&vfv_i);
  }

  // Duplicate vector surface forms, enhance them with 
  // additional external solutions, and anter them as
  // blocks of the stage residual.
  for (unsigned int m = 0; m < vfsurf_base.size(); m++) {
    WeakForm::VectorFormSurf vfs_base = vfsurf_base[m];
    for (int i = 0; i < num_stages; i++) {
      WeakForm::VectorFormSurf vfs_i;
      vfs_i.i = i;
      vfs_i.area = vfs_base.area;
      vfs_i.fn = vfs_base.fn;
      vfs_i.ord = vfs_base.ord;
      vfs_i.ext.copy(vfs_base.ext);
      vfs_i.scaling_factor = -1.0;

      // Add stage_time_sol[i] as an external function to the form.
      vfs_i.ext.push_back(stage_time_sol[i]);

      // Add the matrix form to the corresponding block of the 
      // stage Jacobian matrix.
      stage_wf_right->add_vector_form_surf(&vfs_i);
    }
  }
}

bool HERMES_RESIDUAL_AS_VECTOR_RK = false;
bool rk_time_step(double current_time, double time_step, ButcherTable* const bt, 
                  scalar* coeff_vec, DiscreteProblem* dp, MatrixSolverType matrix_solver, 
                  double newton_tol, int newton_max_iter, bool verbose, 
                  double newton_damping_coeff, double newton_max_allowed_residual_norm)
{
  // Matrix and vector for the time derivative part of 
  // the equation (left-hand side).
  SparseMatrix* matrix_left = create_matrix(matrix_solver);
  Vector* vector_left = create_vector(matrix_solver);

  // Matrix and vector for the stationary part of the 
  // equation (right-hand side).
  SparseMatrix* matrix_right = create_matrix(matrix_solver);
  Vector* vector_right = create_vector(matrix_solver);

  // The two matrices and two vectors will be merged into "matrix" 
  // and "rhs", respectively, so we create just one matrix solver.
  Solver* solver = create_linear_solver(matrix_solver, matrix_right, vector_right);

  // Get number of stages from the Butcher's table.
  int num_stages = bt->get_size();

  // Get original space, mesh, and ndof.
  Space* space = dp->get_space(0);
  Mesh* mesh = space->get_mesh();
  int ndof = space->get_num_dofs();

  // Create spaces for stage solutions. This is necessary 
  // to define a num_stages x num_stages block weak formulation. 
  Hermes::Tuple<Space*> stage_spaces;
  stage_spaces.push_back(dp->get_space(0));
  for (int i = 1; i < num_stages; i++) {
    stage_spaces.push_back(space->dup(mesh));
    stage_spaces[i]->copy_orders(space);
  }

  // Create a multistage weak formulation.
  WeakForm stage_wf_left(num_stages);           // For the time derivative term (written on the left).
  WeakForm stage_wf_right(num_stages);        // For the rest of equation (written on the right).
  create_stage_wf(current_time, time_step, bt, dp, &stage_wf_right, &stage_wf_left); 

  // Create the weak forms for the left- and right-hand sides of 
  // the equation, respectively.
  bool is_linear = dp->get_is_linear();
  DiscreteProblem stage_dp_left(&stage_wf_left, stage_spaces, is_linear);
  DiscreteProblem stage_dp_right(&stage_wf_right, stage_spaces, is_linear);

  // Create stage vector of length num_stages * ndof.
  // It contains coefficients of all stage solutions,
  // to be passed into the time derivative part of the 
  // discrete problem.
  scalar* stage_vec = new scalar[num_stages*ndof];
  memset(stage_vec, 0, num_stages * ndof * sizeof(scalar));

  // Create u_prev vector to be passed to weak forms of
  // the stationary part of the discrete problem.
  scalar* u_prev_vec = new scalar[num_stages*ndof];

  // Prepare residuals of stage solutions.
  Hermes::Tuple<Solution*> residuals;
  Hermes::Tuple<bool> add_dir_lift;
  for (int i = 0; i < num_stages; i++) {
    residuals.push_back(new Solution(mesh));
    add_dir_lift.push_back(false);
  }

  // The Newton's loop.
  double residual_norm;
  int it = 1;
  while (true)
  {
    // Prepare the u_prev vector for weak forms.
    for (int r = 0; r < num_stages; r++) {
      for (int i = 0; i < ndof; i++) {
        double increment_r = 0;
        for (int s = 0; s < num_stages; s++) {
          increment_r += time_step * bt->get_A(r, s) * stage_vec[s*ndof + i]; 
        }
        u_prev_vec[r*ndof + i] = coeff_vec[i] + increment_r;
      }
    } 

    // Assemble the stage Jacobian matrix and residual vector.
    bool rhs_only = false;
    stage_dp_left.assemble(stage_vec, matrix_left, vector_left, rhs_only);
    bool force_diagonal_blocks = true;
    // Sparsity structure is forced so that matrix_left can be added later.
    stage_dp_right.assemble(u_prev_vec, matrix_right, vector_right, 
                            rhs_only, force_diagonal_blocks);

    /*
    // Debug.
    FILE* f = fopen("debug-1.txt", "w");
    matrix_left->dump(f, "tmp", DF_MATLAB_SPARSE); 
    info("Matrix timedep dumped.");
    fclose(f);
    f = fopen("debug-2.txt", "w");
    matrix_right->dump(f, "tmp", DF_MATLAB_SPARSE); 
    info("Matrix stationary dumped.");
    fclose(f);
    */

    // Putting the two parts together into matrix and rhs.
    ((UMFPackMatrix*)matrix_right)->add_umfpack_matrix((UMFPackMatrix*)matrix_left);
    vector_right->add_vector(vector_left);

    // Debug.
    //f = fopen("debug-3.txt", "w");
    //matrix_right->dump(f, "tmp", DF_MATLAB_SPARSE); 
    //info("Merged matrix dumped.");
    //fclose(f);
    //exit(0);

    // Multiply the residual vector with -1 since the matrix 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    vector_right->change_sign();
    
    // Measure the residual norm.
    if (HERMES_RESIDUAL_AS_VECTOR_RK) {
      // Calculate the l2-norm of residual vector.
      residual_norm = get_l2_norm(vector_right);
    }
    else {
      // Translate residual vector into residual functions.
      Solution::vector_to_solutions(vector_right, stage_dp_right.get_spaces(), 
                                    residuals, add_dir_lift);
      residual_norm = calc_norms(residuals);
    }

    // Info for the user.
    if (verbose) info("---- Newton iter %d, ndof %d, residual norm %g", 
                      it, ndof, residual_norm);

    // If maximum allowed residual norm is exceeded, fail.
    if (residual_norm > newton_max_allowed_residual_norm) {
      if (verbose) {
        info("Current residual norm: %g", residual_norm);
        info("Maximum allowed residual norm: %g", newton_max_allowed_residual_norm);
        info("Newton solve not successful, returning false.");
      }
      return false;
    }

    // If residual norm is within tolerance, or the maximum number 
    // of iteration has been reached, then quit.
    if ((residual_norm < newton_tol || it > newton_max_iter) && it > 1) break;

    // Solve the linear system.
    if(!solver->solve()) error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < num_stages*ndof; i++) 
      stage_vec[i] += newton_damping_coeff * solver->get_solution()[i];

    // Increase iteration counter.
    it++;
  }

  // If max number of iterations was exceeded, fail. 
  if (it >= newton_max_iter) {
    if (verbose) info("Maximum allowed number of Newton iterations exceeded, returning false.");
    return false;
  }

  // Calculate new coefficient vector using the stage vector and the Butcher's table.
  for (int i = 0; i < ndof; i++) {
    double increment = 0;
    for (int s = 0; s < num_stages; s++) {
      increment += time_step * bt->get_B(s) * stage_vec[s*ndof + i]; 
    }
    coeff_vec[i] += increment;
  } 

  // Clean up.
  delete matrix_left;
  delete matrix_right;
  delete vector_left;
  delete vector_right;
  delete solver;

  // Delete stage spaces, but not the first (original) one.
  for (int i = 1; i < num_stages; i++) delete stage_spaces[i];

  // Delete residuals.
  for (int i = 0; i < num_stages; i++) delete residuals[i];

  // TODO: Delete stage_wf, in particular its external solutions 
  // stage_time_sol[i], i = 0, 1, ..., num_stages-1.

  // Delete stage_vec and u_prev_vec.
  delete [] stage_vec;
  delete [] u_prev_vec;
  
  return true;
}
