// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Distributed under the terms of the BSD license (see the LICENSE
// file for the exact terms).
// Email: hermes1d@googlegroups.com, home page: http://hpfem.org/
#define HERMES_REPORT_INFO
#include "ogprojection.h"
#include "math.h"
#include "discrete_problem.h"
#include "weakform.h"

int OGProjection::sln = 0;

Space* OGProjection::ref_space = NULL;

ExactFunction OGProjection::fn = OGProjection::ref_mesh_fn;

double OGProjection::L2_projection_biform(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for (int i=0; i<num; i++) {
      val += weights[i] * (u[i] * v[i]);
  }
  return val;
}

double OGProjection::H1_projection_biform(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
  double val = 0;
  for (int i=0; i<num; i++) {
      val += weights[i] * (u[i] * v[i] + dudx[i] * dvdx[i]);
  }
  return val;
}

double OGProjection::L2_projection_liform(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  // Value of the projected function at gauss points:
  double* f = new double[num];
  fn(num, x, f, NULL, NULL);
  double val = 0;
  for (int i=0; i<num; i++) {
    val += weights[i] * (f[i] * v[i]);
  }
  delete [] f;
  return val;
}

double OGProjection::L2_projection_liform_0(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  OGProjection::sln = 0;
  return L2_projection_liform(num, x, weights, u_prev, du_prevdx, v, dvdx, user_data);
}

double OGProjection::L2_projection_liform_1(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  OGProjection::sln = 1;
  return L2_projection_liform(num, x, weights, u_prev, du_prevdx, v, dvdx, user_data);
}

double OGProjection::L2_projection_liform_2(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  OGProjection::sln = 2;
  return L2_projection_liform(num, x, weights, u_prev, du_prevdx, v, dvdx, user_data);
}

double OGProjection::L2_projection_liform_3(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  OGProjection::sln = 3;
  return L2_projection_liform(num, x, weights, u_prev, du_prevdx, v, dvdx, user_data);
}

double OGProjection::L2_projection_liform_4(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  OGProjection::sln = 4;
  return L2_projection_liform(num, x, weights, u_prev, du_prevdx, v, dvdx, user_data);
}



double OGProjection::H1_projection_liform(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  // Value of the projected function at gauss points:
  double* f = new double[num];
  double* dfdx = new double[num];
  fn(num, x, f, dfdx, NULL);
  double val = 0;
  for (int i=0; i<num; i++) {
      val += weights[i] * (f[i] * v[i] + dfdx[i] * dvdx[i]);
  }
  delete [] f;
  delete [] dfdx;
  return val;
}

double OGProjection::H1_projection_liform_0(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  OGProjection::sln = 0;
  return H1_projection_liform(num, x, weights, u_prev, du_prevdx, v, dvdx, user_data);
}

double OGProjection::H1_projection_liform_1(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  OGProjection::sln = 1;
  return H1_projection_liform(num, x, weights, u_prev, du_prevdx, v, dvdx, user_data);
}

double OGProjection::H1_projection_liform_2(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  OGProjection::sln = 2;
  return H1_projection_liform(num, x, weights, u_prev, du_prevdx, v, dvdx, user_data);
}

double OGProjection::H1_projection_liform_3(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  OGProjection::sln = 3;
  return H1_projection_liform(num, x, weights, u_prev, du_prevdx, v, dvdx, user_data);
}

double OGProjection::H1_projection_liform_4(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  OGProjection::sln = 4;
  return H1_projection_liform(num, x, weights, u_prev, du_prevdx, v, dvdx, user_data);
}

void OGProjection::project_internal(Space *space, MatrixSolverType matrix_solver, ProjNormType proj_norm, int sln_to_save)
{
  // Check.
  if(sln_to_save >= space->get_n_sln())
    error("The variable sln_to_save set incorrectly in OGProjection::project_internal.");

  // Setting of weakforms.
  WeakForm* wf = new WeakForm(space->get_n_sln());
  for(int i = 0; i < space->get_n_sln(); i++)
  {
    if (proj_norm == HERMES_L2_NORM) {
        wf->add_matrix_form(i, i, OGProjection::L2_projection_biform);
        switch ( i ) {
          case 0 : 
            wf->add_vector_form(0, OGProjection::L2_projection_liform_0);
            break;
          case 1 : 
            wf->add_vector_form(0, OGProjection::L2_projection_liform_1);
            break;
          case 2 : 
            wf->add_vector_form(0, OGProjection::L2_projection_liform_2);
            break;
          case 3 : 
            wf->add_vector_form(0, OGProjection::L2_projection_liform_3);
            break;
          case 4 : 
            wf->add_vector_form(0, OGProjection::L2_projection_liform_4);
            break;
        }
    } else if (proj_norm == HERMES_H1_NORM) {
        wf->add_matrix_form(i, i, OGProjection::H1_projection_biform);
        switch ( i ) {
          case 0 : 
            wf->add_vector_form(0, OGProjection::H1_projection_liform_0);
            break;
          case 1 : 
            wf->add_vector_form(1, OGProjection::H1_projection_liform_1);
            break;
          case 2 : 
            wf->add_vector_form(2, OGProjection::H1_projection_liform_2);
            break;
          case 3 : 
            wf->add_vector_form(3, OGProjection::H1_projection_liform_3);
            break;
          case 4 : 
            wf->add_vector_form(4, OGProjection::H1_projection_liform_4);
            break;
        }
    } else
        throw std::runtime_error("Unknown proj_norm in project_global.");
  }

  // Initialize the FE problem.
  DiscreteProblem* dp = new DiscreteProblem(wf, space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Get the target space's number of degrees of freedom.
  int ndof = Space::get_num_dofs(space);
  
  // Initialize coefficient vector to copy the solution to ref_space.
  double* coeff_vec = new double[ndof];
  set_zero(coeff_vec, ndof);

  info("Assembling projection linear system. ndofs: %d", ndof);
  dp->assemble(coeff_vec, matrix, rhs);
  
  // Solve the linear system.
  if(!solver->solve())
    error ("Matrix solver failed.\n");
  for (int i = 0; i < ndof; i++) coeff_vec[i] = solver->get_solution()[i];

  // Save the projection.
  if(sln_to_save == -1)
    for(int i = 0; i < space->get_n_sln(); i++)
      set_coeff_vector(coeff_vec + (ndof/(OGProjection::sln+1)) * i, space, i);
  else
      set_coeff_vector(coeff_vec + (ndof/(OGProjection::sln+1)) * sln_to_save, space, sln_to_save);

  // Cleanup.
  delete dp;
  delete matrix;
  delete rhs;
  delete solver;
  OGProjection::ref_space = NULL;
}

void OGProjection::project_global(Space *space, Space* space_ref, MatrixSolverType matrix_solver, ProjNormType proj_norm)
{
  // Check.
  if(space->get_n_sln() != space_ref->get_n_sln())
    error("Number of solutions of reference and coarse spaces differ in OGProjection::project_global.");

  // Settings specific for this type of projection.
  OGProjection::ref_space = space_ref;
  OGProjection::fn = OGProjection::ref_mesh_fn;

  OGProjection::project_internal(space, matrix_solver, proj_norm);

  return;
}

void OGProjection::project_global(Space *space, ExactFunction fun, int sln_to_save, MatrixSolverType matrix_solver, ProjNormType proj_norm)
{
  OGProjection::fn = fun;

  OGProjection::project_internal(space, matrix_solver, proj_norm, sln_to_save);

  return;
}

int OGProjection::ref_mesh_fn(int n, double x[], double f[], double dfdx[],
        void *data)
{
  // Check.
  if(OGProjection::ref_space == NULL)
    error("OGProjection::ref_space is NULL.");
  for (int i = 0; i < n; i++) {
    // Looking for an element in which lies the point x[i]:
    Iterator* I = new Iterator(OGProjection::ref_space);
    Element *e;
    // Traversing all elements of the provided space.
    while ((e = I->next_active_element()) != NULL) {
      // Does x[i] lie within this element?
      if(x[i] <= e->x2 && x[i] > e->x1) {
        // If so, calculate the value and derivative of the
        // correct component of the solution on the provided space
        f[i] = e->get_solution_value(x[i], OGProjection::sln);
        if(dfdx != NULL)
          dfdx[i] = e->get_solution_deriv(x[i], OGProjection::sln);
      }
    } 
  }
  return 0;
}
