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
  double val = 0;
  for (int i=0; i<num; i++) {
    // Looking for an element in which lies the point x[i]:
    Iterator* I = new Iterator((Space*)user_data);
    Element *e;
    double f = 0;
    // Traversing all elements of the provided space.
    while ((e = I->next_active_element()) != NULL) {
      // Does x[i] lie within this element?
      if(x[i] <= e->x2 && x[i] > e->x1)
        // If so, calculate the value of the correct component of 
        // the solution on the provided space
        f = e->get_solution_value(x[i], OGProjection::sln);
    }
    val += weights[i] * f * v[i];
  }
  return val;
}

double OGProjection::H1_projection_liform(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  double val = 0;
  for (int i=0; i<num; i++) {
    // Looking for an element in which lies the point x[i]:
    Iterator* I = new Iterator((Space*)user_data);
    Element *e;
    double f = 0;
    double dfdx = 0;
    // Traversing all elements of the provided space.
    while ((e = I->next_active_element()) != NULL) {
      // Does x[i] lie within this element?
      if(x[i] <= e->x2 && x[i] > e->x1) {
        // If so, calculate the value and derivative of the
        // correct component of the solution on the provided space
        f = e->get_solution_value(x[i], OGProjection::sln);
        dfdx = e->get_solution_deriv(x[i], OGProjection::sln);
      }
    }
    val += weights[i] * (f * v[i] + dfdx * dvdx[i]);
  }
  return val;
}
void OGProjection::project_global(Space *space, Space* space_ref, MatrixSolverType matrix_solver, ProjNormType proj_norm)
{
  for(int i = 0; i < space->get_n_sln(); i++)
    OGProjection::project_global(i, space, space_ref, matrix_solver, proj_norm);
}

void OGProjection::project_global(int sln, Space *space, Space* space_ref, MatrixSolverType matrix_solver, ProjNormType proj_norm)
{
  WeakForm* wf = new WeakForm;
  if (proj_norm == HERMES_L2_NORM) {
      wf->add_matrix_form(OGProjection::L2_projection_biform);
      wf->add_vector_form(OGProjection::L2_projection_liform, space_ref);
  } else if (proj_norm == HERMES_H1_NORM) {
      wf->add_matrix_form(OGProjection::H1_projection_biform);
      wf->add_vector_form(OGProjection::H1_projection_liform, space_ref);
  } else
      throw std::runtime_error("Unknown proj_norm in project_global.");

  // Set the desired solution component.
  OGProjection::sln = sln;

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
  memset(coeff_vec, 0, ndof * sizeof(double));

  info("Assembling projection linear system. ndofs: %d", ndof);
  dp->assemble(matrix, rhs);
  
  // Solve the linear system.
  if(!solver->solve())
    error ("Matrix solver failed.\n");
  for (int i = 0; i < ndof; i++) coeff_vec[i] = solver->get_solution()[i];

  // Save the projection.
  vector_to_solution(coeff_vec, space, sln);

  // Cleanup.
  delete dp;
  delete matrix;
  delete rhs;
  delete solver;
}
