#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This test makes sure that example 09-timedep works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

const int P_INIT = 3;                             // initial polynomial degree in elements
const int INIT_REF_NUM = 0;                       // number of initial uniform refinements
const double TAU = 200.0;                         // time step in seconds
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem constants
const double T_INIT = 10;                         // temperature of the ground (also initial temperature)
const double ALPHA = 10;                          // heat flux coefficient for Newton's boundary condition
const double LAMBDA = 1e5;                        // thermal conductivity of the material
const double HEATCAP = 1e6;                       // heat capacity
const double RHO = 3000;                          // material density
const double FINAL_TIME = 2100;                   // length of time interval (24 hours) in seconds

// Global variable.
double TIME = 0;

// Time-dependent exterior temperature.
double temp_ext(double t) {
  return T_INIT + 10. * sin(2*M_PI*t/FINAL_TIME);
}

// Boundary markers.
const int BDY_GROUND = 1, BDY_AIR = 2;

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, u, v) / TAU +
         LAMBDA * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, ext->fn[0], v) / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return LAMBDA * ALPHA * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return LAMBDA * ALPHA * temp_ext(TIME) * int_v<Real, Scalar>(n, wt, v);
}


int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("cathedral.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(2, 5);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_GROUND);
  bc_types.add_bc_newton(BDY_AIR);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_const(BDY_GROUND, T_INIT);

  // Initialize an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);
 
  // Initialize the solution.
  Solution tsln;

  // Set the initial condition.
  tsln.set_const(&mesh, T_INIT);

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, double>, bilinear_form<Ord, Ord>);
  wf.add_matrix_form_surf(bilinear_form_surf<double, double>, bilinear_form_surf<Ord, Ord>, BDY_AIR);
  wf.add_vector_form(linear_form<double, double>, linear_form<Ord, Ord>, HERMES_ANY, &tsln);
  wf.add_vector_form_surf(linear_form_surf<double, double>, linear_form_surf<Ord, Ord>, BDY_AIR);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Time stepping:
  int nsteps = (int)(FINAL_TIME/TAU + 0.5);
  bool rhs_only = false;
  for(int ts = 1; ts <= nsteps; ts++)
  {
    info("---- Time step %d, time %3.5f, ext_temp %g", ts, TIME, temp_ext(TIME));

    // First time assemble both the stiffness matrix and right-hand side vector,
    // then just the right-hand side vector.
    if (rhs_only == false) info("Assembling the stiffness matrix and right-hand side vector.");
    else info("Assembling the right-hand side vector (only).");
    dp.assemble(matrix, rhs, rhs_only);
    rhs_only = true;

    // Solve the linear system and if successful, obtain the solution.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution::vector_to_solution(solver->get_solution(), &space, &tsln);
    else 
      error ("Matrix solver failed.\n");

    // Update the time variable.
    TIME += TAU;
  }

  double sum = 0;
  for (int i=0; i < ndof; i++) sum += solver->get_solution()[i];
  printf("coefficient sum = %g\n", sum);

  // Actual test. The value of 'sum' depend on the
  // current shapeset. If you change the shapeset,
  // you need to correct this number.
  int success = 1;
  if (fabs(sum - 4737.73) > 1e-1) success = 0;

  if (success == 1) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

