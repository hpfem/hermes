#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Teuchos;

// This test makes sure that example 40-trilinos-linear works correctly.

const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int P_INIT = 3;                             // Initial polynomial degree of all mesh elements.
const bool JFNK = false;                          // true = Jacobian-free method (for NOX),
                                                  // false = Newton (for NOX).
const bool PRECOND = true;                        // Preconditioning by jacobian in case of JFNK (for NOX),
                                                  // default ML preconditioner in case of Newton.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "least-squares";     // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers).
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  //  preconditioner from IFPACK (see solver/aztecoo.h)

// Boundary markers.
const int BDY_DIRICHLET = 1;

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(double x, double y)
{
  return x*x + y*y;
}

// Initial condition.
double init_cond(double x, double y, double &dx, double &dy)
{
	dx = 0;
	dy = 0;
	return 0;
}

// Exact solution.
double exact(double x, double y, double &dx, double &dy)
{
	dx = 2*x;
	dy = 2*y;
	return x*x +y*y;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char **argv)
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
 
  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_function(BDY_DIRICHLET, essential_bc_values);
 
  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof: %d", ndof);

  info("---- Assembling by DiscreteProblem, solving by %s:", MatrixSolverNames[matrix_solver].c_str());

  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Initialize weak formulation.
  WeakForm wf1;
  wf1.add_matrix_form(callback(bilinear_form));
  wf1.add_vector_form(callback(linear_form));

  // Initialize the solution.
  Solution sln1;
  
  // Initialize the linear FE problem.
  bool is_linear = true;
  DiscreteProblem dp1(&wf1, &space, is_linear);
    
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver)->set_solver(iterative_method);
    ((AztecOOSolver*) solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }
  
  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling the stiffness matrix and right-hand side vector.");
  dp1.assemble(matrix, rhs);
  
  // Solve the linear system and if successful, obtain the solution.
  info("Solving the matrix problem.");
  if(solver->solve())
    Solution::vector_to_solution(solver->get_solution(), &space, &sln1);
  else
    error ("Matrix solver failed.\n");

  delete(matrix);
  delete(rhs);
  delete(solver);
  
  // CPU time needed by UMFpack.
  double time1 = cpu_time.tick().last();
  
  // TRILINOS PART:
  info("---- Assembling by DiscreteProblem, solving by NOX:");

  // Initialize the weak formulation for Trilinos.
  WeakForm wf2(1, JFNK ? true : false);
  wf2.add_matrix_form(callback(jacobian_form), HERMES_SYM);
  wf2.add_vector_form(callback(residual_form));
  
  // Initialize DiscreteProblem.
  DiscreteProblem dp2(&wf2, &space);
  
  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Set initial vector for NOX.
  bool projected_ic;    

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the NOX solver.
  
  info("Projecting to obtain initial vector for the Newton's method.");
  projected_ic = true;
  scalar* coeff_vec = new scalar[ndof] ;
  Solution* init_sln = new ExactSolution(&mesh, init_cond);
  OGProjection::project_global(&space, init_sln, coeff_vec);
  delete init_sln;
  
  // Measure the projection time.
  double proj_time = cpu_time.tick().last();
  
  // Initialize the NOX solver with the vector "coeff_vec".
  info("Initializing NOX.");
  NoxSolver nox_solver(&dp2);
  nox_solver.set_init_sln(coeff_vec);
  
  if (!projected_ic)  delete  coeff_vec;
  else  delete [] coeff_vec;

  // Choose preconditioning.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  if (PRECOND)
  {
    if (JFNK) nox_solver.set_precond(pc);
    else nox_solver.set_precond("ML");
  }

  // Assemble and solve using NOX.
  bool solved = nox_solver.solve();

  Solution sln2;
  if (solved)
  {
    scalar *coeffs = nox_solver.get_solution();
    
    Solution::vector_to_solution(coeffs, &space, &sln2);
    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      nox_solver.get_num_iters(), nox_solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      nox_solver.get_num_lin_iters(), nox_solver.get_achieved_tol());
  }
  else error("NOX failed");

  // CPU time needed by NOX.
  double time2 = cpu_time.tick().last();

  // Calculate errors.
  Solution ex;
  ex.set_exact(&mesh, &exact);
  Adapt adaptivity(&space);
  bool solutions_for_adapt = false;
  double rel_err_1 = adaptivity.calc_err_exact(&sln1, &ex, solutions_for_adapt) * 100;
  info("Solution 1 (%s):  exact H1 error: %g (time %g s)", MatrixSolverNames[matrix_solver].c_str(), rel_err_1, time1);
  double rel_err_2 = adaptivity.calc_err_exact(&sln2, &ex, solutions_for_adapt) * 100;
  info("Solution 2 (NOX): exact H1 error: %g (time %g + %g = %g [s])", rel_err_2, proj_time, time2, proj_time+time2);

  info("Coordinate (-0.6, -0.6) hermes value = %lf", sln1.get_pt_value(-0.6, -0.6));
  info("Coordinate ( 0.4, -0.6) hermes value = %lf", sln1.get_pt_value( 0.4, -0.6));
  info("Coordinate ( 0.4,  0.4) hermes value = %lf", sln1.get_pt_value( 0.4,  0.4));
  info("Coordinate (-0.6,  0.0) hermes value = %lf", sln1.get_pt_value(-0.6,  0.0));
  info("Coordinate ( 0.0,  0.0) hermes value = %lf", sln1.get_pt_value( 0.0,  0.0));

  info("Coordinate (-0.6, -0.6) nox value = %lf", sln2.get_pt_value(-0.6, -0.6));
  info("Coordinate ( 0.4, -0.6) nox value = %lf", sln2.get_pt_value( 0.4, -0.6));
  info("Coordinate ( 0.4,  0.4) nox value = %lf", sln2.get_pt_value( 0.4,  0.4));
  info("Coordinate (-0.6,  0.0) nox value = %lf", sln2.get_pt_value(-0.6,  0.0));
  info("Coordinate ( 0.0,  0.0) nox value = %lf", sln2.get_pt_value( 0.0,  0.0));

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  int success = 1;
  double eps = 1e-5;
  if (fabs(sln2.get_pt_value(-0.6, -0.6) - 0.720000) > eps) {
    printf("Coordinate (-0.6, -0.6) nox value is %g\n", sln2.get_pt_value(-0.6, -0.6));
    success = 0;
  }
  if (fabs(sln2.get_pt_value( 0.4, -0.6) - 0.520000) > eps) {
    printf("Coordinate ( 0.4, -0.6) nox value is %g\n", sln2.get_pt_value( 0.4, -0.6));
    success = 0;
  }
  if (fabs(sln2.get_pt_value( 0.4,  0.4) - 0.320000) > eps) {
    printf("Coordinate ( 0.4,  0.4) nox value is %g\n", sln2.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if (fabs(sln2.get_pt_value(-0.6,  0.0) - 0.360000) > eps) {
    printf("Coordinate (-0.6,  0.0) nox value is %g\n", sln2.get_pt_value(-0.6,  0.0));
    success = 0;
  }
  if (fabs(sln2.get_pt_value( 0.0,  0.0) - 0.000000) > eps) {
    printf("Coordinate ( 0.0,  0.0) nox value is %g\n", sln2.get_pt_value( 0.0,  0.0));
    success = 0;
  }

  if (fabs(sln1.get_pt_value(-0.6, -0.6) - 0.720000) > eps) {
    printf("Coordinate (-0.6, -0.6) hermes value is %g\n", sln1.get_pt_value(-0.6, -0.6));
    success = 0;
  }
  if (fabs(sln1.get_pt_value( 0.4, -0.6) - 0.520000) > eps) {
    printf("Coordinate ( 0.4, -0.6) hermes value is %g\n", sln1.get_pt_value( 0.4, -0.6));
    success = 0;
  }
  if (fabs(sln1.get_pt_value( 0.4,  0.4) - 0.320000) > eps) {
    printf("Coordinate ( 0.4,  0.4) hermes value is %g\n", sln1.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if (fabs(sln1.get_pt_value(-0.6,  0.0) - 0.360000) > eps) {
    printf("Coordinate (-0.6,  0.0) hermes value is %g\n", sln1.get_pt_value(-0.6,  0.0));
    success = 0;
  }
  if (fabs(sln1.get_pt_value( 0.0,  0.0) - 0.000000) > eps) {
    printf("Coordinate ( 0.0,  0.0) hermes value is %g\n", sln1.get_pt_value( 0.0,  0.0));
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
