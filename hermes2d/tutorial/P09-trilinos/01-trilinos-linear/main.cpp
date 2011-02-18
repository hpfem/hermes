#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Teuchos;

//  The purpose of this example is to show how to use Trilinos
//  for linear PDE problems. It compares performance of the LinearProblem 
//  class in Hermes using the UMFpack matrix solver with the performance
//  of the Trilinos NOX solver (using Newton's method or JFNK, with or 
//  without preconditioning).
//
//  PDE: Poisson equation.
//
//  Domain: Square (-1, 1)^2.
//
//  BC: Nonhomogeneous Dirichlet, see the function essential_bc_values() below.
//
//  Exact solution: x*x + y*y.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
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
                                                  // preconditioner from IFPACK (see solver/aztecoo.h)

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
  wf1.add_matrix_form(callback(bilinear_form), HERMES_SYM);
  wf1.add_vector_form(callback(linear_form));

  // Initialize the solution.
  Solution sln1;
  
  // Initialize the linear discrete problem.
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

  // View the solution and mesh.
  ScalarView sview("Solution", new WinGeom(0, 0, 440, 350));
  sview.show(&sln1);
  //OrderView  oview("Polynomial orders", new WinGeom(450, 0, 400, 350));
  //oview.show(&space);
  
  // TRILINOS PART:

  info("---- Assembling by DiscreteProblem, solving by NOX:");

  // Initialize the weak formulation for Trilinos.
  WeakForm wf2(1, JFNK ? true : false);
  wf2.add_matrix_form(callback(jacobian_form), HERMES_SYM);
  wf2.add_vector_form(callback(residual_form));
  
  // Initialize DiscreteProblem.
  is_linear = false;
  DiscreteProblem dp2(&wf2, &space, is_linear);
  
  // Time measurement.
  cpu_time.tick(HERMES_SKIP);

  // Set initial vector for NOX.
  info("Projecting to obtain initial vector for the Newton's method.");
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
  
  delete coeff_vec;

  // Choose preconditioning.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  if (PRECOND)
  {
    if (JFNK) nox_solver.set_precond(pc);
    else nox_solver.set_precond("ML");
  }

  // Assemble and solve using NOX.
  Solution sln2;
  if (nox_solver.solve())
  {
    Solution::vector_to_solution(nox_solver.get_solution(), &space, &sln2);

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      nox_solver.get_num_iters(), nox_solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      nox_solver.get_num_lin_iters(), nox_solver.get_achieved_tol());
  }
  else error("NOX failed");

  // CPU time needed by NOX.
  double time2 = cpu_time.tick().last();

  // Show the NOX solution.
  ScalarView view2("Solution 2", new WinGeom(450, 0, 440, 350));
  view2.show(&sln2);

  // Calculate errors.
  Solution ex;
  ex.set_exact(&mesh, &exact);
  double rel_err_1 = calc_rel_error(&sln1, &ex, HERMES_H1_NORM) * 100;
  info("Solution 1 (%s):  exact H1 error: %g (time %g s)", MatrixSolverNames[matrix_solver].c_str(), rel_err_1, time1);
  double rel_err_2 = calc_rel_error(&sln2, &ex, HERMES_H1_NORM) * 100;
  info("Solution 2 (NOX): exact H1 error: %g (time %g + %g = %g [s])", rel_err_2, proj_time, time2, proj_time+time2);
 
  // Wait for all views to be closed.
  View::wait();
  
  return 0;
}
