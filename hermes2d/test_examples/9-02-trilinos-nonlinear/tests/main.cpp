#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

using namespace Teuchos;

// This test makes sure that example 41-trilinos-nonlinear works correctly.

const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int P_INIT = 3;                             // Initial polynomial degree of all mesh elements.
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.

const bool JFNK = false;                          // true = jacobian-free method,
                                                  // false = Newton.
const int PRECOND = 2;                            // Preconditioning by jacobian (1) or approximation of jacobian (2)
                                                  // in case of JFNK,
                                                  // Default ML proconditioner in case of Newton.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers).
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "least-squares";     // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers).
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  //  preconditioner from IFPACK (see solver/aztecoo.h)
// NOX parameters.
unsigned message_type = NOX::Utils::Error | NOX::Utils::Warning | NOX::Utils::OuterIteration | NOX::Utils::InnerIteration | NOX::Utils::Parameters | NOX::Utils::LinearSolverDetails;
                                                  // NOX error messages, see NOX_Utils.h.

double ls_tolerance = 1e-5;                       // Tolerance for linear system.
unsigned flag_absresid = 0;                       // Flag for absolute value of the residuum.
double abs_resid = 1.0e-3;                        // Tolerance for absolute value of the residuum.
unsigned flag_relresid = 1;                       // Flag for relative value of the residuum.
double rel_resid = 1.0e-2;                        // Tolerance for relative value of the residuum.
int max_iters = 100;                              // Max number of iterations.

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("../square.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc);

  // Create an H1 space with default shapeset.
  H1Space<double> space1(&mesh, &bcs, P_INIT);
  H1Space<double> space2(&mesh, &bcs, P_INIT);
  int ndof = Space<double>::get_num_dofs(&space1);

  // Initialize weak formulation,
  CustomWeakForm wf1;

  // Initialize the discrete problem.
  DiscreteProblem<double> dp1(&wf1, &space1);

  // Initialize the solution.
  Solution<double> sln1;

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  double* coeff_vec = new double[ndof];
  // We can start with a zero vector.
  memset(coeff_vec, 0, ndof * sizeof(double));
  // Or we can project the initial condition to obtain the initial
  // coefficient vector.
  ////CustomInitialSolution sln_tmp(&mesh);
  //OGProjection::project_global(&space, &sln_tmp, coeff_vec);

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::Solution<double> sln;
  Hermes::Hermes2D::NewtonSolver<double> newton(&dp1);
  newton.set_verbose_output(true);
  try
  {
    newton.solve(coeff_vec);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.printMsg();
  }
  Hermes::Hermes2D::Solution<double>::vector_to_solution(newton.get_sln_vector(), &space1, &sln);

  // Translate the resulting coefficient vector into the Solution sln1.
  Solution<double>::vector_to_solution(coeff_vec, &space1, &sln1);

  // TRILINOS PART:

  // Project the initial condition to obtain the initial
  // coefficient vector.
  ZeroSolution<double> sln_tmp(&mesh);
  OGProjection<double> ogProjection;
  ogProjection.project_global(&space2, &sln_tmp, coeff_vec);

  // Initialize the weak formulation for Trilinos.
  CustomWeakForm wf2(JFNK, PRECOND == 1, PRECOND == 2);

  // Initialize DiscreteProblem.
  DiscreteProblem<double> dp2(&wf2, &space2);

  // Initialize the NOX solver with the vector "coeff_vec".
  NewtonSolverNOX<double> nox_solver(&dp2);
  nox_solver.set_output_flags(message_type);
  nox_solver.set_ls_tolerance(ls_tolerance);
  nox_solver.set_conv_rel_resid(rel_resid);
  nox_solver.set_conv_iters(max_iters);

  // Choose preconditioning.
  MlPrecond<double> pc("sa");
  if(PRECOND)
  {
    if(JFNK) nox_solver.set_precond(pc);
    else nox_solver.set_precond("ML");
  }

  // Solve the nonlinear problem using NOX.
  Solution<double> sln2;
  try{
    nox_solver.solve(coeff_vec);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.printMsg();
  }
  Solution<double>::vector_to_solution(nox_solver.get_sln_vector(), &space2, &sln2);

  // Calculate errors.
  CustomExactSolution ex(&mesh);
  double rel_err_1 = Global<double>::calc_rel_error(&sln1, &ex, HERMES_H1_NORM) * 100;
  double rel_err_2 = Global<double>::calc_rel_error(&sln2, &ex, HERMES_H1_NORM) * 100;

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  int success = 1;
  double eps = 1e-5;
  if(fabs(sln2.get_pt_value(0.6, 0.6) - 0.057600) > eps) {
    printf("Coordinate (0.6, 0.6) sln2 value is %g\n", sln2.get_pt_value(0.6, 0.6));
    success = 0;
  }
  if(fabs(sln2.get_pt_value( 0.4, 0.6) - 0.057600) > eps) {
    printf("Coordinate ( 0.4, 0.6) sln2 value is %g\n", sln2.get_pt_value( 0.4, 0.6));
    success = 0;
  }
  if(fabs(sln2.get_pt_value( 0.4,  0.4) - 0.057600) > eps) {
    printf("Coordinate ( 0.4,  0.4) sln2 value is %g\n", sln2.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if(fabs(sln2.get_pt_value(0.6,  0.0) - 0.000000) > eps) {
    printf("Coordinate (0.6,  0.0) sln2 value is %g\n", sln2.get_pt_value(0.6,  0.0));
    success = 0;
  }
  if(fabs(sln2.get_pt_value( 0.5,  0.5) - 0.062500) > eps) {
    printf("Coordinate ( 0.5,  0.5) sln2 value is %g\n", sln2.get_pt_value( 0.5,  0.5));
    success = 0;
  }

  if(fabs(sln1.get_pt_value(0.6, 0.6) - 0.057600) > eps) {
    printf("Coordinate (0.6, 0.6) sln1 value is %g\n", sln1.get_pt_value(0.6, 0.6));
    success = 0;
  }
  if(fabs(sln1.get_pt_value( 0.4, 0.6) - 0.057600) > eps) {
    printf("Coordinate ( 0.4, 0.6) sln1 value is %g\n", sln1.get_pt_value( 0.4, 0.6));
    success = 0;
  }
  if(fabs(sln1.get_pt_value( 0.4,  0.4) - 0.057600) > eps) {
    printf("Coordinate ( 0.4,  0.4) sln1 value is %g\n", sln1.get_pt_value( 0.4,  0.4));
    success = 0;
  }
  if(fabs(sln1.get_pt_value(0.6,  0.0) - 0.000000) > eps) {
    printf("Coordinate (0.6,  0.0) sln1 value is %g\n", sln1.get_pt_value(0.6,  0.0));
    success = 0;
  }
  if(fabs(sln1.get_pt_value( 0.5,  0.5) - 0.062500) > eps) {
    printf("Coordinate ( 0.5,  0.5) sln1 value is %g\n", sln1.get_pt_value( 0.5,  0.5));
    success = 0;
  }

  if(success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}