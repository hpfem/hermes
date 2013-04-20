#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"
#include "function/function.h"

using namespace RefinementSelectors;

//  This example uses the Picard's method to solve a nonlinear problem.
//  Try to run this example with PICARD_NUM_LAST_ITER_USED = 1 for 
//  comparison (Anderson acceleration turned off).)
//
//  PDE: Stationary heat transfer equation with nonlinear thermal
//  conductivity, -div[lambda(u) grad u] + src(x, y) = 0.
//
//  Nonlinearity: lambda(u) = 1 + Hermes::pow(u, 4).
//
//  Picard's linearization: -div[lambda(u^n) grad u^{n+1}] + src(x, y) = 0.
//
//  Domain: square (-10, 10)^2.
//
//  BC: Nonconstant Dirichlet.
//
//  The following parameters can be changed:

// Initial polynomial degree.
const int P_INIT = 2;                             
// Number of initial uniform mesh refinements.
const int INIT_GLOB_REF_NUM = 3;                  
// Number of initial refinements towards boundary.
const int INIT_BDY_REF_NUM = 5;                   
// Value for custom constant initial condition.
const double INIT_COND_CONST = 3.0;               
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Picard's method.
// Number of last iterations used. 
// 1... standard fixed point.
// >1... Anderson acceleration.
const int PICARD_NUM_LAST_ITER_USED = 4;          
// 0 <= beta <= 1... parameter for the Anderson acceleration. 
const double PICARD_ANDERSON_BETA = 0.2;          
// Stopping criterion for the Picard's method.
const double PICARD_TOL = 1e-3;
// Maximum allowed number of Picard iterations.
const int PICARD_MAX_ITER = 100;                  

// Problem parameters.
double heat_src = 1.0;
double alpha = 4.0;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh->refine_all_elements();
  mesh->refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  CustomEssentialBCNonConst bc_essential("Bdy");
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();

  // Initialize previous iteration solution for the Picard's method.
  MeshFunctionSharedPtr<double> sln_prev_iter(new ConstantSolution<double>(mesh, INIT_COND_CONST));

  // Initialize the weak formulation.
  CustomNonlinearity lambda(alpha);
  Hermes2DFunction<double> src(-heat_src);
  CustomWeakFormPicard wf(sln_prev_iter, &lambda, &src);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, space);
  dp.set_linear();

  // Initialize the Picard solver.
  PicardSolver<double> picard(&dp);
  picard.use_Anderson_acceleration(true);

  // Perform the Picard's iteration (Anderson acceleration on by default).
  picard.set_tolerance(PICARD_TOL);
  picard.set_max_allowed_iterations(PICARD_MAX_ITER);
  picard.set_num_last_vector_used(PICARD_NUM_LAST_ITER_USED);
  picard.set_anderson_beta(PICARD_ANDERSON_BETA);
  try
  {
    picard.solve(sln_prev_iter);
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  // Translate the coefficient vector into a Solution. 
  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  Solution<double>::vector_to_solution(picard.get_sln_vector(), space, sln);
  
  // Visualise the solution and mesh.
  ScalarView s_view("Solution", new WinGeom(0, 0, 440, 350));
  s_view.show_mesh(false);
  s_view.show(sln);
  OrderView o_view("Mesh", new WinGeom(450, 0, 420, 350));
  o_view.show(space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

