#pragma region notInterestingForThisPresentation
#include "definitions.h"
using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// Parameters
const double omega_c = 3.0 * M_PI / 2.0;
const double x_w = 0.0;
const double y_w = -3.0 / 4.0;
const double r_0 = 3.0 / 4.0;
const double alpha_w = 200.0;
const double x_p = -Hermes::sqrt(5.0) / 4.0;
const double y_p = -1.0 / 4.0;
const double alpha_p = 1000.0;
const double epsilon = 1.0 / 100.0;
#pragma endregion

// Initial polynomial degree of mesh elements.
const int P_INIT = 1;        

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;

// Adaptivity type to use.
AdaptivityType adaptivityType = hpAdaptivity;

// Global adaptivity criterion - either fixed (number of adaptive steps) or 
AdaptSolverCriterionErrorThreshold global_criterion(1.);

// Error calculator.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(CalculatedErrorType::RelativeErrorToGlobalNorm, 1);

// Criterion for stopping to refine elements within one adaptive step.
AdaptStoppingCriterionCumulative<double> criterion(0.45);

// Refinement selector.
H1ProjBasedSelector<double> selector(adaptivityType == pAdaptivity ? H2D_P_ANISO : (adaptivityType == hAdaptivity ? H2D_H_ANISO : H2D_HP_ANISO));

////// The main function ///////
int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  // Quadrilaterals.
  mloader.load("square_quad.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();
  
  // Set exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh, alpha_w, alpha_p, x_w, y_w, r_0, omega_c, epsilon, x_p, y_p));

  // Define right-hand side.
  CustomRightHandSide f(alpha_w, alpha_p, x_w, y_w, r_0, omega_c, epsilon, x_p, y_p);

  // Initialize the weak formulation.
  Hermes1DFunction<double> lambda(1.0);
  WeakFormSharedPtr<double> wf(new WeakFormsH1::DefaultWeakFormPoisson<double>(HERMES_ANY, &lambda, &f));
  
  // Initialize boundary conditions
  DefaultEssentialBCNonConst<double> bc_essential("Bdy", exact_sln);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  
  // Initialize approximate solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>());

  AdaptSolver<double, NewtonSolver<double> >::view_size = 600;
  AdaptSolver<double, NewtonSolver<double> >::scalar_views_switch = true;
  AdaptSolver<double, NewtonSolver<double> >::order_views_switch = true;
  AdaptSolver<double, NewtonSolver<double> >::base_views_switch = false;
  AdaptSolver<double, NewtonSolver<double> > adaptSolver(space, wf, &errorCalculator, &criterion, &selector, &global_criterion);
  adaptSolver.set_exact_solution(exact_sln);

  adaptSolver.switch_visualization(true);
  adaptSolver.set_verbose_output(true);
  adaptSolver.solve(adaptivityType);
  return 0;
}