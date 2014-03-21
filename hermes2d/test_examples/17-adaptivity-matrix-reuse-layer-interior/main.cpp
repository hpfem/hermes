#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;


// Initial polynomial degree of mesh elements.
const int P_INIT = 2;                             
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;

// Problem parameters.
// Slope of the layer.
double slope = 60;                                

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
  
  // Define exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh, slope));

  // Define custom function f.
  CustomFunction f(slope);

  // Initialize the weak formulation.
  Hermes::Hermes1DFunction<double> lambda(1.0);
  WeakFormSharedPtr<double> wf(new DefaultWeakFormPoisson<double>(HERMES_ANY, &lambda, &f));
  
  // Initialize boundary conditions
  DefaultEssentialBCNonConst<double> bc_essential("Bdy", exact_sln);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  
  // Initialize approximate solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>());

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(H2D_HP_ANISO);

  DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(CalculatedErrorType::RelativeErrorToGlobalNorm, 1);
  AdaptStoppingCriterionSingleElement<double> criterion(0.);
  
  AdaptSolverCriterionFixed global_criterion(2);

  AdaptSolver<double, LinearSolver<double> > adaptSolver(space, wf, &errorCalculator, &criterion, &selector, &global_criterion);

  adaptSolver.switch_visualization(false);
  adaptSolver.set_verbose_output(true);
  adaptSolver.solve(hpAdaptivity);
  return 0;
}
