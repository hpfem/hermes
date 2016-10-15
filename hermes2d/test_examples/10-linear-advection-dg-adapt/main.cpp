#include "definitions.h"

//  This example solves a linear advection equation using Dicontinuous Galerkin (DG) method.
//  It is intended to show how evalutation of surface matrix forms that take basis functions defined
//  on different elements work. It is the same example as linear-advection-dg, but with automatic adaptivity.
//
//  PDE: \nabla \cdot (\Beta u) = 0, where \Beta = (-x_2, x_1) / |x| represents a circular counterclockwise flow field.
//
//  Domain: Square (0, 1) x (0, 1).
//
//  BC:    Dirichlet, u = 1 where \Beta(x) \cdot n(x) < 0, that is on[0,0.5] x {0}, and g = 0 anywhere else.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF = 1;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
int P_INIT = 1;
// Use Taylor shapeset - which does not have order > 2 implemented.
// This switches to h-adaptivity & turns on Vertex-based limiting.
bool USE_TAYLOR_SHAPESET = true;

// Error calculation & adaptivity.
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.75;
// Error calculator
DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = USE_TAYLOR_SHAPESET ? H2D_H_ANISO : H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-2;

int main(int argc, char* args[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i = 0; i < INIT_REF; i++)
    mesh->refine_all_elements();

  // Create an L2 space.
  // This is both coarse and fine for Taylor shapeset, coarse for standard Legendre.
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, USE_TAYLOR_SHAPESET ? std::max(P_INIT, 2) : P_INIT, (USE_TAYLOR_SHAPESET ? (Shapeset*)(new L2ShapesetTaylor) : (Shapeset*)(new L2ShapesetLegendre))));
  // Set this space to the AMR processor.
  adaptivity.set_space(space);

  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST);
  selector.set_error_weights(1., 1., 1.);

  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  MeshFunctionSharedPtr<double> refsln(new Solution<double>);

  // Initialize the weak formulation.
  WeakFormSharedPtr<double> wf(new CustomWeakForm("Bdy_bottom_left", mesh));
  ScalarView view1("Solution", new WinGeom(900, 0, 450, 350));
  view1.fix_scale_width(60);
  view1.set_min_max_range(0., 1.);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver;
  linear_solver.set_weak_formulation(wf);

  int adaptivity_step = 1; bool done = false;
  do
  {
    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<double>::ReferenceSpaceCreator refspace_creator(space, ref_mesh, USE_TAYLOR_SHAPESET ? 0 : 1);
    SpaceSharedPtr<double> refspace = refspace_creator.create_ref_space();

    // Solve the problem on the reference space.
    linear_solver.set_space(refspace);
    linear_solver.solve();

    // Get the Hermes2D Solution object from the solution vector.
    Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), refspace, refsln);

    // Project the reference solution to the coarse space in L2 norm for error calculation.
    OGProjection<double>::project_global(space, refsln, sln, HERMES_L2_NORM);

    // Calculate element errors and total error estimate.
    errorCalculator.calculate_errors(sln, refsln);
    double total_error_estimate = errorCalculator.get_total_error_squared() * 100;

	view1.show(refsln);

    std::cout << "Elements: " << ref_mesh->get_num_active_elements() << ", Error: " << total_error_estimate << "%." << std::endl;

    // If error is too large, adapt the (coarse) space.
    if (total_error_estimate > ERR_STOP)
      adaptivity.adapt(&selector);
    else
      done = true;

    adaptivity_step++;
  } while (done == false);

  // Wait for keyboard or mouse input.
  View::wait();
  return 0;
}