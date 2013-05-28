#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
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
const int INIT_REF = 2;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 1;
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.9;

// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-1;

int main(int argc, char* args[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF; i++)
    mesh->refine_all_elements();

  // Create an L2 space->
  SpaceSharedPtr<double> space(new L2Space<double>(mesh, P_INIT));

  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Display the mesh->
  OrderView oview("Coarse mesh", new WinGeom(0, 0, 440, 350));
  oview.show(space);

  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  MeshFunctionSharedPtr<double> ref_sln(new Solution<double>);

  // Initialize the weak formulation.
  CustomWeakForm wf("Bdy_bottom_left", mesh);

  ScalarView view1("Solution", new WinGeom(900, 0, 450, 350));
  view1.fix_scale_width(60);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);

  int as = 1; bool done = false;
  do
  {
    // Construct globally refined reference mesh
    // and setup reference space->
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
    SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();

    linear_solver.set_space(ref_space);

    // Solve the linear system. If successful, obtain the solution.
    try
    {
      linear_solver.solve();
      Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), ref_space, ref_sln);
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }
    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<double> ogProjection;
    ogProjection.project_global(space, ref_sln, sln, HERMES_L2_NORM);

    MeshFunctionSharedPtr<double> val_filter(new ValFilter(ref_sln, 0.0, 1.0));

    // View the coarse mesh solution.
    view1.show(val_filter);
    oview.show(space);

    // Calculate element errors and total error estimate.
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    adaptivity.set_space(space);

    std::cout << "Error: " << err_est_rel << "%." << std::endl;

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space<double>::get_num_dofs(space), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");

    // If err_est_rel too large, adapt the mesh->
    if(err_est_rel < ERR_STOP)
      done = true;
    else
      done = adaptivity.adapt(&selector);
    as++;
  }
  while (done == false);

  // Wait for keyboard or mouse input.
  View::wait();
  return 0;
}
