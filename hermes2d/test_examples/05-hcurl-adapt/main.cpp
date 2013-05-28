#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::RefinementSelectors;

typedef std::complex<double> complex;

//  This example comes with an exact solution, and it describes the diffraction
//  of an electromagnetic wave from a re-entrant corner. Convergence graphs saved
//  (both exact error and error estimate, and both wrt. dof number and cpu time).
//
//  PDE: time-harmonic Maxwell's equations
//
//  Known exact solution, see functions exact_sol_val(), exact_sol(), exact()
//
//  Domain: L-shape domain
//
//  Meshes: you can use either "lshape3q.mesh" (quadrilateral mesh) or
//          "lshape3t.mesh" (triangular mesh). See the mesh->load(...) command below.
//
//  BC: perfect conductor on boundary markers 1 and 6 (essential BC)
//      impedance boundary condition on rest of boundary (natural BC)
//
//  The following parameters can be changed:

// Set to "false" to suppress Hermes OpenGL visualization.
const bool HERMES_VISUALIZATION = true;
// Initial polynomial degree. NOTE: The meaning is different from
// standard continuous elements in the space H1. Here, P_INIT refers
// to the maximum poly order of the tangential component, and polynomials
// of degree P_INIT + 1 are present in element interiors. P_INIT = 0
// is for Whitney elements.
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;

// Error calculation & adaptivity.
DefaultErrorCalculator<complex, HERMES_HCURL_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<complex> stoppingCriterion(0.75);
// Adaptivity processor class.
Adapt<complex> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-3;

// Problem parameters.
const double MU_R   = 1.0;
const double KAPPA  = 1.0;
const double LAMBDA = 1.0;

// Bessel functions, exact solution, and weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("lshape3q.mesh", mesh);    // quadrilaterals
  //mloader.load("lshape3t.mesh", mesh);  // triangles

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++)  mesh->refine_all_elements();

  // Initialize boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<complex > bc_essential(Hermes::vector<std::string>("Corner_horizontal",
    "Corner_vertical"), 0);
  EssentialBCs<complex > bcs(&bc_essential);

  // Create an Hcurl space with default shapeset.
  SpaceSharedPtr<complex > space(new HcurlSpace<complex >(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();

  // Initialize the weak formulation.
  CustomWeakForm wf(MU_R, KAPPA);

  // Initialize coarse and reference mesh solutions.
  MeshFunctionSharedPtr<complex > sln(new Hermes::Hermes2D::Solution<complex >());
  MeshFunctionSharedPtr<complex > ref_sln(new Hermes::Hermes2D::Solution<complex >());

  // Initialize exact solution.
  MeshFunctionSharedPtr<complex > sln_exact(new CustomExactSolution(mesh));

  // Initialize refinement selector.
  HcurlProjBasedSelector<complex > selector(CAND_LIST);

  // Initialize views.
  Views::VectorView v_view("Solution (magnitude)", new Views::WinGeom(0, 0, 460, 350));
  v_view.set_min_max_range(0, 1.5);
  Views::OrderView  o_view("Polynomial orders", new Views::WinGeom(470, 0, 400, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est,
    graph_dof_exact, graph_cpu_exact;

  DiscreteProblem<complex > dp(&wf, space);

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::NewtonSolver<complex > newton(&dp);

  Views::Linearizer lin;
  Views::Orderizer ord;
  Views::Vectorizer vec;

  // Adaptivity loop:
  int as = 1; bool done = false;
  adaptivity.set_space(space);
  do
  {
    // Construct globally refined reference mesh and setup reference space->
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<complex >::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
    SpaceSharedPtr<complex > ref_space = ref_space_creator.create_ref_space();

    newton.set_space(ref_space);
    int ndof_ref = ref_space->get_num_dofs();

    // Initial coefficient vector for the Newton's method.
    complex* coeff_vec = new complex[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(complex));
    
    try
    {
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.print_msg();
    }
    Hermes::Hermes2D::Solution<complex >::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<complex > ogProjection;
    ogProjection.project_global(space, ref_sln, sln);

    // View the coarse mesh solution and polynomial orders.
    if(HERMES_VISUALIZATION)
    {
      MeshFunctionSharedPtr<double> real_filter(new RealFilter(sln));
      v_view.show(real_filter);
      o_view.show(space);
      lin.save_solution_vtk(real_filter, "sln.vtk", "a");
      ord.save_mesh_vtk(space, "mesh.vtk");
      lin.free();
    }

    // Calculate element errors and total error estimate.
    errorCalculator.calculate_errors(sln, sln_exact, false);
    double err_exact_rel = errorCalculator.get_total_error_squared() * 100;
    // Calculate exact error.
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    Hermes::Mixins::Loggable::Static::info("\nError estimate: %f%%.\n", err_est_rel);
    Hermes::Mixins::Loggable::Static::info("\nError exact: %f%%.\n", err_exact_rel);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(space->get_num_dofs(), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_dof_exact.add_values(space->get_num_dofs(), err_exact_rel);
    graph_dof_exact.save("conv_dof_exact.dat");

    // If err_est_rel too large, adapt the mesh->
    if(err_est_rel < ERR_STOP)
      done = true;
    else
    {
      done = adaptivity.adapt(&selector);
      // Increase the counter of performed adaptivity steps.
      if(done == false)  as++;
    }

    // Clean up.
    delete [] coeff_vec;
  }
  while (done == false);

  // Show the reference solution - the final result.
  if(HERMES_VISUALIZATION)
  {
    v_view.set_title("Fine mesh solution (magnitude)");
    MeshFunctionSharedPtr<double> real_filter(new RealFilter(ref_sln));
    v_view.show(real_filter);

    // Wait for all views to be closed.
    Views::View::wait();
  }
  return 0;
}