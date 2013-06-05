#include "definitions.h"

using namespace Hermes::Hermes2D::RefinementSelectors;

typedef std::complex<double> complex;

//  This problem describes the distribution of the vector potential in
//  a 2D domain comprising a wire carrying electrical current, air, and
//  an iron which is not under voltage.
//
//  PDE: -(1/mu)Laplace A + ii*omega*gamma*A - J_ext = 0.
//
//  Domain: Rectangle of height 0.003 and width 0.004. Different
//  materials for the wire, air, and iron (see mesh file domain2.mesh).
//
//  BC: Zero Dirichlet on the top and right edges, zero Neumann
//  elsewhere.
//
//  The following parameters can be changed:
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.95;                    // This is a quantitative parameter of Adaptivity.

// Error calculation & adaptivity.
DefaultErrorCalculator<complex, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<complex> stoppingCriterion(0.9);
// Adaptivity processor class.
Adapt<complex> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-3;

// Problem parameters.
const double MU_0 = 4.0*M_PI*1e-7;
const double MU_IRON = 1e3 * MU_0;
const double GAMMA_IRON = 6e6;
const double J_EXT = 1e6;
const double FREQ = 5e3;
const double OMEGA = 2 * M_PI * FREQ;


int main(int argc, char* argv[])
{
  Hermes::Mixins::TimeMeasurable m;
  m.tick();

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

  // Initialize boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<complex> bc_essential("Dirichlet", complex(0.0, 0.0));
  EssentialBCs<complex> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<complex> space(new H1Space<complex>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();

  // Initialize the weak formulation.
  CustomWeakForm wf("Air", MU_0, "Iron", MU_IRON, GAMMA_IRON,
    "Wire", MU_0, complex(J_EXT, 0.0), OMEGA);

  // Initialize coarse and reference mesh solution.
  MeshFunctionSharedPtr<complex> sln(new Hermes::Hermes2D::Solution<complex>());
  MeshFunctionSharedPtr<complex> ref_sln(new Hermes::Hermes2D::Solution<complex>());

  // Initialize refinement selector.
  H1ProjBasedSelector<complex> selector(CAND_LIST);

  // Initialize views.
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 600, 350));
  Views::ScalarView sview2("Ref. Solution", new Views::WinGeom(0, 0, 600, 350));
  Views::OrderView oview("Polynomial orders", new Views::WinGeom(610, 0, 520, 350));

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  DiscreteProblem<complex> dp(&wf, space);

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::NewtonSolver<complex> newton(&dp);
  SimpleGraph graph_dof_est;
    
  Views::MeshView m1, m2;
  Views::OrderView o1, o2;
  // Adaptivity loop:
  int as = 1; bool done = false;
  adaptivity.set_space(space);
  do
  {
    // Construct globally refined reference mesh and setup reference space->
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<complex>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
    SpaceSharedPtr<complex> ref_space = ref_space_creator.create_ref_space();
    
    newton.set_space(ref_space);

    int ndof_ref = ref_space->get_num_dofs();

    // Initialize reference problem.

    // Initial coefficient vector for the Newton's method.
    complex* coeff_vec = new complex[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(complex));

    // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
    try
    {
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.print_msg();
    }

    Hermes::Hermes2D::Solution<complex>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<complex> ogProjection;
    ogProjection.project_global(space, ref_sln, sln);

    // View the coarse mesh solution and polynomial orders.
    MeshFunctionSharedPtr<double> real_filter(new RealFilter(sln));
    MeshFunctionSharedPtr<double> rreal_filter(new RealFilter(ref_sln));
    sview2.show(rreal_filter);

    oview.show(space);

    // Calculate element errors and total error estimate.
    errorCalculator.calculate_errors(sln, ref_sln);

    std::cout << (std::string)"Relative error: " << errorCalculator.get_total_error_squared() * 100. << '%' << std::endl;

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(space->get_num_dofs(), errorCalculator.get_total_error_squared() * 100.);
    sview.show(errorCalculator.get_errorMeshFunction());

    // If err_est too large, adapt the mesh->
    if(errorCalculator.get_total_error_squared()  * 100. < ERR_STOP)
      done = true;
    else
    {
      std::cout << (std::string)"Adapting..." << std::endl << std::endl;
      adaptivity.adapt(&selector);
    }

    // Clean up.
    delete [] coeff_vec;

    // Increase counter.
    as++;
  }
  while (done == false);

    graph_dof_est.save("conv_dof_est.dat");

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");

  MeshFunctionSharedPtr<double> real_filter(new RealFilter(ref_sln));
  sview.show(real_filter);

  m.tick();
  std::cout << m.accumulated();

  // Wait for all views to be closed.
  Views::View::wait();
  return 0;
}