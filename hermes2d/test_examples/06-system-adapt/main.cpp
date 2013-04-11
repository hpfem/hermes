#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// This example explains how to use the multimesh adaptive hp-FEM,
// where different physical fields (or solution components) can be
// approximated using different meshes and equipped with mutually
// independent adaptivity mechanisms. For the tutorial purposes,
// we manufactured an exact solution for a simplified version of
// the FitzHugh-Nagumo equation. This equation, in its full form,
// is a prominent example of activator-inhibitor systems in two-component
// reaction-diffusion equations, It describes a prototype of an
// excitable system (e.g., a neuron).
//
// PDE: Linearized FitzHugh-Nagumo equation
//      -d_u^2 \Delta u - f(u) + \sigma v - g1 = 0,
//      -d_v^2 \Delta v - u + v - g2 = 0.
// In the original equation, f(u) = \lambda u - u^3 - \kappa. For
// simplicity, here we just take f(u) = u.
//
// Domain: Square (-1,1)^2.
//
// BC: Both solution components are zero on the boundary.
//
// Exact solution: The functions g1 and g2 were calculated so that
//                 the exact solution is:
//        u(x,y) = U(x)*U(y) where U(t) = Hermes::cos(M_PI*t/2)
//        v(x,y) = V(x)V(y) where V(t) = 1 - (exp(K*t)+exp(-K*t))/(exp(K) + exp(-K))
// Note: V(t) is the exact solution of the 1D singularly perturbed equation
//       -u'' + K*K*u = K*K in (-1, 1) with zero Dirichlet BC.
//
// The following parameters can be changed: In particular, compare hp- and
// h-adaptivity via the CAND_LIST option, and compare the multi-mesh vs.
// single-mesh using the MULTI parameter.

// Initial polynomial degree for u.
const int P_INIT_U = 2;
// Initial polynomial degree for v.
const int P_INIT_V = 1;
// Number of initial boundary refinements
const int INIT_REF_BDY = 5;
// MULTI = true  ... use multi-mesh,
// MULTI = false ... use single-mesh->
// Note: In the single mesh option, the meshes are
// forced to be geometrically the same but the
// polynomial degrees can still vary.
const bool MULTI = true;
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 1.3;
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
// More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int STRATEGY = 0;
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;
// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0.
const double CONV_EXP = 1;
// Stopping criterion for adaptivity.
const double ERR_STOP = 15.0;
// Adaptivity process stops when the number of degrees of freedom grows over
// this limit. This is mainly to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;
// Newton's method.
double NEWTON_TOL_FINE = 1e-0;
int max_allowed_iterations = 10;
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;

// Problem parameters.
const double D_u = 1;
const double D_v = 1;
const double SIGMA = 1;
const double LAMBDA = 1;
const double KAPPA = 1;
const double K = 100.;

int main(int argc, char* argv[])
{
  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;
  cpu_time.tick();

  // Load the mesh.
  MeshSharedPtr u_mesh(new Mesh), v_mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", u_mesh);
  if (MULTI == false) u_mesh->refine_towards_boundary("Bdy", INIT_REF_BDY);

  // Create initial mesh (master mesh).
  v_mesh->copy(u_mesh);

  // Initial mesh refinements in the v_mesh towards the boundary.
  if (MULTI == true) v_mesh->refine_towards_boundary("Bdy", INIT_REF_BDY);

  // Set exact solutions.
  MeshFunctionSharedPtr<double> exact_u(new ExactSolutionFitzHughNagumo1 (u_mesh));
  MeshFunctionSharedPtr<double> exact_v(new ExactSolutionFitzHughNagumo2 (v_mesh, K));

  // Define right-hand sides.
  CustomRightHandSide1 g1(K, D_u, SIGMA);
  CustomRightHandSide2 g2(K, D_v);

  // Initialize the weak formulation.
  CustomWeakForm wf(&g1, &g2);

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_u("Bdy", 0.0);
  EssentialBCs<double> bcs_u(&bc_u);
  DefaultEssentialBCConst<double> bc_v("Bdy", 0.0);
  EssentialBCs<double> bcs_v(&bc_v);

  // Create H1 spaces with default shapeset for both displacement components.
  SpaceSharedPtr<double> u_space(new H1Space<double>(u_mesh, &bcs_u, P_INIT_U));
  SpaceSharedPtr<double> v_space(new H1Space<double>(MULTI ? v_mesh : u_mesh, &bcs_v, P_INIT_V));

  // Initialize coarse and reference mesh solutions.
  MeshFunctionSharedPtr<double> u_sln(new Solution<double>()), v_sln(new Solution<double>()), u_ref_sln(new Solution<double>()), v_ref_sln(new Solution<double>());

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  Views::ScalarView s_view_0("Solution[0]", new Views::WinGeom(0, 0, 440, 350));
  s_view_0.show_mesh(false);
  Views::OrderView  o_view_0("Mesh[0]", new Views::WinGeom(450, 0, 420, 350));
  Views::ScalarView s_view_1("Solution[1]", new Views::WinGeom(880, 0, 440, 350));
  s_view_1.show_mesh(false);
  Views::OrderView o_view_1("Mesh[1]", new Views::WinGeom(1330, 0, 420, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;
  SimpleGraph graph_dof_exact, graph_cpu_exact;

  NewtonSolver<double> newton;

  newton.set_weak_formulation(&wf);

  newton.set_tolerance(1e-1);

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space->
    Mesh::ReferenceMeshCreator u_ref_mesh_creator(u_mesh);
    MeshSharedPtr u_ref_mesh = u_ref_mesh_creator.create_ref_mesh();
    Mesh::ReferenceMeshCreator v_ref_mesh_creator(v_mesh);
    MeshSharedPtr v_ref_mesh = v_ref_mesh_creator.create_ref_mesh();
    Space<double>::ReferenceSpaceCreator u_ref_space_creator(u_space, u_ref_mesh);
    SpaceSharedPtr<double> u_ref_space = u_ref_space_creator.create_ref_space();
    Space<double>::ReferenceSpaceCreator v_ref_space_creator(v_space, v_ref_mesh);
    SpaceSharedPtr<double> v_ref_space = v_ref_space_creator.create_ref_space();

    Hermes::vector<SpaceSharedPtr<double> > ref_spaces_const(u_ref_space, v_ref_space);

    newton.set_spaces(ref_spaces_const);

    int ndof_ref = Space<double>::get_num_dofs(ref_spaces_const);

    // Initialize reference problem.
    Hermes::Mixins::Loggable::Static::info("Solving on reference mesh.");

    // Time measurement.
    cpu_time.tick();

    // Perform Newton's iteration.
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      std::cout << e.what();
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), ref_spaces_const, Hermes::vector<MeshFunctionSharedPtr<double> >(u_ref_sln, v_ref_sln));

    // Project the fine mesh solution onto the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Projecting reference solution on coarse mesh.");
    OGProjection<double> ogProjection; ogProjection.project_global(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space),
                                                                   Hermes::vector<MeshFunctionSharedPtr<double> >(u_ref_sln, v_ref_sln),
                                                                   Hermes::vector<MeshFunctionSharedPtr<double> >(u_sln, v_sln));

    cpu_time.tick();

    // View the coarse mesh solution and polynomial orders.
    s_view_0.show(u_sln);
    o_view_0.show(u_space);
    s_view_1.show(v_sln);
    o_view_1.show(v_space);

    // Calculate element errors.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate and exact error.");
    Adapt<double>* adaptivity = new Adapt<double>(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space));

    // Calculate error estimate for each solution component and the total error estimate.
    Hermes::vector<double> err_est_rel;
    double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<MeshFunctionSharedPtr<double> >(u_sln, v_sln),
                                                        Hermes::vector<MeshFunctionSharedPtr<double> >(u_ref_sln, v_ref_sln),
                                                        &err_est_rel) * 100;

    // Calculate exact error for each solution component and the total exact error.
    Hermes::vector<double> err_exact_rel;
    bool solutions_for_adapt = false;
    double err_exact_rel_total = adaptivity->calc_err_exact(Hermes::vector<MeshFunctionSharedPtr<double> >(u_sln, v_sln),
                                                            Hermes::vector<MeshFunctionSharedPtr<double> >(exact_u, exact_v),
                                                            &err_exact_rel, solutions_for_adapt) * 100;

    // Time measurement.
    cpu_time.tick();

    // Report results.
    Hermes::Mixins::Loggable::Static::info("ndof_coarse[0]: %d, ndof_fine[0]: %d",
                                           u_space->get_num_dofs(), u_ref_space->get_num_dofs());
    Hermes::Mixins::Loggable::Static::info("err_est_rel[0]: %g%%, err_exact_rel[0]: %g%%", err_est_rel[0]*100, err_exact_rel[0]*100);
    Hermes::Mixins::Loggable::Static::info("ndof_coarse[1]: %d, ndof_fine[1]: %d",
                                           v_space->get_num_dofs(), v_ref_space->get_num_dofs());
    Hermes::Mixins::Loggable::Static::info("err_est_rel[1]: %g%%, err_exact_rel[1]: %g%%", err_est_rel[1]*100, err_exact_rel[1]*100);
    Hermes::Mixins::Loggable::Static::info("ndof_coarse_total: %d, ndof_fine_total: %d",
                                           Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space)),
                                           Space<double>::get_num_dofs(ref_spaces_const));
    Hermes::Mixins::Loggable::Static::info("err_est_rel_total: %g%%, err_est_exact_total: %g%%", err_est_rel_total, err_exact_rel_total);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space)),
                             err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");

    graph_dof_exact.add_values(Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space)),
                               err_exact_rel_total);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel_total);
    graph_cpu_exact.save("conv_cpu_exact.dat");

    // If err_est too large, adapt the mesh->
    if (err_est_rel_total < ERR_STOP)
      done = true;
    else
    {
      Hermes::Mixins::Loggable::Static::info("Adapting coarse mesh.");
      done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector<double> *>(&selector, &selector),
                               THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(u_space, v_space)) >= NDOF_STOP) done = true;

    // Clean up.
    delete adaptivity;

    // Increase counter.
    as++;
  }
  while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  // Views::View::wait();
  return 0;
}
