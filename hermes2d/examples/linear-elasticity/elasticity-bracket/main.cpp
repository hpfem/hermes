#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example uses adaptive multimesh hp-FEM to solve a simple problem
// of linear elasticity. Note that since both displacement components
// have similar qualitative behavior, the advantage of the multimesh 
// discretization is less striking than for example in the tutorial 
// example P04-linear-adapt/02-system-adapt.
//
// PDE: Lame equations of linear elasticity, treated as a coupled system
//      of two PDEs.
//
// BC: u_1 = u_2 = 0 on Gamma_1 (right edge)
//     du_1/dn = f0 on Gamma_2 (top edge)
//     du_2/dn = f1 on Gamma_2 (top edge)
//     du_1/dn = du_2/dn = 0 elsewhere
//
// The following parameters can be changed: 

const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const bool MULTI = true;                          // MULTI = true  ... use multi-mesh,
                                                  // MULTI = false ... use single-mesh.
                                                  // Note: In the single mesh option, the meshes are
                                                  // forced to be geometrically the same but the
                                                  // polynomial degrees can still vary.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.3;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const std::string BDY_RIGHT = "1", BDY_TOP = "2";

// Problem parameters.
const double E  = 200e9;                          // Young modulus for steel: 200 GPa.
const double nu = 0.3;                            // Poisson ratio.
const double rho = 8000.0;                        // Density.
const double g1 = -9.81;                          // Gravitational acceleration.
const double f0  = 0;                             // Surface force in x-direction.
const double f1  = 1e3;                           // Surface force in y-direction.

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh u1_mesh, u2_mesh;
  H2DReader mloader;
  mloader.load("bracket.mesh", &u1_mesh);

  // Initial mesh refinements.
  u1_mesh.refine_element_id(1);
  u1_mesh.refine_element_id(4);

  // Create initial mesh for the vertical displacement component.
  // This also initializes the multimesh hp-FEM.
  u2_mesh.copy(&u1_mesh);

  // Initialize boundary conditions.
  DefaultEssentialBCConst zero_disp(BDY_RIGHT, 0.0);
  EssentialBCs bcs(&zero_disp);

  // Create x- and y- displacement space using the default H1 shapeset.
  H1Space u1_space(&u1_mesh, &bcs, P_INIT);
  H1Space u2_space(&u2_mesh, &bcs, P_INIT);
  info("ndof = %d.", Space::get_num_dofs(Hermes::vector<Space *>(&u1_space, &u2_space)));

  // Initialize the weak formulation.
  // NOTE; These weak forms are identical to those in example P01-linear/08-system.
  CustomWeakForm wf(E, nu, rho*g1, BDY_TOP, f0, f1);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&u1_space, &u2_space), is_linear);

  // Initialize coarse and reference mesh solutions.
  Solution u1_sln, u2_sln, u1_ref_sln, u2_ref_sln;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView s_view_0("Solution (x-displacement)", new WinGeom(0, 0, 400, 350));
  s_view_0.show_mesh(false);
  ScalarView s_view_1("Solution (y-displacement)", new WinGeom(760, 0, 400, 350));
  s_view_1.show_mesh(false);
  OrderView  o_view_0("Mesh (x-displacement)", new WinGeom(410, 0, 340, 350));
  OrderView  o_view_1("Mesh (y-displacement)", new WinGeom(1170, 0, 340, 350));
  ScalarView mises_view("Von Mises stress [Pa]", new WinGeom(0, 405, 400, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Hermes::vector<Space *>* ref_spaces = Space::construct_refined_spaces(Hermes::vector<Space *>(&u1_space, &u2_space));

    // Initialize matrix solver.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Assemble the reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, *ref_spaces, is_linear);
    dp->assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();
    
    // Solve the linear system of the reference problem. If successful, obtain the solutions.
    if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), *ref_spaces, 
                                            Hermes::vector<Solution *>(&u1_ref_sln, &u2_ref_sln));
    else error ("Matrix solver failed.\n");
  
    // Time measurement.
    cpu_time.tick();

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(Hermes::vector<Space *>(&u1_space, &u2_space), 
                                 Hermes::vector<Solution *>(&u1_ref_sln, &u2_ref_sln), 
                                 Hermes::vector<Solution *>(&u1_sln, &u2_sln), matrix_solver); 
   
    // View the coarse mesh solution and polynomial orders.
    s_view_0.show(&u1_sln); 
    o_view_0.show(&u1_space);
    s_view_1.show(&u2_sln); 
    o_view_1.show(&u2_space);
    // For von Mises stress Filter.
    double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
    double mu = E / (2*(1 + nu));
    VonMisesFilter stress(Hermes::vector<MeshFunction *>(&u1_sln, &u2_sln), lambda, mu);
    mises_view.show(&stress, HERMES_EPS_HIGH, H2D_FN_VAL_0, &u1_sln, &u2_sln, 1e4);

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Initialize adaptivity.
    Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&u1_space, &u2_space));

    /* 
    // Register custom forms for error calculation.
    adaptivity->set_error_form(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
    adaptivity->set_error_form(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
    adaptivity->set_error_form(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
    adaptivity->set_error_form(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
    */

    // Calculate error estimate for each solution component and the total error estimate.
    info("Calculating error estimate and exact error."); 
    Hermes::vector<double> err_est_rel;
    double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&u1_sln, &u2_sln), 
                               Hermes::vector<Solution *>(&u1_ref_sln, &u2_ref_sln), &err_est_rel) * 100;

    // Time measurement.
    cpu_time.tick();

    // Report results.
    info("ndof_coarse[0]: %d, ndof_fine[0]: %d, err_est_rel[0]: %g%%", 
         u1_space.Space::get_num_dofs(), Space::get_num_dofs((*ref_spaces)[0]), err_est_rel[0]*100);
    info("ndof_coarse[1]: %d, ndof_fine[1]: %d, err_est_rel[1]: %g%%",
         u2_space.Space::get_num_dofs(), Space::get_num_dofs((*ref_spaces)[1]), err_est_rel[1]*100);
    info("ndof_coarse_total: %d, ndof_fine_total: %d, err_est_rel_total: %g%%",
         Space::get_num_dofs(Hermes::vector<Space *>(&u1_space, &u2_space)), 
         Space::get_num_dofs(*ref_spaces), err_est_rel_total);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space::get_num_dofs(Hermes::vector<Space *>(&u1_space, &u2_space)), 
                             err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < ERR_STOP) 
      done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&selector, &selector), 
                               THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (Space::get_num_dofs(Hermes::vector<Space *>(&u1_space, &u2_space)) >= NDOF_STOP) done = true;

    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false)
      for(unsigned int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i]->get_mesh();
    delete ref_spaces;
    delete dp;
    
    // Increase counter.
    as++;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  s_view_0.set_title("Fine mesh solution (x-displacement)");
  s_view_0.show(&u1_ref_sln);
  s_view_1.set_title("Fine mesh solution (y-displacement)");
  s_view_1.show(&u2_ref_sln);
  // For von Mises stress Filter.
  double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
  double mu = E / (2*(1 + nu));
  VonMisesFilter stress(Hermes::vector<MeshFunction *>(&u1_ref_sln, &u2_ref_sln), lambda, mu);
  mises_view.show(&stress, HERMES_EPS_HIGH, H2D_FN_VAL_0, &u1_ref_sln, &u2_ref_sln, 1e4);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

