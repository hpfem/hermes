#define H2D_REPORT_INFO
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This is a long version of example 10-adapt: function solve_linear_adapt() is not used.

const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.2;                     // This is a quantitative parameter of the adapt(...) function and
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
const CandList CAND_LIST = H2D_HP_ANISO_H;        // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                                  // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
const int OMEGA_1 = 1;
const int OMEGA_2 = 2;
const int STATOR_BDY = 2;
const double EPS_1 = 1.0;       // Relative electric permittivity in Omega_1.
const double EPS_2 = 10.0;      // Relative electric permittivity in Omega_2.
const double VOLTAGE = 50.0;    // Voltage on the stator.

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return (ess_bdy_marker == STATOR_BDY) ? VOLTAGE : 0.0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("motor.mesh", &mesh);

  // Create an H1 space with default shapeset.
  H1Space* space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(biform1), H2D_SYM, OMEGA_1);
  wf.add_matrix_form(callback(biform2), H2D_SYM, OMEGA_2);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, get_num_dofs(space), mat, rhs, solver);

  // Initialize views.
  ScalarView s_view("Solution", new WinGeom(0, 0, 400, 600));
  OrderView  o_view("Mesh", new WinGeom(410, 0, 400, 600));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Adaptivity loop.
  Solution *sln = new Solution();
  Solution *ref_sln = new Solution();
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    info("Solving on reference mesh.");

    // Construct globally refined reference mesh
    // and setup reference space.
    Mesh *ref_mesh = new Mesh();
    ref_mesh->copy(space->get_mesh());
    ref_mesh->refine_all_elements();
    Space* ref_space = space->dup(ref_mesh);
    int order_increase = 1;
    ref_space->copy_orders(space, order_increase);

    // Solve the reference problem.
    solve_linear(ref_space, &wf, matrix_solver, ref_sln);

    // Project the reference solution on the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    // NULL means that we do not want to know the resulting coefficient vector.
    project_global(space, H2D_H1_NORM, ref_sln, sln, NULL); 

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution.
    s_view.show(sln); 
    o_view.show(space);

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors.
    info("Calculating error (est).");
    Adapt hp(space, H2D_H1_NORM);
    hp.set_solutions(sln, ref_sln);
    double err_est_rel = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;
 
    // Report results.
    info("ndof: %d, ref_ndof: %d, err_est_rel_total: %g%%", 
         get_num_dofs(space), get_num_dofs(ref_space), err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(get_num_dofs(space), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else {
      info("Adapting the coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (get_num_dofs(space) >= NDOF_STOP) done = true;
    }

    // Free reference space and mesh.
    ref_space->free();

    as++;
  }
  while (done == false);

  info("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

