#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example shows how to adapt the mesh to match an arbitrary 
// given function (no PDE solved). The definition of the function 
// f() is below.
//
// The following parameters can be changed:

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
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
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

// This function can be modified.
scalar f(double x, double y, double& dx, double& dy)
{
  dx = 0.25 * pow(x*x + y*y, -0.75) * 2 * x;
  dy = 0.25 * pow(x*x + y*y, -0.75) * 2 * y;
  return pow(x*x + y*y, 0.25);
}

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, NULL, NULL, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf_dummy;

  // Initialize views.
  ScalarView sview("Scalar potential Phi", 0, 0, 600, 300);
  OrderView  oview("Mesh", 620, 0, 600, 300);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, get_num_dofs(&space), mat, rhs, solver);

  // Adaptivity loop:
  Solution sln_coarse, sln_fine;
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh
    // and setup reference space.
    Mesh *ref_mesh = new Mesh();
    ref_mesh->copy(space.get_mesh());
    ref_mesh->refine_all_elements();
    Space* ref_space = space.dup(ref_mesh);
    int order_increase = 1;
    ref_space->copy_orders(&space, order_increase);


    // Assign the function f() to the fine mesh.
    sln_fine.set_exact(ref_mesh, f);

    // Project the function f() on the coarse mesh.
    Solution sln_tmp;
    sln_tmp.set_exact(&mesh, f);
    project_global(&space, H2D_H1_NORM, &sln_tmp, &sln_coarse);
 
    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution.
    sview.show(&sln_coarse);
    oview.show(&space);

    // Time measurement.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error.");
    Adapt hp(&space, H2D_H1_NORM);
    hp.set_solutions(&sln_coarse, &sln_fine);
    double err_est_rel = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;

    // Report results.
    info("ndof_coarse: %d, err_est: %g%%", get_num_dofs(&space), err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(get_num_dofs(&space), err_est_rel);
    graph_dof.save("conv_dof.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else {
      info("Adapting coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (get_num_dofs(&space) >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the fine mesh solution - the final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(&sln_fine);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

