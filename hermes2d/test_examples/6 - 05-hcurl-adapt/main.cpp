#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

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
//          "lshape3t.mesh" (triangular mesh). See the mesh.load(...) command below.
//
//  BC: perfect conductor on boundary markers 1 and 6 (essential BC)
//      impedance boundary condition on rest of boundary (natural BC)
//
//  The following parameters can be changed:

const int P_INIT = 2;                             // Initial polynomial degree. NOTE: The meaning is different from
                                                  // standard continuous elements in the space H1. Here, P_INIT refers
                                                  // to the maximum poly order of the tangential component, and polynomials
                                                  // of degree P_INIT + 1 are present in element interiors. P_INIT = 0
                                                  // is for Whitney elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
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
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double MU_R   = 1.0;
const double KAPPA  = 1.0;
const double LAMBDA = 1.0;

// Bessel functions, exact solution, and weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Time measurement
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("lshape3q.mesh", &mesh);    // quadrilaterals
  //mloader.load("lshape3t.mesh", &mesh);  // triangles

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_essential(Hermes::vector<std::string>("Corner horizontal",
                                                                   "Corner vertical"), 0);
  EssentialBCs bcs(&bc_essential);

  // Create an Hcurl space with default shapeset.
  HcurlSpace space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakForm wf(MU_R, KAPPA);

  // Initialize coarse and reference mesh solutions.
  Solution sln, ref_sln;

  // Initialize exact solution.
  CustomExactSolution sln_exact(&mesh);

  // Initialize refinement selector.
  HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  VectorView v_view("Solution (magnitude)", new WinGeom(0, 0, 460, 350));
  v_view.set_min_max_range(0, 1.5);
  OrderView  o_view("Polynomial orders", new WinGeom(470, 0, 400, 350));
  
  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, 
              graph_dof_exact, graph_cpu_exact;
  
  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = Space::construct_refined_space(&space);
    int ndof_ref = Space::get_num_dofs(ref_space);

    // Initialize matrix solver.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Initialize reference problem.
    info("Solving on reference mesh.");
    DiscreteProblem dp(&wf, ref_space);

    // Time measurement.
    cpu_time.tick();

    // Initial coefficient vector for the Newton's method.  
    scalar* coeff_vec = new scalar[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(scalar));

    // Perform Newton's iteration.
    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs)) error("Newton's iteration failed.");

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solution(coeff_vec, ref_space, &ref_sln);

    // Time measurement.
    cpu_time.tick();

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 
   
    // View the coarse mesh solution and polynomial orders.
    v_view.show(&sln);
    o_view.show(&space);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate and exact error."); 
    Adapt* adaptivity = new Adapt(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // Calculate exact error.
    bool solutions_for_adapt = false;
    double err_exact_rel = adaptivity->calc_err_exact(&sln, &sln_exact, solutions_for_adapt) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d", 
      Space::get_num_dofs(&space), Space::get_num_dofs(ref_space));
    info("err_est_rel: %g%%, err_exact_rel: %g%%", err_est_rel, err_exact_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");
    graph_dof_exact.add_values(Space::get_num_dofs(&space), err_exact_rel);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel);
    graph_cpu_exact.save("conv_cpu_exact.dat");

    // If err_est_rel too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      
      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;

    // Clean up.
    delete [] coeff_vec;
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false) delete ref_space->get_mesh();
    delete ref_space;
  }
  while (done == false);
  
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  v_view.set_title("Fine mesh solution (magnitude)");
  v_view.show(&ref_sln);
  
  // Wait for all views to be closed.
  View::wait();
  return 0;
}

