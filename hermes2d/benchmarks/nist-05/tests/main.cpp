#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

/** \addtogroup t_bench_nist-05 Benchmarks/nist-05
 *  \{
 *  \brief This test makes sure that the benchmark "nist-05" works correctly.
 *
 *  \section s_params Parameters
 *   - INIT_REF_NUM = 0
 *   - P_INIT = 3
 *   - THRESHOLD = 0.3
 *   - STRATEGY = 0
 *   - CAND_LIST = H2D_HP_ANISO;
 *   - MESH_REGULARITY = -1
 *   - CONV_EXP = 0.5
 *   - ERR_STOP = 0.1
 *   - NDOF_STOP = 60000
 *   - matrix_solver = SOLVER_UMFPACK
 *
 *   - DOFs: 3581
 *   - Adaptivity steps: 32
 */

const int P_INIT = 3;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
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
const double ERR_STOP = 0.1;                      // Stopping criterion for adaptivity (rel. error tolerance between the
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const int OMEGA_1 = 1;
const int OMEGA_2 = 2;
const int OMEGA_3 = 3;
const int OMEGA_4 = 4;
const int OMEGA_5 = 5;

const double P_1 = 25.0;
const double P_2 = 7.0;
const double P_3 = 5.0;
const double P_4 = 0.2;
const double P_5 = 0.05;

const double Q_1 = 25.0;
const double Q_2 = 0.8;
const double Q_3 = 0.0001;
const double Q_4 = 0.2;
const double Q_5 = 0.05;

const double F_1 = 0.0;
const double F_2 = 1.0;
const double F_3 = 1.0;
const double F_4 = 0.0;
const double F_5 = 0.0;

// Boundary markers.
const int BDY_LEFT = 1;
const int BDY_TOP = 2;
const int BDY_RIGHT = 3;
const int BDY_BOTTOM = 4;

// Boundary condition coefficients for the four sides.
const double C_LEFT = 0.0;
const double C_TOP = 1.0;
const double C_RIGHT = 2.0;
const double C_BOTTOM = 3.0;

const double G_N_LEFT = 0.0;
const double G_N_TOP = 3.0;
const double G_N_RIGHT = 2.0;
const double G_N_BOTTOM = 1.0;

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("battery.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_neumann(BDY_LEFT);
  bc_types.add_bc_newton(Hermes::vector<int>(BDY_TOP, BDY_RIGHT, BDY_BOTTOM));

  // Enter Dirichlet boudnary values.
  BCValues bc_values;

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(biform1), HERMES_SYM, OMEGA_1);
  wf.add_matrix_form(callback(biform2), HERMES_SYM, OMEGA_2);
  wf.add_matrix_form(callback(biform3), HERMES_SYM, OMEGA_3);
  wf.add_matrix_form(callback(biform4), HERMES_SYM, OMEGA_4);
  wf.add_matrix_form(callback(biform5), HERMES_SYM, OMEGA_5);

  wf.add_matrix_form_surf(bilinear_form_surf_right, bilinear_form_ord, BDY_RIGHT);
  wf.add_matrix_form_surf(bilinear_form_surf_top, bilinear_form_ord, BDY_TOP);
  wf.add_matrix_form_surf(bilinear_form_surf_bottom, bilinear_form_ord, BDY_BOTTOM);

  wf.add_vector_form_surf(callback(linear_form_surf_left), BDY_LEFT);
  wf.add_vector_form_surf(callback(linear_form_surf_right), BDY_RIGHT);
  wf.add_vector_form_surf(callback(linear_form_surf_top), BDY_TOP);
  wf.add_vector_form_surf(callback(linear_form_surf_bottom), BDY_BOTTOM);

  wf.add_vector_form(callback(linear_form_1), OMEGA_1);
  wf.add_vector_form(callback(linear_form_2), OMEGA_2);
  wf.add_vector_form(callback(linear_form_3), OMEGA_3);
  wf.add_vector_form(callback(linear_form_4), OMEGA_4);
  wf.add_vector_form(callback(linear_form_5), OMEGA_5);

  // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln;
  
  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = construct_refined_space(&space);

    // Initialize matrix solver.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Assemble reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space, is_linear);
    dp->assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();
    
    // Solve the linear system of the reference problem. 
    // If successful, obtain the solution.
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), ref_space, &ref_sln);
    else error ("Matrix solver failed.\n");

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

    // Time measurement.
    cpu_time.tick();
   
    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt* adaptivity = new Adapt(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
      Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
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
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false) delete ref_space->get_mesh();
    delete ref_space;
    delete dp;
    
  }
  while (done == false);
  
  verbose("Total running time: %g s", cpu_time.accumulated());

  int ndof = Space::get_num_dofs(&space);

  int n_dof_allowed = 3600;
  printf("n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);
  if (ndof <= n_dof_allowed) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }

}

