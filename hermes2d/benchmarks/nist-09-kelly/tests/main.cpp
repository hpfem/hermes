#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

/** \addtogroup t_bench_nist-09-kelly benchmarks/nist-09-kelly
 *  \{
 *  \brief This test makes sure that the benchmark nist-09-kelly works correctly.
 *  (Testing version NIST-09 "well")
 *
 *  \section s_params Parameters
 *   - PROB_PARAM = 3
 *   - INIT_REF_NUM=1
 *   - P_INIT=1
 *   - THRESHOLD=0.3
 *   - STRATEGY=0
 *   - MESH_REGULARITY=-1
 *   - USE_RESIDUAL_ESTIMATOR = false
 *   - ERR_STOP=1.0
 *   - NDOF_STOP=60000
 *   - matrix_solver = SOLVER_UMFPACK
 *
 *
 *  \section s_res Results
 *   - DOFs: 15089
 *   - Adaptivity steps: 14
 *   - Exact error: 4.07387%
 */

int PROB_PARAM = 3;    // PROB_PARAM determines which parameter values you wish to use for the steepness and location of the wave front. 
                       //    name               ALPHA   X_LOC   Y_LOC   R_ZERO
                       // 0: mild               20      -0.05   -0.05   0.7
                       // 1: steep              1000    -0.05   -0.05   0.7
                       // 2: asymmetric         1000     1.5     0.25   0.92
                       // 3: well               50       0.5     0.5    0.25

const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
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
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const bool USE_RESIDUAL_ESTIMATOR = false;        // Add also the norm of the residual to the error estimate of each element.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 10000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
double ALPHA;          // (X_LOC, Y_LOC) is the center of the circular wave front, R_ZERO is the distance from the 
double X_LOC;          // wave front to the center of the circle, and ALPHA gives the steepness of the wave front.
double Y_LOC;
double R_ZERO;

// Exact solution.
#include "exact_solution.cpp"

// Boundary markers.
const int BDY_DIRICHLET = 1;

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(double x, double y) { return fn(x, y);}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh);     // quadrilaterals
  // mloader.load("square_tri.mesh", &mesh);   // triangles

  // Perform initial mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boudnary values.
  BCValues bc_values;
  bc_values.add_function(BDY_DIRICHLET, essential_bc_values);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form), HERMES_SYM);
  wf.add_vector_form(callback(linear_form));
  
  // Initialize approximate solution.
  Solution sln;

  // Set exact solution.
  ExactSolution exact(&mesh, fndd);

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
  
  // Initialize the matrix solver.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  double err_exact_rel;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Assemble the discrete problem.
    info("Solving.");
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, &space, is_linear);
    dp->assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();

    // Solve the linear system. If successful, obtain the solution.
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
    else error ("Matrix solver failed.\n");

    // Calculate element errors and total error estimate.
    info("Calculating error estimate and exact error.");
    BasicKellyAdapt* adaptivity = new BasicKellyAdapt(&space);
    
    if (USE_RESIDUAL_ESTIMATOR)
      adaptivity->add_error_estimator_vol(callback(residual_estimator));
    
    double err_est_rel = adaptivity->calc_err_est(&sln) * 100;  
    err_exact_rel = adaptivity->calc_err_exact(&sln, &exact) * 100;   
    
    // Time measurement.
    cpu_time.tick();
    
    // Report results.
    info("ndof_coarse: %d", Space::get_num_dofs(&space));
    info("err_est_rel: %1.15g%%, err_exact_rel: %1.15g%%", err_est_rel, err_exact_rel);
    
    // This is to ensure that the two possible approaches to interface error estimators accumulation give
    // same results.
    KellyTypeAdapt* adaptivity2 = new KellyTypeAdapt(&space, Hermes::vector<ProjNormType>(), false);
    adaptivity2->disable_aposteriori_interface_scaling();
    adaptivity2->add_error_estimator_surf(callback(interface_estimator));

    if (USE_RESIDUAL_ESTIMATOR)
      adaptivity2->add_error_estimator_vol(callback(residual_estimator));
    
    double err_est_rel2 = adaptivity2->calc_err_est(&sln) * 100;  
    double err_exact_rel2 = adaptivity2->calc_err_exact(&sln, &exact) * 100;
    
    info("err_est_rel_2: %1.15g%%, err_exact_rel_2: %1.15g%%", err_est_rel2, err_exact_rel2);
    
    if (fabs(err_est_rel2 - err_est_rel) >= 1e-13 || fabs(err_exact_rel2 - err_exact_rel) >= 1e-13)
    {
      info("err_est_rel diff: %g, err_exact_rel diff: %g", 
           std::abs(err_est_rel2 - err_est_rel), std::abs(err_exact_rel2 - err_exact_rel));
      err_est_rel = 0;      // to immediately exit the adaptivity loop
      err_exact_rel = 1e20; // to fail the test
    }
     
    // Time measurement.
    cpu_time.tick(HERMES_SKIP);
    
    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else
    {
      info("Adapting the mesh.");
      done = adaptivity->adapt(THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;
    
    // Increase the counter of performed adaptivity steps.
    if (done == false)  as++;

    // Clean up.
    delete adaptivity;
    delete adaptivity2;
    delete dp;
  }
  while (done == false);
  
  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  verbose("Total running time: %g s", cpu_time.accumulated());

  int ndof = Space::get_num_dofs(&space);

  int n_dof_allowed = 15555;
  double err_exact_rel_allowed = 4.075;
  printf("n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);
  printf("err_exact_rel_actual = %g%%\n", err_exact_rel);
  printf("err_exact_rel_allowed = %g%%\n", err_exact_rel_allowed);
  if (ndof <= n_dof_allowed && err_exact_rel <= err_exact_rel_allowed) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
  
  // Wait for all views to be closed.
  View::wait();
  return 0;
}
