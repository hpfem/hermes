#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

//  This is the ninth in the series of NIST benchmarks with known exact solutions. This benchmark
//  has four different versions, use the global variable PROB_PARAM below to switch among them.
//  It differs from 09-wave-front in the mesh adaptation method.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  PDE: -Laplace u - f = 0
//
//  Known exact solution; atan(ALPHA * (sqrt(pow(x - X_LOC, 2) + pow(y - Y_LOC, 2)) - R_ZERO));
//  See the class CustomExactSolution.
//
//  Domain: unit square (0, 1) x (0, 1), see the file square.mesh.
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

int PARAM = 3;         // PARAM determines which parameter values you wish to use 
                       //            for the steepness and location of the wave front. 
                       // #| name   |   ALPHA | X_LOC	| Y_LOC | R_ZERO
                       // 0: mild		    20      -0.05  -0.05    0.7
                       // 1: steep      1000    -0.05  -0.05    0.7
                       // 2: asymmetric 1000     1.5    0.25    0.92
                       // 3: well       50       0.5    0.5     0.25

const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
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
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const bool USE_RESIDUAL_ESTIMATOR = true;         // Add also the norm of the residual to the error estimate of each element.

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Define problem parameters: (x_loc, y_loc) is the center of the circular wave front, R_ZERO is the distance from the 
  // wave front to the center of the circle, and alpha gives the steepness of the wave front.
  double alpha, x_loc, y_loc, r_zero;
  switch(PARAM) {
  case 0:
    alpha = 20;
    x_loc = -0.05;
    y_loc = -0.05;
    r_zero = 0.7;
    break;
  case 1:
    alpha = 1000;
    x_loc = -0.05;
    y_loc = -0.05;
    r_zero = 0.7;
    break;
  case 2:
    alpha = 1000;
    x_loc = 1.5;
    y_loc = 0.25;
    r_zero = 0.92;
    break;
  case 3:
    alpha = 50;
    x_loc = 0.5;
    y_loc = 0.5;
    r_zero = 0.25;
    break;
  default:   // The same as 0.
    alpha = 20;
    x_loc = -0.05;
    y_loc = -0.05;
    r_zero = 0.7;
    break;
  }

  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh);     // quadrilaterals
  // mloader.load("square_tri.mesh", &mesh);   // triangles

  // Perform initial mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  
  // Set exact solution.
  CustomExactSolution exact(&mesh, alpha, x_loc, y_loc, r_zero);

  // Define right-hand side.
  CustomRightHandSide rhs(alpha, x_loc, y_loc, r_zero);

  // Initialize the weak formulation.
  DefaultWeakFormPoisson wf(&rhs);

  // Initialize boundary conditions.
  DefaultEssentialBCNonConst bc("Bdy", &exact);
  EssentialBCs bcs(&bc);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bcs, P_INIT);
  
  // Initialize approximate solution.
  Solution sln;

  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 440, 350));
  sview.show_mesh(false);
  sview.fix_scale_width(50);
  OrderView  oview("Polynomial orders", new WinGeom(450, 0, 420, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu, graph_dof_exact, graph_cpu_exact;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs_vec = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs_vec);

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Assemble the discrete problem.
    info("Solving.");
    DiscreteProblem* dp = new DiscreteProblem(&wf, &space);

    // Time measurement.
    cpu_time.tick();

    // Initial coefficient vector for the Newton's method.  
    int ndof = Space::get_num_dofs(&space);
    scalar* coeff_vec = new scalar[ndof];
    memset(coeff_vec, 0, ndof * sizeof(scalar));

    // Perform Newton's iteration.
    if (!hermes2d.solve_newton(coeff_vec, dp, solver, matrix, rhs_vec)) error("Newton's iteration failed.");

    // Translate the resulting coefficient vector into the Solution sln.
    Solution ref_sln;
    Solution::vector_to_solution(coeff_vec, &space, &sln);

    // View the approximate solution and polynomial orders.
    sview.show(&sln);
    oview.show(&space);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate and exact error.");
    BasicKellyAdapt* adaptivity = new BasicKellyAdapt(&space);
    
    // Use energy norm for error estimate normalization and measuring of exact error.    
    adaptivity->set_error_form(new EnergyErrorForm(&wf));

    if (USE_RESIDUAL_ESTIMATOR) 
      adaptivity->add_error_estimator_vol(new ResidualErrorForm(&rhs));
    
    double err_est_rel = adaptivity->calc_err_est(&sln) * 100;  
    double err_exact_rel = adaptivity->calc_err_exact(&sln, &exact, false) * 100;

    // Time measurement.
    cpu_time.tick();
    
    // Report results.
    info("ndof_coarse: %d", Space::get_num_dofs(&space));
    info("err_est_rel: %g%%, err_exact_rel: %g%%", err_est_rel, err_exact_rel);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
    graph_dof_exact.add_values(Space::get_num_dofs(&space), err_exact_rel);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel);
    graph_cpu_exact.save("conv_cpu_exact.dat");

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
    delete [] coeff_vec;
    delete adaptivity;
    delete dp;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
