#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

//  This example solves a system of two linear second-order equations.
//
//  PDE: u' + k^2 v = 0
//       u - v' = 0,
//  which is equivalent to u'' + k^2 u = 0.
//
//  Interval: (0, 2*pi).
//
//  BC: Dirichlet, u(0) = 0, v(0) = k.
//
//  Exact solution: u(x) = sin(k*x), v(x) = k*cos(k*x).
//
//  The following parameters can be changed:

static int NEQ = 2;
int NELEM = 3;                          // Number of elements.
int P_INIT = 1;                         // Initial polynomial degree.
const int ADAPT_TYPE = 0;               // 0... hp-adaptivity.
                                        // 1... h-adaptivity.
                                        // 2... p-adaptivity.
const double THRESHOLD = 0.7;           // Refined will be all elements whose error 
                                        // is greater than THRESHOLD*max_elem_error.
const double TOL_ERR_REL = 1e-1;        // Tolerance for the relative error between 
                                        // the coarse mesh and fine solutions.
const int NORM = 1;                     // To measure errors.
                                        // 1... H1 norm.
                                        // 0... L2 norm.
double NEWTON_TOL_COARSE = 1e-6;        // Newton tolerance on coarse mesh.
double NEWTON_TOL_REF = 1e-6;           // Newton tolerance on fine mesh.
int NEWTON_MAX_ITER = 150;              // Maximum number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Equation parameters.
double K = 1.0;                         // Equation parameter.
double A = 0, B = 2*M_PI;               // Domain end points.

// Boundary conditions.
Hermes::vector<BCSpec *>DIR_BC_LEFT =  Hermes::vector<BCSpec *>(new BCSpec(0,0), new BCSpec(1,K));
Hermes::vector<BCSpec *> DIR_BC_RIGHT = Hermes::vector<BCSpec *>();

// Exact solution.
// When changing exact solution, do not 
// forget to update interval accordingly.
const int EXACT_SOL_PROVIDED = 1;
void exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) 
{
  u[0] = sin(K*x);
  dudx[0] = K*cos(K*x);
  u[1] = K*cos(K*x);
  dudx[1] = -K*K*sin(K*x);
}

// Weak forms for Jacobi matrix and residual.
#include "definitions.cpp"

int main() 
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Create coarse mesh, set Dirichlet BC, enumerate basis functions.
  Space* space = new Space(A, B, NELEM, DIR_BC_LEFT, DIR_BC_RIGHT, P_INIT, NEQ, NEQ);

  // Enumerate basis functions, info for user.
  int ndof = Space::get_num_dofs(space);
  info("ndof: %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, jacobian_0_0);
  wf.add_matrix_form(0, 1, jacobian_0_1);
  wf.add_matrix_form(1, 0, jacobian_1_0);
  wf.add_matrix_form(1, 1, jacobian_1_1);
  wf.add_vector_form(0, residual_0);
  wf.add_vector_form(1, residual_1);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem *dp_coarse = new DiscreteProblem(&wf, space, is_linear);

  // Set zero initial condition.
  double *coeff_vec_coarse = new double[Space::get_num_dofs(space)];
  set_zero(coeff_vec_coarse, Space::get_num_dofs(space));

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
  Vector* rhs_coarse = create_vector(matrix_solver);
  Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

  // Newton's loop on coarse mesh.
  int it = 1;
  while (1) 
  {
    // Obtain the number of degrees of freedom.
    int ndof_coarse = Space::get_num_dofs(space);

    // Assemble the Jacobian matrix and residual vector.
    dp_coarse->assemble(coeff_vec_coarse, matrix_coarse, rhs_coarse);

    // Calculate the l2-norm of residual vector.
    double res_l2_norm = get_l2_norm(rhs_coarse);

    // Info for user.
    info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(space), res_l2_norm);

    // If l2 norm of the residual vector is within tolerance, then quit.
    // NOTE: at least one full iteration forced here because sometimes the initial
    //       residual on fine mesh is already too small.
    if(res_l2_norm < NEWTON_TOL_COARSE && it > 1) break;

    // Multiply the residual vector with -1 since the matrix 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    for(int i = 0; i < ndof_coarse; i++) rhs_coarse->set(i, -rhs_coarse->get(i));

    // Solve the linear system.
    if(!solver_coarse->solve()) error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof_coarse; i++) coeff_vec_coarse[i] += solver_coarse->get_solution()[i];

    // If the maximum number of iteration has been reached, then quit.
    if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
        
    it++;
  }
  
  // Cleanup.
  delete matrix_coarse;
  delete rhs_coarse;
  delete solver_coarse;
  delete [] coeff_vec_coarse;
  delete dp_coarse;

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;
  SimpleGraph graph_dof_exact, graph_cpu_exact;

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as); 

    // Construct globally refined reference mesh, create reference space on it,
    // and transfer the coarse mesh solution there.
    Space* ref_space = construct_refined_space(space);
    
    // Initialize the FE problem.
    bool is_linear = false;
    DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space, is_linear);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Fill vector coeff_vec using dof and coeffs arrays in elements.
    double *coeff_vec = new double[Space::get_num_dofs(ref_space)];
    get_coeff_vector(ref_space, coeff_vec);

    // Newton's loop on the fine mesh.
    info("Solving on fine mesh:");
    int it = 1;
    while (1) 
    {
      // Obtain the number of degrees of freedom.
      int ndof = Space::get_num_dofs(ref_space);

      // Assemble the Jacobian matrix and residual vector.
      dp->assemble(coeff_vec, matrix, rhs);

      // Calculate the l2-norm of residual vector.
      double res_l2_norm = get_l2_norm(rhs);

      // Info for user.
      info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(ref_space), res_l2_norm);

      // If l2 norm of the residual vector is within tolerance, then quit.
      // NOTE: at least one full iteration forced here because sometimes the initial
      //       residual on fine mesh is already too small.
      if(res_l2_norm < NEWTON_TOL_REF && it > 1) break;

      // Multiply the residual vector with -1 since the matrix 
      // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
      for(int i=0; i<ndof; i++) rhs->set(i, -rhs->get(i));

      // Solve the linear system.
      if(!solver->solve()) error ("Matrix solver failed.\n");

      // Add \deltaY^{n+1} to Y^n.
      for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];

      // If the maximum number of iteration has been reached, then quit.
      if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
      
      it++;
    }
    
    // Starting with second adaptivity step, obtain new coarse 
    // mesh solution via projecting the fine mesh solution.
    if(as > 1)
    {
      info("Projecting the fine mesh solution onto the coarse mesh.");
      // Project the fine mesh solution (defined on space_ref) onto the coarse mesh (defined on space).
      OGProjection::project_global(space, ref_space, matrix_solver);
    }

    // Calculate element errors and total error estimate.
    info("Calculating error estimate.");
    double err_est_array[MAX_ELEM_NUM];
    double err_est_rel = calc_err_est(NORM, space, ref_space, err_est_array) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
      Space::get_num_dofs(space), Space::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // If exact solution available, also calculate exact error.
    if (EXACT_SOL_PROVIDED) 
    {
      // Calculate element errors wrt. exact solution.
      double err_exact_rel = calc_err_exact(NORM, space, exact_sol, NEQ, A, B) * 100;
     
      // Info for user.
      info("Relative error (exact) = %g %%", err_exact_rel);

      // Add entry to DOF and CPU convergence graphs.
      graph_dof_exact.add_values(Space::get_num_dofs(space), err_exact_rel);
      graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel);
    }

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space::get_num_dofs(space), err_est_rel);
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);

    // If err_est_rel too large, adapt the mesh.
    if (err_est_rel < NEWTON_TOL_REF) done = true;
    else 
    {
      info("Adapting the coarse mesh.");
      adapt(NORM, ADAPT_TYPE, THRESHOLD, err_est_array, space, ref_space);
    }

    as++;
    
    // Plot meshes, results, and errors.
    adapt_plotting(space, ref_space, NORM, EXACT_SOL_PROVIDED, exact_sol);

    // Cleanup.
    delete solver;
    delete matrix;
    delete rhs;
    delete ref_space;
    delete dp;
    delete [] coeff_vec;
  }
  while (done == false);

  info("Total running time: %g s", cpu_time.accumulated());

  // Save convergence graphs.
  graph_dof_est.save("conv_dof_est.dat");
  graph_cpu_est.save("conv_cpu_est.dat");
  graph_dof_exact.save("conv_dof_exact.dat");
  graph_cpu_exact.save("conv_cpu_exact.dat");

  return 0;
}
