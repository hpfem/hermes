#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This test makes sure that an exact function 
// sin(x) is approximated adaptively with the 
// relative error 1e-3% with less than 20 DOF.
// Adaptivity is done in H1 norm.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

//  The following parameters can be changed:
static int NEQ = 1;
int NELEM = 2;                          // Number of elements.
double A = -M_PI, B = M_PI;             // Domain end points.
int P_INIT = 1;                         // Initial polynomial degree.

// Newton's method.
double NEWTON_TOL_COARSE = 1e-6;        // Coarse mesh.
double NEWTON_TOL_REF = 1e-6;           // Fine mesh.
int NEWTON_MAX_ITER = 150;

// Adaptivity.
const int ADAPT_TYPE = 0;               // 0... hp-adaptivity.
                                        // 1... h-adaptivity.
                                        // 2... p-adaptivity.
const double THRESHOLD = 0.7;           // Refined will be all elements whose error 
                                        // is greater than THRESHOLD*max_elem_error.
const double TOL_ERR_REL = 1e-6;        // Tolerance for the relative error between 
                                        // the coarse mesh and fine solutions.
const int NORM = 1;                     // To measure errors.
                                        // 1... H1 norm.
                                        // 0... L2 norm.

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary conditions.
Hermes::vector<BCSpec *>DIR_BC_LEFT =  Hermes::vector<BCSpec *>(new BCSpec(0,0));
Hermes::vector<BCSpec *>DIR_BC_RIGHT = Hermes::vector<BCSpec *>(new BCSpec(0,0));

// Function f(x).
double f(double x) 
{
  return sin(x);
}

// Exact solution.
const int EXACT_SOL_PROVIDED = 1;
void exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) 
{
  u[0] = sin(x);
  dudx[0] = cos(x);
}


// bilinear form for the Jacobi matrix .
// num...number of Gauss points in element.
// x[]...Gauss points.
// weights[]...Gauss weights for points in x[].
// u...basis function.
// v...test function.
// u_prev...previous solution (all solution components).
double jacobian(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) 
  {
    val += dudx[i]*dvdx[i]*weights[i];
  }
  return val;
};

// (nonlinear) form for the residual vector.
// num...number of Gauss points in element.
// x[]...Gauss points.
// weights[]...Gauss weights for points in x[].
// u...approximate solution.
// v...test function.
// u_prev...previous solution (all solution components).
double residual(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  double val = 0;
  // Solution index (only 0 is relevant for this example).
  int si = 0;
  for(int i = 0; i<num; i++) 
  {
    val += (du_prevdx[si][0][i]*dvdx[i] - f(x[i])*v[i])*weights[i];
  }
  return val;
};


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
  WeakForm wf;
  wf.add_matrix_form(jacobian);
  wf.add_vector_form(residual);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem *dp_coarse = new DiscreteProblem(&wf, space, is_linear);

  // Newton's loop on coarse mesh.
  // Fill vector coeff_vec using dof and coeffs arrays in elements.
  double *coeff_vec_coarse = new double[Space::get_num_dofs(space)];
  get_coeff_vector(space, coeff_vec_coarse);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
  Vector* rhs_coarse = create_vector(matrix_solver);
  Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

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
    // NOTE: at least one full iteration forced
    //       here because sometimes the initial
    //       residual on fine mesh is too small.
    if(res_l2_norm < NEWTON_TOL_COARSE && it > 1) break;

    // Multiply the residual vector with -1 since the matrix 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    for(int i=0; i<ndof_coarse; i++) rhs_coarse->set(i, -rhs_coarse->get(i));

    // Solve the linear system.
    if(!solver_coarse->solve())
      error ("Matrix solver failed.\n");

      // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof_coarse; i++) coeff_vec_coarse[i] += solver_coarse->get_solution()[i];

    // If the maximum number of iteration has been reached, then quit.
    if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
    
    // Copy coefficients from vector y to elements.
    set_coeff_vector(coeff_vec_coarse, space);
    
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

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = construct_refined_space(space);

    // Initialize the FE problem.
    bool is_linear = false;
    DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space, is_linear);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Newton's loop on the fine mesh.
    info("Solving on fine mesh:");
    // Fill vector coeff_vec using dof and coeffs arrays in elements.
    double *coeff_vec = new double[Space::get_num_dofs(ref_space)];
    get_coeff_vector(ref_space, coeff_vec);

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
      // NOTE: at least one full iteration forced
      //       here because sometimes the initial
      //       residual on fine mesh is too small.
      if(res_l2_norm < NEWTON_TOL_REF && it > 1) break;

      // Multiply the residual vector with -1 since the matrix 
      // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
      for(int i=0; i<ndof; i++) rhs->set(i, -rhs->get(i));

      // Solve the linear system.
      if(!solver->solve())
        error ("Matrix solver failed.\n");

        // Add \deltaY^{n+1} to Y^n.
      for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];

      // If the maximum number of iteration has been reached, then quit.
      if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
      
      // Copy coefficients from vector y to elements.
      set_coeff_vector(coeff_vec, ref_space);

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
    if (EXACT_SOL_PROVIDED) {
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

    // Decide whether the relative error is sufficiently small.
    if(err_est_rel < TOL_ERR_REL) break;

    // Returns updated coarse and fine meshes, with the last 
    // coarse and fine mesh solutions on them, respectively. 
    // The coefficient vectors and numbers of degrees of freedom 
    // on both meshes are also updated. 
    adapt(NORM, ADAPT_TYPE, THRESHOLD, err_est_array, space, ref_space);

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

  // Test variable.
  int success_test = 1; 
  info("N_dof = %d.", Space::get_num_dofs(space));
  if (Space::get_num_dofs(space) > 40) success_test = 0;

  if (success_test) 
  {
    info("Success!");
    return ERROR_SUCCESS;
  }
  else 
  {
    info("Failure!");
    return ERROR_FAILURE;
  }
}
