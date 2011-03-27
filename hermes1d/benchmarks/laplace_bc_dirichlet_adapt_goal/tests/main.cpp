#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This test makes sure that the benchmark "laplace_bc_dirichlet_adapt_goal" works correctly.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

const int NEQ = 1;
const int NELEM = 30;                     // Number of elements
const double A = -M_PI, B = M_PI;         // Domain end points
const double EPSILON = 0.01;              // Equation parameter
const int P_INIT = 1;                     // Initial polynomial degree

// JFNK or classical Newton.
const double MATRIX_SOLVER_TOL = 1e-5;
const int MATRIX_SOLVER_MAXITER = 150;
const int JFNK = 0;                       // 0... classical Newton
                                          // 1... JFNK  
                                          // Only relevant for JFNK:
const double JFNK_EPSILON = 1e-4;         // Parameter in the JFNK finite difference
 
// Newton's method.
double NEWTON_TOL_COARSE = 1e-8;          // Coarse mesh
double NEWTON_TOL_REF = 1e-8;             // Fine mesh.
int NEWTON_MAX_ITER = 150;

// Adaptivity.
const int GOAL_ORIENTED = 0;              // 0... standard adaptivity in norm.
                                          // 1... goal-oriented adaptivity.
const double X_QOI = 0.9*B;               // Value of u[0] at X_QOI is the 
                                          // Quantity of interest.
const double TOL_ERR_QOI = 1e-8;          // Tolerance for the maximum FTR error.
const int ADAPT_TYPE = 0;                 // 0... hp-adaptivity.
                                          // 1... h-adaptivity.
                                          // 2... p-adaptivity.
const double THRESHOLD = 0.7;             // Refined will be all elements whose FTR error 
                                          // is greater than THRESHOLD*max_ftr_error.
const int NORM = 1;                       // To measure errors.
                                          // 1... H1 norm.
                                          // 0... L2 norm.

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary conditions.
BCSpec DIR_BC_LEFT(0, 0);
BCSpec DIR_BC_RIGHT(0, 0);


// Function f(x).
double f(double x) 
{
  return sin(x);
}

// Exact solution.
// When changing exact solution, do not 
// forget to update interval accordingly.
const int EXACT_SOL_PROVIDED = 1;
void exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) 
{
  u[0] = sin(x);
  dudx[0] = cos(x);
}

// Weak forms for Jacobi matrix and residual.
#include "../definitions.cpp"

// Sample quantity of interest - solution value u[0] at 'x'.
// NOTE: quantity of interest is any linear functional of solution.
double quantity_of_interest(Space *space, double x)
{
  // Traversal of the space->
  Iterator *I = new Iterator(space);
  Element *e;
  double val[MAX_EQN_NUM];
  double der[MAX_EQN_NUM];
  while ((e = I->next_active_element()) != NULL) 
  {
    if (e->x1 <= x && x <= e->x2) 
    {
      e->get_solution_point(x, val, der);
      return val[0];
    }
  }

  error("computation of quantity of interest failed.");
  return 0.0;
}

int main() 
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Create coarse mesh, set Dirichlet BC, enumerate basis functions.
  Space* space = new Space(A, B, NELEM, DIR_BC_LEFT, DIR_BC_RIGHT, P_INIT, NEQ);

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
  if(JFNK == 0)
  {
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
      for(int i = 0; i < ndof_coarse; i++) rhs_coarse->set(i, -rhs_coarse->get(i));

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
  }
  else
    jfnk_cg(dp_coarse, space, MATRIX_SOLVER_TOL, MATRIX_SOLVER_MAXITER,
            JFNK_EPSILON, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, matrix_solver);

  // Cleanup.
  delete dp_coarse;

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;
  SimpleGraph graph_dof_exact, graph_cpu_exact;

  // Adaptivity loop:
  int as = 1;
  double ftr_errors[MAX_ELEM_NUM];        // This array decides what 
                                          // elements will be refined.

  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as); 

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = construct_refined_space(space);
 
    // Initialize the FE problem. 
    bool is_linear = false;
    DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space, is_linear);
      
    if(JFNK == 0)
    {
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
          for(int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));

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
      // Cleanup.
      delete matrix;
      delete rhs;
      delete solver;
      delete [] coeff_vec;
    }
    else
      jfnk_cg(dp, ref_space, MATRIX_SOLVER_TOL, MATRIX_SOLVER_MAXITER,
              JFNK_EPSILON, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, matrix_solver);
 
    // Cleanup.
    delete dp;
    
    // Starting with second adaptivity step, obtain new coarse 
    // mesh solution via projecting the fine mesh solution.
    if(as > 1)
    {
      info("Projecting the fine mesh solution onto the coarse mesh.");
      // Project the fine mesh solution (defined on space_ref) onto the coarse mesh (defined on space).
      OGProjection::project_global(space, ref_space, matrix_solver);
    }

    double max_qoi_err_est = 0;
    for (int i=0; i < space->get_n_active_elem(); i++)
    {
      if (GOAL_ORIENTED == 1) 
      {
        // Use quantity of interest.
        double qoi_est = quantity_of_interest(space, X_QOI);
        double qoi_ref_est = quantity_of_interest(ref_space, X_QOI);
        ftr_errors[i] = fabs(qoi_ref_est - qoi_est);
      }
      else 
      {
        // Use global norm
        double err_est_array[MAX_ELEM_NUM];
        ftr_errors[i] = calc_err_est(NORM, space, ref_space, err_est_array);
      }
      // Info for user.
      info("Elem [%d]: absolute error (est) = %g%%", i, ftr_errors[i]);

      // Time measurement.
      cpu_time.tick();

      // Calculating maximum of QOI FTR error for plotting purposes
      if (GOAL_ORIENTED == 1) 
      {
        if (ftr_errors[i] > max_qoi_err_est)
          max_qoi_err_est = ftr_errors[i];
      }
      else 
      {
        double qoi_est = quantity_of_interest(space, X_QOI);
        double qoi_ref_est = quantity_of_interest(ref_space, X_QOI);
        double err_est = fabs(qoi_ref_est - qoi_est);
        if (err_est > max_qoi_err_est)
          max_qoi_err_est = err_est;
      }
    }

    // Add entries to convergence graphs.
    if (EXACT_SOL_PROVIDED) 
    {
      double qoi_est = quantity_of_interest(space, X_QOI);
      double u[MAX_EQN_NUM], dudx[MAX_EQN_NUM];
      exact_sol(X_QOI, u, dudx);
      double err_qoi_exact = fabs(u[0] - qoi_est);
      // Info for user.
      info("Relative error (exact) = %g %%", err_qoi_exact);
      // Add entry to DOF and CPU convergence graphs.
      graph_dof_exact.add_values(Space::get_num_dofs(space), err_qoi_exact);
      graph_cpu_exact.add_values(cpu_time.accumulated(), err_qoi_exact);
    }
    
    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space::get_num_dofs(space), max_qoi_err_est);
    graph_cpu_est.add_values(cpu_time.accumulated(), max_qoi_err_est);

    // Decide whether the max. FTR error in the quantity of interest 
    // is sufficiently small.
    if(max_qoi_err_est < TOL_ERR_QOI) break;

    // Returns updated coarse and fine meshes, with the last 
    // coarse and fine mesh solutions on them, respectively. 
    // The coefficient vectors and numbers of degrees of freedom 
    // on both meshes are also updated. 
    adapt(NORM, ADAPT_TYPE, THRESHOLD, ftr_errors, space, ref_space);

    as++;

    // Plot meshes, results, and errors.
    adapt_plotting(space, ref_space, NORM, EXACT_SOL_PROVIDED, exact_sol);

    // Cleanup.
    delete ref_space;
  }
  while (done == false);

  info("Total running time: %g s", cpu_time.accumulated());

  // Save convergence graphs.
  graph_dof_est.save("conv_dof_est.dat");
  graph_cpu_est.save("conv_cpu_est.dat");
  graph_dof_exact.save("conv_dof_exact.dat");
  graph_cpu_exact.save("conv_cpu_exact.dat");

  // Test variable.
  bool success = true;
  info("ndof = %d.", Space::get_num_dofs(space));
  if (Space::get_num_dofs(space) > 150) success = false;

  if (success)
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
