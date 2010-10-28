#define HERMES_REPORT_ALL
#include "hermes1d.h"

// This example solves a system of two linear second-order equations 
// u' + k^2 v = 0
// u - v' = 0
// which is equivalent to u'' + k^2 u = 0
// in an interval (0, 2*pi) equipped with Dirichlet bdy conditions 
// u(0) = 0, v(0) = k
// The exact solution is u(x) = sin(k*x), v(x) = k*cos(k*x)

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// General input.
static int NEQ = 2;
int NELEM = 3;                          // Number of elements.
double A = 0, B = 2*M_PI;               // Domain end points.
int P_init = 1;                         // Initial polynomial degree.
double K = 1.0;                         // Equation parameter.

// Newton's method.
double NEWTON_TOL_COARSE = 1e-6;        // Coarse mesh.
double NEWTON_TOL_REF = 1e-6;           // Fine mesh.
int NEWTON_MAX_ITER = 150;

// Adaptivity.
const int ADAPT_TYPE = 0;         // 0... hp-adaptivity.
                                  // 1... h-adaptivity.
                                  // 2... p-adaptivity.
const double THRESHOLD = 0.7;     // Refined will be all elements whose error 
                                  // is greater than THRESHOLD*max_elem_error.
const double TOL_ERR_REL = 1e-1;  // Tolerance for the relative error between 
                                  // the coarse mesh and fine solutions.
const int NORM = 1;               // To measure errors.
                                  // 1... H1 norm.
                                  // 0... L2 norm.

// Boundary conditions.
double Val_dir_left_0 = 0;
double Val_dir_left_1 = K;

// Exact solution.
// When changing exact solution, do not 
// forget to update interval accordingly.
const int EXACT_SOL_PROVIDED = 1;
void exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) {
  u[0] = sin(K*x);
  dudx[0] = K*cos(K*x);
  u[1] = K*cos(K*x);
  dudx[1] = -K*K*sin(K*x);
}

// Weak forms for Jacobi matrix and residual.
#include "forms.cpp"


int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate basis functions.
  Space *space = new Space(A, B, NELEM, P_init, NEQ);
  space->set_bc_left_dirichlet(0, Val_dir_left_0);
  space->set_bc_left_dirichlet(1, Val_dir_left_1);
  info("N_dof = %d", space->assign_dofs());

  // Initialize the FE problem.
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_0_0);
  dp->add_matrix_form(0, 1, jacobian_0_1);
  dp->add_matrix_form(1, 0, jacobian_1_0);
  dp->add_matrix_form(1, 1, jacobian_1_1);
  dp->add_vector_form(0, residual_0);
  dp->add_vector_form(1, residual_1);

  // Newton's loop on coarse mesh.
  // Obtain the number of degrees of freedom.
  int ndof = Space::get_num_dofs(space);

  // Fill vector y using dof and coeffs arrays in elements.
  double *coeff_vec = new double[ndof];
  solution_to_vector(space, coeff_vec);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  int it = 1;
  while (1)
  {
    // Assemble the Jacobian matrix and residual vector.
    dp->assemble_matrix_and_vector(space, matrix, rhs);

    // Calculate the l2-norm of residual vector.
    double res_norm_squared = 0;
    for(int i=0; i<ndof; i++) res_norm_squared += rhs->get(i)*rhs->get(i);

    // Info for user.
    info("---- Newton iter %d, residual norm: %.15f", it, sqrt(res_norm_squared));

    // If l2 norm of the residual vector is within tolerance, then quit.
    // NOTE: at least one full iteration forced
    //       here because sometimes the initial
    //       residual on fine mesh is too small.
    if(res_norm_squared < NEWTON_TOL_COARSE*NEWTON_TOL_COARSE && it > 1) break;

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
    vector_to_solution(coeff_vec, space);
    
    it++;
  }
  
  // Cleanup.
  delete matrix;
  delete rhs;
  delete solver;

  // Replicate coarse mesh including solution.
  Space *space_ref = space->replicate();

  // Refine entire space_ref uniformly in 'h' and 'p'.
  int start_elem_id = 0; 
  int num_to_ref = space_ref->get_n_active_elem();
  space_ref->reference_refinement(start_elem_id, num_to_ref);
  info("Fine mesh created (%d DOF).", Space::get_num_dofs(space_ref));

  // Convergence graph wrt. the number of degrees of freedom.
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Convergence History", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");

  // Main adaptivity loop.
  int adapt_iterations = 1;
  while(1) {
    info("============ Adaptivity step %d ============", adapt_iterations); 

    // Newton's loop on fine mesh.
    // Obtain the number of degrees of freedom.
    int ndof = Space::get_num_dofs(space_ref);

    // Fill vector y using dof and coeffs arrays in elements.
    double *coeff_vec = new double[ndof];
    solution_to_vector(space_ref, coeff_vec);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    int it = 1;
    while (1)
    {
      // Assemble the Jacobian matrix and residual vector.
      dp->assemble_matrix_and_vector(space_ref, matrix, rhs);

      // Calculate the l2-norm of residual vector.
      double res_norm_squared = 0;
      for(int i=0; i<ndof; i++) res_norm_squared += rhs->get(i)*rhs->get(i);

      // Info for user.
      info("---- Newton iter %d, residual norm: %.15f", it, sqrt(res_norm_squared));

      // If l2 norm of the residual vector is within tolerance, then quit.
      // NOTE: at least one full iteration forced
      //       here because sometimes the initial
      //       residual on fine mesh is too small.
      if(res_norm_squared < NEWTON_TOL_REF*NEWTON_TOL_REF && it > 1) break;

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
      vector_to_solution(coeff_vec, space_ref);

      it++;
    }
    
    // Cleanup.
    delete matrix;
    delete rhs;
    delete solver;

    // Starting with second adaptivity step, obtain new coarse 
    // space solution via Newton's method. Initial condition is 
    // the last coarse mesh solution.
    if (adapt_iterations > 1) {

      // Newton's loop on coarse mesh.
      // Obtain the number of degrees of freedom.
      int ndof = Space::get_num_dofs(space);

      // Fill vector y using dof and coeffs arrays in elements.
      double *coeff_vec = new double[ndof];
      solution_to_vector(space, coeff_vec);

      // Set up the solver, matrix, and rhs according to the solver selection.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      int it = 1;
      while (1)
      {
        // Assemble the Jacobian matrix and residual vector.
        dp->assemble_matrix_and_vector(space, matrix, rhs);

        // Calculate the l2-norm of residual vector.
        double res_norm_squared = 0;
        for(int i=0; i<ndof; i++) res_norm_squared += rhs->get(i)*rhs->get(i);

        // Info for user.
        info("---- Newton iter %d, residual norm: %.15f", it, sqrt(res_norm_squared));

        // If l2 norm of the residual vector is within tolerance, then quit.
        // NOTE: at least one full iteration forced
        //       here because sometimes the initial
        //       residual on fine mesh is too small.
        if(res_norm_squared < NEWTON_TOL_COARSE*NEWTON_TOL_COARSE && it > 1) break;

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
        vector_to_solution(coeff_vec, space);

        it++;
      }
      
      // Cleanup.
      delete matrix;
      delete rhs;
      delete solver;
    }

    // In the next step, estimate element errors based on 
    // the difference between the fine mesh and coarse mesh solutions. 
    double err_est_array[MAX_ELEM_NUM]; 
    double err_est_total = calc_error_estimate(NORM, 
              space, space_ref, err_est_array);

    // Calculate the norm of the fine mesh solution.
    double ref_sol_norm = calc_solution_norm(NORM, space_ref);

    // Calculate an estimate of the global relative error.
    double err_est_rel = err_est_total/ref_sol_norm;
    info("Relative error (est) = %g %%", 100.*err_est_rel);

    // If exact solution available, also calculate exact error.
    if (EXACT_SOL_PROVIDED) {
      // Calculate element errors wrt. exact solution.
      double err_exact_total = calc_error_exact(NORM, 
         space, exact_sol);
     
      // Calculate the norm of the exact solution
      // (using a fine subdivision and high-order quadrature).
      int subdivision = 500; // heuristic parameter
      int order = 20;        // heuristic parameter
      double exact_sol_norm = calc_solution_norm(NORM, 
           exact_sol, NEQ, A, B, subdivision, order);
      // Calculate an estimate of the global relative error.
      double err_exact_rel = err_exact_total/exact_sol_norm;
      info("Relative error (exact) = %g %%", 100.*err_exact_rel);
      graph.add_values(0, Space::get_num_dofs(space), 100 * err_exact_rel);
    }

    // Add entry to DOF convergence graph.
    graph.add_values(1, Space::get_num_dofs(space), 100 * err_est_rel);

    // Decide whether the relative error is sufficiently small.
    if(err_est_rel*100 < TOL_ERR_REL) break;

    // Returns updated coarse and fine meshes, with the last 
    // coarse and fine mesh solutions on them, respectively. 
    // The coefficient vectors and numbers of degrees of freedom 
    // on both meshes are also updated. 
    adapt(NORM, ADAPT_TYPE, THRESHOLD, err_est_array,
          space, space_ref);

    adapt_iterations++;
  }

  // Plot meshes, results, and errors.
  adapt_plotting(space, space_ref, 
                 NORM, EXACT_SOL_PROVIDED, exact_sol);

  // Save convergence graph.
  graph.save("conv_dof.gp");

  info("Done.");
  return 1;
}
