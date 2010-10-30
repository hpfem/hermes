#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This example uses automatic hp-adaptivity to solve a general 
// first-order equation.
//
//  PDE: y' = f(y, x).
//
//  Interval: (A, B).
//
//  BC: Dirichlet, y(A) = YA.
//
//  Exact solution: y(x) = sin(x).
//
//  The following parameters can be changed:
const int NEQ = 1;                      // Number of equations.
const int NELEM = 5;                    // Number of elements.
const double A = 0, B = 10;             // Domain end points.
const double YA = 1;                    // Equation parameter.
const int P_init = 1;                   // Initial polynomial degree.

// Newton's method.
double NEWTON_TOL_COARSE = 1e-8;        // Coarse mesh.
double NEWTON_TOL_REF = 1e-8;           // Fine mesh.
int NEWTON_MAX_ITER = 150;

// Adaptivity.
const int ADAPT_TYPE = 0;         // 0... hp-adaptivity.
                                  // 1... h-adaptivity.
                                  // 2... p-adaptivity.
const double THRESHOLD = 0.7;     // Refined will be all elements whose error 
                                  // is greater than THRESHOLD*max_elem_error.
const double TOL_ERR_REL = 1e-5;  // Tolerance for the relative error between 
                                  // the coarse and fine mesh solutions.
const int NORM = 0;               // To measure errors.
                                  // 1... H1 norm.
                                  // 0... L2 norm.
 
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Right-hand side function f(y, x).
double f(double y, double x) {
  // With y(0)=1, exact solution is y=1/(x+1).
  return -y*y; 
}

// y-derivative of dfdy(y, x).
double dfdy(double y, double x) {
  return -2*y;
}

// Exact solution.
const int EXACT_SOL_PROVIDED = 1;
void exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) {
  u[0] = 1./(x+1);
  dudx[0] = -1/((x+1)*(x+1));
}

// Weak forms for Jacobi matrix and residual.
#include "forms.cpp"


int main() {
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Create coarse mesh, set Dirichlet BC, enumerate basis functions.
  Space* space = new Space(A, B, NELEM, P_init, NEQ);
  space->set_bc_left_dirichlet(0, YA);
  info("N_dof = %d", space->assign_dofs());

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(jacobian);
  wf.add_vector_form(residual);

  // Initialize the FE problem.
  DiscreteProblem *dp_coarse = new DiscreteProblem(&wf, space);

  // Newton's loop on coarse mesh.
  // Fill vector coeff_vec using dof and coeffs arrays in elements.
  double *coeff_vec_coarse = new double[Space::get_num_dofs(space)];
  solution_to_vector(space, coeff_vec_coarse);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
  Vector* rhs_coarse = create_vector(matrix_solver);
  Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

  int it = 1;
  while (1) {
    // Obtain the number of degrees of freedom.
    int ndof_coarse = Space::get_num_dofs(space);

    // Assemble the Jacobian matrix and residual vector.
    dp_coarse->assemble(matrix_coarse, rhs_coarse);

    // Calculate the l2-norm of residual vector.
    double res_norm_squared = 0;
    for(int i=0; i<ndof_coarse; i++) res_norm_squared += rhs_coarse->get(i)*rhs_coarse->get(i);

    // Info for user.
    info("---- Newton iter %d, residual norm: %.15f", it, sqrt(res_norm_squared));

    // If l2 norm of the residual vector is within tolerance, then quit.
    // NOTE: at least one full iteration forced
    //       here because sometimes the initial
    //       residual on fine mesh is too small.
    if(res_norm_squared < NEWTON_TOL_COARSE*NEWTON_TOL_COARSE && it > 1) break;

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
    vector_to_solution(coeff_vec_coarse, space);

    it++;
  }
  
  // Cleanup.
  delete matrix_coarse;
  delete rhs_coarse;
  delete solver_coarse;

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;
  SimpleGraph graph_dof_exact, graph_cpu_exact;

  // Main adaptivity loop.
  int as = 1;
  while(1) {
    info("============ Adaptivity step %d ============", as); 

    // Construct globally refined reference mesh and setup reference space->
    Space* ref_space = construct_refined_space(space);

    // Initialize the FE problem. 
    DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Newton's loop on fine mesh.
    // Fill vector coeff_vec using dof and coeffs arrays in elements.
    double *coeff_vec = new double[Space::get_num_dofs(ref_space)];
    solution_to_vector(ref_space, coeff_vec);

    int it = 1;
    while (1) {
      // Obtain the number of degrees of freedom.
      int ndof = Space::get_num_dofs(ref_space);

      // Assemble the Jacobian matrix and residual vector.
      dp->assemble(matrix, rhs);

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
      vector_to_solution(coeff_vec, ref_space);

      it++;
    }
    
    // Cleanup.
    delete matrix;
    delete rhs;
    delete solver;

    // Starting with second adaptivity step, obtain new coarse 
    // space solution via Newton's method. Initial condition is 
    // the last coarse mesh solution.
    if (as > 1) {
      //Info for user.
      info("Solving on coarse mesh");

      // Initialize the FE problem.
      DiscreteProblem* dp_coarse = new DiscreteProblem(&wf, space);

      // Newton's loop on coarse mesh.
      // Fill vector coeff_vec using dof and coeffs arrays in elements.
      double *coeff_vec_coarse = new double[Space::get_num_dofs(space)];
      solution_to_vector(space, coeff_vec_coarse);

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
        dp_coarse->assemble(matrix_coarse, rhs_coarse);

        // Calculate the l2-norm of residual vector.
        double res_norm_squared = 0;
        for(int i=0; i<ndof_coarse; i++) res_norm_squared += rhs_coarse->get(i)*rhs_coarse->get(i);

        // Info for user.
        info("---- Newton iter %d, residual norm: %.15f", it, sqrt(res_norm_squared));

        // If l2 norm of the residual vector is within tolerance, then quit.
        // NOTE: at least one full iteration forced
        //       here because sometimes the initial
        //       residual on fine mesh is too small.
        if(res_norm_squared < NEWTON_TOL_COARSE*NEWTON_TOL_COARSE && it > 1) break;

        // Multiply the residual vector with -1 since the matrix 
        // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
        for(int i=0; i<ndof_coarse; i++) rhs_coarse->set(i, -rhs_coarse->get(i));

        // Solve the linear system.
        if(!solver->solve())
          error ("Matrix solver failed.\n");

        // Add \deltaY^{n+1} to Y^n.
        for (int i = 0; i < ndof_coarse; i++) coeff_vec[i] += solver->get_solution()[i];

        // If the maximum number of iteration has been reached, then quit.
        if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
        
        // Copy coefficients from vector y to elements.
        vector_to_solution(coeff_vec_coarse, space);

        it++;
      }
      
      // Cleanup.
      delete matrix_coarse;
      delete rhs_coarse;
      delete solver_coarse;
      delete dp;
      delete [] coeff_vec_coarse;
    }

    // In the next step, estimate element errors based on 
    // the difference between the fine mesh and coarse mesh solutions. 
    double err_est_array[MAX_ELEM_NUM]; 
    double err_est_rel = calc_err_est(NORM, 
              space, ref_space, err_est_array) * 100;

    // Info for user.
    info("Relative error (est) = %g %%", err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // If exact solution available, also calculate exact error.
    if (EXACT_SOL_PROVIDED) {
      // Calculate element errors wrt. exact solution.
      double err_exact_rel = calc_err_exact(NORM, 
         space, exact_sol, NEQ, A, B) * 100;
     
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
    adapt(NORM, ADAPT_TYPE, THRESHOLD, err_est_array,
          space, ref_space);

    as++;

  // Plot meshes, results, and errors.
    adapt_plotting(space, ref_space, 
           NORM, EXACT_SOL_PROVIDED, exact_sol);

    // Cleanup.
    delete ref_space;
    delete dp;
    delete [] coeff_vec;

  }

  // Save convergence graphs.
  graph_dof_est.save("conv_dof_est.dat");
  graph_cpu_est.save("conv_cpu_est.dat");
  graph_dof_exact.save("conv_dof_exact.dat");
  graph_cpu_exact.save("conv_cpu_exact.dat");

  info("Done.");
  return 1;
}
