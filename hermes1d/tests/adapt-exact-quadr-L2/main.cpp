#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This test makes sure that an exact function 
// 1-x^2 is found after one step of hp-adaptivity, 
// and that both refinements were p-refinements.
// Adaptivity is done in L2 norm.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

// General input:
static int NEQ = 1;
int NELEM = 2;                         // Number of elements
double A = -1, B = 1;               // Domain end points
int P_init = 1;                         // Initial polynomal degree

// Newton's method
double NEWTON_TOL_COARSE = 1e-6;        // Coarse mesh
double NEWTON_TOL_REF = 1e-6;           // Reference mesh
int NEWTON_MAX_ITER = 150;

// Adaptivity
const int ADAPT_TYPE = 0;               // 0... hp-adaptivity
                                        // 1... h-adaptivity
                                        // 2... p-adaptivity
const double THRESHOLD = 0.7;           // Refined will be all elements whose error
                                        // is greater than THRESHOLD*max_elem_error
const double TOL_ERR_REL = 1e-3;        // Tolerance for the relative error between 
                                        // the coarse mesh and reference solutions
const int NORM = 0;                     // 1... H1 norm
                                        // 0... L2 norm

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Boundary conditions
double Val_dir_left = 0;                // Dirichlet condition left
double Val_dir_right = 0;               // Dirichlet condition right

// Function f(x)
double f(double x) {
  //return sin(x);
  return 2;
}

// Exact solution
const int EXACT_SOL_PROVIDED = 1;
void exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) {
  //u[0] = sin(x);
  //dudx[0] = cos(x);
  u[0] = 1. - x*x;
  dudx[0] = -2.*x;
}

// ********************************************************************

// bilinear form for the Jacobi matrix 
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...basis function
// v...test function
// u_prev...previous solution (all solution components)
double jacobian(int num, double *x, double *weights, 
                double *u, double *dudx, double *v, double *dvdx, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                void *user_data)
{
  double val = 0;
  for(int i = 0; i<num; i++) {
    val += dudx[i]*dvdx[i]*weights[i];
  }
  return val;
};

// (nonlinear) form for the residual vector
// num...number of Gauss points in element
// x[]...Gauss points
// weights[]...Gauss weights for points in x[]
// u...approximate solution
// v...test function
// u_prev...previous solution (all solution components)
double residual(int num, double *x, double *weights, 
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM], 
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],  
                double *v, double *dvdx, void *user_data)
{
  double val = 0;
  int si = 0;      // solution index (only 0 is relevant for this example)
  for(int i = 0; i<num; i++) {
    val += (du_prevdx[si][0][i]*dvdx[i] - f(x[i])*v[i])*weights[i];
  }
  return val;
};

/******************************************************************************/
int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate basis functions
  Mesh *mesh = new Mesh(A, B, NELEM, P_init, NEQ);
  mesh->set_bc_left_dirichlet(0, Val_dir_left);
  mesh->set_bc_right_dirichlet(0, Val_dir_right);
  info("N_dof = %d\n", mesh->assign_dofs());

  // Create discrete problem on coarse mesh
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian);
  dp->add_vector_form(0, residual);

  // Newton's loop on coarse mesh
  // Obtain the number of degrees of freedom.
  int ndof = mesh->get_num_dofs();

  // Fill vector y using dof and coeffs arrays in elements.
  double *y = new double[ndof];
  solution_to_vector(mesh, y);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  int it = 1;
  while (1)
  {
    // Construct matrix and residual vector.
    dp->assemble_matrix_and_vector(mesh, matrix, rhs);

    // Calculate L2 norm of residual vector.
    double res_norm_squared = 0;
    for(int i=0; i<ndof; i++) res_norm_squared += rhs->get(i)*rhs->get(i);

    info("---- Newton iter %d, residual norm: %.15f\n", it, sqrt(res_norm_squared));

    // If residual norm less than 'NEWTON_TOL', quit
    // latest solution is in the vector y.
    // NOTE: at least one full iteration forced
    //       here because sometimes the initial
    //       residual on fine mesh is too small
    if(res_norm_squared < NEWTON_TOL_COARSE*NEWTON_TOL_COARSE && it > 1) break;

    // Changing sign of vector res.
    for(int i=0; i<ndof; i++) rhs->set(i, -rhs->get(i));

    // Calculate the coefficient vector.
    bool solved = solver->solve();
    if (solved) 
    {
      double* solution_vector = new double[ndof];
      solution_vector = solver->get_solution();
      for(int i=0; i<ndof; i++) y[i] += solution_vector[i];
      // No need to deallocate the solution_vector here, it is done later by the call to ~Solver.
      solution_vector = NULL;
    }
    it++;

    if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
    
    // copy coefficients from vector y to elements
    vector_to_solution(y, mesh);
  }
  
  delete matrix;
  delete rhs;
  delete solver;

  // Replicate coarse mesh including dof arrays
  Mesh *mesh_ref = mesh->replicate();

  // Refine entire mesh_ref uniformly in 'h' and 'p'
  int start_elem_id = 0; 
  int num_to_ref = mesh_ref->get_n_active_elem();
  mesh_ref->reference_refinement(start_elem_id, num_to_ref);
  info("Fine mesh created (%d DOF).\n", mesh_ref->get_num_dofs());

  // Convergence graph wrt. the number of degrees of freedom
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Convergence History", "Degrees of Freedom", "Error [%]");
  graph.add_row("exact error", "k", "-", "o");
  graph.add_row("error estimate", "k", "--");

  // Main adaptivity loop
  int adapt_iterations = 1;
  int success_test = 1;
  while(1) {
    info("============ Adaptivity step %d ============\n", adapt_iterations); 

    // Newton's loop on fine mesh
    // Obtain the number of degrees of freedom.
    int ndof = mesh_ref->get_num_dofs();

    // Fill vector y using dof and coeffs arrays in elements.
    double *y = new double[ndof];
    solution_to_vector(mesh_ref, y);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    int it = 1;
    while (1)
    {
      // Construct matrix and residual vector.
      dp->assemble_matrix_and_vector(mesh_ref, matrix, rhs);

      // Calculate L2 norm of residual vector.
      double res_norm_squared = 0;
      for(int i=0; i<ndof; i++) res_norm_squared += rhs->get(i)*rhs->get(i);

      info("---- Newton iter %d, residual norm: %.15f\n", it, sqrt(res_norm_squared));

      // If residual norm less than 'NEWTON_TOL', quit
      // latest solution is in the vector y.
      // NOTE: at least one full iteration forced
      //       here because sometimes the initial
      //       residual on fine mesh is too small
      if(res_norm_squared < NEWTON_TOL_REF*NEWTON_TOL_REF && it > 1) break;

      // Changing sign of vector res.
      for(int i=0; i<ndof; i++) rhs->set(i, -rhs->get(i));

      // Calculate the coefficient vector.
      bool solved = solver->solve();
      if (solved) 
      {
        double* solution_vector = new double[ndof];
        solution_vector = solver->get_solution();
        for(int i=0; i<ndof; i++) y[i] += solution_vector[i];
        // No need to deallocate the solution_vector here, it is done later by the call to ~Solver.
        solution_vector = NULL;
      }
      it++;

      if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
      
      // copy coefficients from vector y to elements
      vector_to_solution(y, mesh_ref);
    }
    
    delete matrix;
    delete rhs;
    delete solver;

    // Starting with second adaptivity step, obtain new coarse 
    // mesh solution via Newton's method. Initial condition is 
    // the last coarse mesh solution.
    if (adapt_iterations > 1) {

      // Newton's loop on coarse mesh
      // Obtain the number of degrees of freedom.
      int ndof = mesh->get_num_dofs();

      // Fill vector y using dof and coeffs arrays in elements.
      double *y = new double[ndof];
      solution_to_vector(mesh, y);

      // Set up the solver, matrix, and rhs according to the solver selection.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      int it = 1;
      while (1)
      {
        // Construct matrix and residual vector.
        dp->assemble_matrix_and_vector(mesh, matrix, rhs);

        // Calculate L2 norm of residual vector.
        double res_norm_squared = 0;
        for(int i=0; i<ndof; i++) res_norm_squared += rhs->get(i)*rhs->get(i);

        info("---- Newton iter %d, residual norm: %.15f\n", it, sqrt(res_norm_squared));

        // If residual norm less than 'NEWTON_TOL', quit
        // latest solution is in the vector y.
        // NOTE: at least one full iteration forced
        //       here because sometimes the initial
        //       residual on fine mesh is too small
        if(res_norm_squared < NEWTON_TOL_COARSE*NEWTON_TOL_COARSE && it > 1) break;

        // Changing sign of vector res.
        for(int i=0; i<ndof; i++) rhs->set(i, -rhs->get(i));

        // Calculate the coefficient vector.
        bool solved = solver->solve();
        if (solved) 
        {
          double* solution_vector = new double[ndof];
          solution_vector = solver->get_solution();
          for(int i=0; i<ndof; i++) y[i] += solution_vector[i];
          // No need to deallocate the solution_vector here, it is done later by the call to ~Solver.
          solution_vector = NULL;
        }
        it++;

        if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
        
        // copy coefficients from vector y to elements
        vector_to_solution(y, mesh);
      }
      
      delete matrix;
      delete rhs;
      delete solver;
    }

    // In the next step, estimate element errors based on 
    // the difference between the fine mesh and coarse mesh solutions. 
    double err_est_array[MAX_ELEM_NUM]; 
    double err_est_total = calc_error_estimate(NORM, mesh, mesh_ref, 
                           err_est_array);

    // Calculate the norm of the fine mesh solution
    double ref_sol_norm = calc_solution_norm(NORM, mesh_ref);

    // Calculate an estimate of the global relative error
    double err_est_rel = err_est_total/ref_sol_norm;
    info("Relative error (est) = %g %%\n", 100.*err_est_rel);

    // If exact solution available, also calculate exact error
    double err_exact_rel;  
    if (EXACT_SOL_PROVIDED) {
      // Calculate element errors wrt. exact solution
      double err_exact_total = calc_error_exact(NORM, mesh, exact_sol);
     
      // Calculate the norm of the exact solution
      // (using a fine subdivision and high-order quadrature)
      int subdivision = 500; // heuristic parameter
      int order = 20;        // heuristic parameter
      double exact_sol_norm = calc_solution_norm(NORM, exact_sol, NEQ, A, B,
                                                  subdivision, order);
      // Calculate an estimate of the global relative error
      err_exact_rel = err_exact_total/exact_sol_norm;
      info("Relative error (exact) = %g %%\n", 100.*err_exact_rel);
      graph.add_values(0, mesh->get_num_dofs(), 100 * err_exact_rel);
    }

    // add entry to DOF convergence graph
    graph.add_values(1, mesh->get_num_dofs(), 100 * err_est_rel);

    // Decide whether the relative error is sufficiently small
    //if(err_est_rel*100 < TOL_ERR_REL) break;

    // extra code for this test:
    if (adapt_iterations == 2) {
      if (err_est_rel > 1e-10) success_test = 0;
      if (err_exact_rel > 1e-10) success_test = 0;
      if (mesh->get_n_active_elem() != 2) success_test = 0;
      Element *e = mesh->first_active_element();
      if (e->p != 2) success_test = 0;
      e = mesh->last_active_element();
      if (e->p != 2) success_test = 0;
      break;
    }

    // Returns updated coarse and fine meshes, with the last 
    // coarse and fine mesh solutions on them, respectively. 
    // The coefficient vectors and numbers of degrees of freedom 
    // on both meshes are also updated. 
    adapt(NORM, ADAPT_TYPE, THRESHOLD, err_est_array,
          mesh, mesh_ref);

    adapt_iterations++;
  }

  // Plot meshes, results, and errors
  adapt_plotting(mesh, mesh_ref, 
                 NORM, EXACT_SOL_PROVIDED, exact_sol);

  // Save convergence graph
  graph.save("conv_dof.gp");

  if (success_test) {
    info("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    info("Failure!\n");
    return ERROR_FAILURE;
  }
}

