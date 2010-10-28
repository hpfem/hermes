#define HERMES_REPORT_ALL
#include "hermes1d.h"


// This example uses hp-adaptivity to solve the Poisson equation 
// -u'' - f = 0 in an interval (A, B), equipped with Dirichlet 
// boundary conditions on both end points. A series of small reference 
// solutions (we call them fast trial refinements, FTR) is used both 
// to decide what elements will be refined, and how they will be refined. 

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// General input.
static int NEQ = 1;
int NELEM = 3;                          // Number of elements.
double A = -M_PI, B = M_PI;             // Domain end points.
int P_init = 1;                         // Initial polynomal degree.

// Matrix solver.
const int MATRIX_SOLVER = 0;            // 0... default (LU decomposition).
                                        // 1... UMFPACK.
                                        // 2... CG (no preconditioning).
                                        // Only relevant for iterative matrix solvers:
const double MATRIX_SOLVER_TOL = 1e-7;  // Tolerance for residual in L2 norm.
const int MATRIX_SOLVER_MAXITER = 150;  // Max. number of iterations.

// Newton's method.
double NEWTON_TOL_COARSE = 1e-10;        // Coarse mesh.
double NEWTON_TOL_REF = 1e-10;           // Reference mesh.
int NEWTON_MAX_ITER = 150;

// Adaptivity.
const int ADAPT_TYPE = 0;               // 0... hp-adaptivity.
                                        // 1... h-adaptivity.
                                        // 2... p-adaptivity.
const double THRESHOLD = 0.7;           // Refined will be all elements whose FTR error 
                                        // is greater than THRESHOLD*max_ftr_error.
const double TOL_ERR_FTR = 1e-8;        // Tolerance for the maximum FTR error.
const int NORM = 1;                     // To measure errors.
                                        // 1... H1 norm.
                                        // 0... L2 norm.

// Boundary conditions.
double Val_dir_left = 0;
double Val_dir_right = 0;

// Function f(x).
double f(double x) {
  return sin(x);
}

// Exact solution.
// When changing exact solution, do not 
// forget to update interval accordingly.
const int EXACT_SOL_PROVIDED = 1;
void exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) {
  u[0] = sin(x);
  dudx[0] = cos(x);
}

// Weak forms for Jacobi matrix and residual.
#include "forms.cpp"

int main() {
  // Create coarse mesh, set Dirichlet BC, enumerate basis functions.
  Mesh *mesh = new Mesh(A, B, NELEM, P_init, NEQ);
  mesh->set_bc_left_dirichlet(0, Val_dir_left);
  mesh->set_bc_right_dirichlet(0, Val_dir_right);
  mesh->assign_dofs();

  // Initialize the FE problem.
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian);
  dp->add_vector_form(0, residual);

  // Convergence graph wrt. the number of degrees of freedom.
  GnuplotGraph graph;
  graph.set_log_y();
  graph.set_captions("Convergence History", "Degrees of Freedom", "Error");
  graph.add_row("exact error [%]", "k", "-", "o");
  graph.add_row("max FTR error", "k", "--");

  // Main adaptivity loop.
  int adapt_iterations = 1;
  double ftr_errors[MAX_ELEM_NUM];        // This array decides what 
                                          // elements will be refined.
  ElemPtr2 ref_ftr_pairs[MAX_ELEM_NUM];   // To store element pairs from the 
                                          // FTR solution. Decides how 
                                          // elements will be hp-refined. 
  for (int i=0; i < MAX_ELEM_NUM; i++) {
    ref_ftr_pairs[i][0] = new Element();
    ref_ftr_pairs[i][1] = new Element();
  }
  while(1) {
    info("============ Adaptivity step %d ============", adapt_iterations); 

    info("N_dof = %d", mesh->get_num_dofs());
 
    // Newton's loop on coarse mesh.
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
      
      // Copy coefficients from vector y to elements.
      vector_to_solution(y, mesh);
    }
    
    delete matrix;
    delete rhs;
    delete solver;

    // For every element perform its fast trial refinement (FTR),
    // calculate the norm of the difference between the FTR
    // solution and the coarse mesh solution, and store the
    // error in the ftr_errors[] array.
    int n_elem = mesh->get_n_active_elem();
    for (int i=0; i < n_elem; i++) {

      info("=== Starting FTR of Elem [%d]", i);

      // Replicate coarse mesh including solution.
      Mesh *mesh_ref_local = mesh->replicate();

      // Perform FTR of element 'i'
      mesh_ref_local->reference_refinement(i, 1);
      info("Elem [%d]: fine mesh created (%d DOF).", 
             i, mesh_ref_local->assign_dofs());

      // Obtain the number of degrees of freedom.
      int ndof = mesh_ref_local->get_num_dofs();

      // Fill vector y using dof and coeffs arrays in elements.
      double *y = new double[ndof];
      solution_to_vector(mesh_ref_local, y);
    
      // Set up the solver, matrix, and rhs according to the solver selection.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
    
      int it = 1;
      while (1)
      {
        // Construct matrix and residual vector.
        dp->assemble_matrix_and_vector(mesh_ref_local, matrix, rhs);

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
        
        // Copy coefficients from vector y to elements.
        vector_to_solution(y, mesh_ref_local);
      }
      
      delete matrix;
      delete rhs;
      delete solver;

      // Print FTR solution (enumerated).
      Linearizer *lxx = new Linearizer(mesh_ref_local);
      char out_filename[255];
      sprintf(out_filename, "solution_ref_%d.gp", i);
      lxx->plot_solution(out_filename);
      delete lxx;

      // Calculate norm of the difference between the coarse mesh 
      // and FTR solutions.
      // NOTE: later we want to look at the difference in some quantity 
      // of interest rather than error in global norm.
      double err_est_array[MAX_ELEM_NUM];
      ftr_errors[i] = calc_error_estimate(NORM, mesh, mesh_ref_local, 
                      err_est_array);

      // Copy the reference element pair for element 'i'
      // into the ref_ftr_pairs[i][] array.
      Iterator *I = new Iterator(mesh);
      Iterator *I_ref = new Iterator(mesh_ref_local);
      Element *e, *e_ref;
      while (1) {
        e = I->next_active_element();
        e_ref = I_ref->next_active_element();
        if (e->id == i) {
  	  e_ref->copy_into(ref_ftr_pairs[e->id][0]);
          // Coarse element 'e' was split in space.
          if (e->level != e_ref->level) {
            e_ref = I_ref->next_active_element();
            e_ref->copy_into(ref_ftr_pairs[e->id][1]);
          }
          break;
        }
      }

      delete I;
      delete I_ref;
      delete mesh_ref_local;
    }  

    // If exact solution available, also calculate exact error.
    if (EXACT_SOL_PROVIDED) {
      // Calculate element errors wrt. exact solution.
      double err_exact_total = calc_error_exact(NORM, mesh, exact_sol);
     
      // Calculate the norm of the exact solution
      // (using a fine subdivision and high-order quadrature).
      int subdivision = 500; // heuristic parameter
      int order = 20;        // heuristic parameter
      double exact_sol_norm = calc_solution_norm(NORM, exact_sol, NEQ, A, B,
                                                  subdivision, order);
      // Calculate an estimate of the global relative error.
      double err_exact_rel = err_exact_total/exact_sol_norm;
      graph.add_values(0, mesh->get_num_dofs(), 100 * err_exact_rel);
    }

    // Calculate max FTR error.
    double max_ftr_error = 0;
    for (int i=0; i < mesh->get_n_active_elem(); i++) {
      if (ftr_errors[i] > max_ftr_error) max_ftr_error = ftr_errors[i];
    }
    info("Max FTR error = %g", max_ftr_error);

    // Add entry to DOF convergence graph.
    graph.add_values(1, mesh->get_num_dofs(), max_ftr_error);

    // Decide whether the max. FTR error is sufficiently small.
    if(max_ftr_error < TOL_ERR_FTR) break;

    // Returns updated coarse mesh with the last solution on it. 
    adapt(NORM, ADAPT_TYPE, THRESHOLD, ftr_errors,
          mesh, ref_ftr_pairs);

    adapt_iterations++;
  }

  // Plot meshes, results, and errors.
  adapt_plotting(mesh, ref_ftr_pairs,
                 NORM, EXACT_SOL_PROVIDED, exact_sol);

  // Save convergence graph.
  graph.save("conv_dof.gp");

  info("Done.");
  return 1;
}









