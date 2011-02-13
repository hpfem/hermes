#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This test makes sure that the benchmark "first_order_general_adapt_ftr" works correctly.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

const int NEQ = 1;                      // Number of equations.
const int NELEM = 5;                    // Number of elements.
const double A = 0, B = 10;             // Domain end points.
const double YA = 1;                    // Equation parameter.
const int P_INIT = 1;                   // Initial polynomal degree.

// Newton's method.
const double NEWTON_TOL_COARSE = 1e-8;  // Coarse space
const double NEWTON_TOL_REF = 1e-8;     // Fine space
const int NEWTON_MAX_ITER = 150;

// Adaptivity.
const int ADAPT_TYPE = 0;         // 0... hp-adaptivity.
                                  // 1... h-adaptivity.
                                  // 2... p-adaptivity.
const double THRESHOLD = 0.7;     // Refined will be all elements whose error
                                  // is greater than THRESHOLD*max_elem_error.
const double TOL_ERR_FTR = 1e-2;  // Tolerance for the maximum FTR error.
const int NORM = 0;               // To measure errors.
                                  // 1... H1 norm.
                                  // 0... L2 norm.

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Boundary conditions.
Hermes::vector<BCSpec *>DIR_BC_LEFT =  Hermes::vector<BCSpec *>(new BCSpec(0, YA));
Hermes::vector<BCSpec *> DIR_BC_RIGHT = Hermes::vector<BCSpec *>();
 
// Right-hand side function f(y, x).
double f(double y, double x) 
{
  return -y*y;
}

// y-derivative of dfdy(y, x).
double dfdy(double y, double x) 
{
  return -2*y;
}

// Exact solution.
// When changing exact solution, do not 
// forget to update interval accordingly.
const int EXACT_SOL_PROVIDED = 1;
void exact_sol(double x, double u[MAX_EQN_NUM], double dudx[MAX_EQN_NUM]) 
{
  u[0] = 1./(x+1);
  dudx[0] = -1/((x+1)*(x+1));
}

// Weak forms for Jacobi matrix and residual.
#include "forms.cpp"

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

  double elem_errors[MAX_ELEM_NUM];      // This array decides what 
                                         // elements will be refined.
  ElemPtr2 ref_elem_pairs[MAX_ELEM_NUM]; // To store element pairs from the 
                                         // FTR solution. Decides how 
                                         // elements will be hp-refined. 
  for (int i=0; i < MAX_ELEM_NUM; i++) 
  {
    ref_elem_pairs[i][0] = new Element();
    ref_elem_pairs[i][1] = new Element();
  }

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_exact, graph_cpu_exact;

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

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

    // For every element perform its fast trial refinement (FTR),
    // calculate the norm of the difference between the FTR
    // solution and the coarse space solution, and store the
    // error in the elem_errors[] array.
    int n_elem = space->get_n_active_elem();
    for (int i=0; i < n_elem; i++) 
    {

      info("=== Starting FTR of Elem [%d].", i);

      // Replicate coarse space including solution.
      Space *space_ref_local = space->replicate();

      // Perform FTR of element 'i'
      space_ref_local->reference_refinement(i, 1);
      info("Elem [%d]: fine space created (%d DOF).", 
             i, space_ref_local->assign_dofs());

      // Initialize the FE problem. 
      bool is_linear = false;
      DiscreteProblem* dp = new DiscreteProblem(&wf, space_ref_local, is_linear);

      // Set up the solver, matrix, and rhs according to the solver selection.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Newton's loop on the FTR space.
      // Fill vector coeff_vec using dof and coeffs arrays in elements.
      double *coeff_vec = new double[Space::get_num_dofs(space_ref_local)];
      get_coeff_vector(space_ref_local, coeff_vec);
      memset(coeff_vec, 0, Space::get_num_dofs(space_ref_local)*sizeof(double));

      int it = 1;
      while (1) 
      {
        // Obtain the number of degrees of freedom.
        int ndof = Space::get_num_dofs(space_ref_local);

        // Assemble the Jacobian matrix and residual vector.
        dp->assemble(coeff_vec, matrix, rhs);

        // Calculate the l2-norm of residual vector.
        double res_l2_norm = get_l2_norm(rhs);

        // Info for user.
        info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(space_ref_local), res_l2_norm);

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
        set_coeff_vector(coeff_vec, space_ref_local);

        it++;
      }
      
      // Cleanup.
      delete matrix;
      delete rhs;
      delete solver;
      delete dp;
      delete [] coeff_vec;

      // Print FTR solution (enumerated). 
      //Linearizer *lxx = new Linearizer(space_ref_local);
      //char out_filename[255];
      //sprintf(out_filename, "solution_ref_%d.gp", i);
      //lxx->plot_solution(out_filename);
      //delete lxx;

      // Calculate norm of the difference between the coarse space 
      // and FTR solutions.
      // NOTE: later we want to look at the difference in some quantity 
      // of interest rather than error in global norm.
      double err_est_array[MAX_ELEM_NUM];
      elem_errors[i] = calc_err_est(NORM, space, space_ref_local, err_est_array) * 100;
      info("Elem [%d]: absolute error (est) = %g%%", i, elem_errors[i]);

      // Copy the reference element pair for element 'i'.
      // into the ref_elem_pairs[i][] array
      Iterator *I = new Iterator(space);
      Iterator *I_ref = new Iterator(space_ref_local);
      Element *e, *e_ref;
      while (1) 
      {
        e = I->next_active_element();
        e_ref = I_ref->next_active_element();
        if (e->id == i) 
        {
  	  e_ref->copy_into(ref_elem_pairs[e->id][0]);
          // coarse element 'e' was split in space.
          if (e->level != e_ref->level) 
          {
            e_ref = I_ref->next_active_element();
            e_ref->copy_into(ref_elem_pairs[e->id][1]);
          }
          break;
        }
      }

      delete I;
      delete I_ref;
      delete space_ref_local;
    }  

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

    // Calculate max FTR error.
    double max_ftr_error = 0;
    for (int i=0; i < space->get_n_active_elem(); i++) 
    {
      if (elem_errors[i] > max_ftr_error) max_ftr_error = elem_errors[i];
    }
    info("Max FTR error = %g%%.", max_ftr_error);

    // Decide whether the max. FTR error is sufficiently small.
    if(max_ftr_error < TOL_ERR_FTR) break;

    // debug
    //if (as == 4) break;

    // Returns updated coarse space with the last solution on it. 
    adapt(NORM, ADAPT_TYPE, THRESHOLD, elem_errors, space, ref_elem_pairs);

    // Plot spaces, results, and errors.
    //adapt_plotting(space, ref_elem_pairs, NORM, EXACT_SOL_PROVIDED, exact_sol);

    as++;
  }
  while (done == false);

  info("Total running time: %g s", cpu_time.accumulated());
  
  // Save convergence graphs.
  graph_dof_exact.save("conv_dof_exact.dat");
  graph_cpu_exact.save("conv_cpu_exact.dat");

  // Test variable.
  bool success = true;
  info("ndof = %d.", Space::get_num_dofs(space));
  if (Space::get_num_dofs(space) > 400) success = false;

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
