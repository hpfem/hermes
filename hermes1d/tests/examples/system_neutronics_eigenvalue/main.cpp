#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This test makes sure that example "system_neutronics_eigenvalue" works correctly.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

// Problem specification (core geometry, material properties, initial FE space).
#include "neutronics_problem_def.cpp"

// Common functions for neutronics problems (requires variable declarations from
// "neutronics_problem_def.cpp").
#include "neutronics_common.cpp"

// Weak forms for the problem (requires variable declarations from
// "neutronics_problem_def.cpp").
#include "forms.cpp"

// General input (eigenvalue problem).

bool flag = false;					// Flag for debugging purposes.
bool verbose = true;

// Power method initialization.
int Max_SI = 1000;          // Max. number of eigenvalue iterations.
int N_SLN = 2;              // Number of solutions.

// Newton's method.
double NEWTON_TOL = 1e-5;               // Tolerance.
int NEWTON_MAX_ITER = 150;              // Max. number of Newton iterations.
double TOL_SI = 1e-6;                   // Tolerance for the source (eigenvalue) iteration.

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.


int main() 
{		
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
       
  // Create space.
  // Transform input data to the format used by the "Space" constructor.
  SpaceData *md = new SpaceData(verbose);		
  Space* space = new Space(md->N_macroel, md->interfaces, md->poly_orders, md->material_markers, md->subdivisions, N_GRP, N_SLN);  
  delete md;
  
  // Enumerate basis functions, info for user.
  int ndof = Space::get_num_dofs(space);
  info("ndof: %d", ndof);

  // Plot the space.
  space->plot("space.gp");
  
  // Initial approximation of the dominant eigenvalue.
  double K_EFF = 1.0;
  // Initial approximation of the dominant eigenvector.
  double init_val = 1.0;

  for (int g = 0; g < N_GRP; g++)  
  {
    set_vertex_dofs_constant(space, init_val, g);
    space->set_bc_right_dirichlet(g, flux_right_surf[g]);
  }
  
  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, jacobian_fuel_0_0, NULL, fuel);
  wf.add_matrix_form(0, 0, jacobian_water_0_0, NULL, water);

  wf.add_matrix_form(0, 1, jacobian_fuel_0_1, NULL, fuel);
  wf.add_matrix_form(0, 1, jacobian_water_0_1, NULL, water);  

  wf.add_matrix_form(1, 0, jacobian_fuel_1_0, NULL, fuel);
  wf.add_matrix_form(1, 0, jacobian_water_1_0, NULL, water);

  wf.add_matrix_form(1, 1, jacobian_fuel_1_1, NULL, fuel);
  wf.add_matrix_form(1, 1, jacobian_water_1_1, NULL, water);
    
  wf.add_vector_form(0, residual_fuel_0, NULL, fuel);
  wf.add_vector_form(0, residual_water_0, NULL, water);  
  
  wf.add_vector_form(1, residual_fuel_1, NULL, fuel);
  wf.add_vector_form(1, residual_water_1, NULL, water); 

  wf.add_vector_form_surf(0, residual_surf_left_0, BOUNDARY_LEFT);
  wf.add_vector_form_surf(1, residual_surf_left_1, BOUNDARY_LEFT);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem *dp = new DiscreteProblem(&wf, space, is_linear);
  
  Linearizer l(space);
  char solution_file[32];

  // Source iteration
  int i;
  int current_solution = 0, previous_solution = 1;
  double K_EFF_old;
  bool success = false;
  for (i = 0; i < Max_SI; i++)
  {	
    // Plot the critical (i.e. steady-state) flux in the actual iteration.
    sprintf(solution_file, "solution_%d.gp", i);
    l.plot_solution(solution_file); 		
	  
    // Store the previous solution (used at the right-hand side).
    for (int g = 0; g < N_GRP; g++)
      copy_dofs(current_solution, previous_solution, space, g);

    // Obtain the number of degrees of freedom.
    int ndof = Space::get_num_dofs(space);

    // Fill vector coeff_vec using dof and coeffs arrays in elements.
    double *coeff_vec = new double[Space::get_num_dofs(space)];
    get_coeff_vector(space, coeff_vec);
  
    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
    int it = 1;
    while (1) 
    {
      // Obtain the number of degrees of freedom.
      int ndof = Space::get_num_dofs(space);

      // Assemble the Jacobian matrix and residual vector.
      dp->assemble(coeff_vec, matrix, rhs);

      // Calculate the l2-norm of residual vector.
      double res_l2_norm = get_l2_norm(rhs);

      // Info for user.
      info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(space), res_l2_norm);

      // If l2 norm of the residual vector is within tolerance, then quit.
      // NOTE: at least one full iteration forced
      //       here because sometimes the initial
      //       residual on fine mesh is too small.
      if(res_l2_norm < NEWTON_TOL && it > 1) break;

      // Multiply the residual vector with -1 since the matrix 
      // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
      for(int i=0; i<ndof; i++) rhs->set(i, -rhs->get(i));

      // Solve the linear system.
      if(!(success = solver->solve()))
        error ("Matrix solver failed.\n");

      // Add \deltaY^{n+1} to Y^n.
      for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];

      // If the maximum number of iteration has been reached, then quit.
      if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
      
      // Copy coefficients from vector y to elements.
      set_coeff_vector(coeff_vec, space);

      it++;
    }
    
    // Cleanup.
    delete matrix;
    delete rhs;
    delete solver;
    delete [] coeff_vec;
			
    // Update the eigenvalue.
    K_EFF_old = K_EFF;
    K_EFF = calc_total_reaction_rate(space, nSf, 0., 40.); 
    
    // Convergence test.
    if (fabs(K_EFF - K_EFF_old)/K_EFF < TOL_SI) break;
    
    // Normalize total neutron flux to one fission neutron.
    multiply_dofs_with_constant(space, 1./K_EFF, current_solution);
    
    if (verbose) info("K_EFF_%d = %.8f", i+1, K_EFF);
  }
  
  // Print the converged eigenvalue.
  info("K_EFF = %.8f, err= %.8f%%", K_EFF, 100*(K_EFF-1));

  // Comparison with analytical results (see the reference above).
  double flux[N_GRP], J[N_GRP], R;

  get_solution_at_point(space, 0.0, flux, J);
  R = flux[0]/flux[1];
  info("phi_fast/phi_therm at x=0 : %.4f, err = %.2f%%", R, 100*(R-2.5332)/2.5332);
	
  get_solution_at_point(space, 40.0, flux, J);
  R = flux[0]/flux[1];
  info("phi_fast/phi_therm at x=40 : %.4f, err = %.2f%%", R, 100*(R-1.5162)/1.5162);

  info("Total running time: %g s", cpu_time.accumulated());

  // Test variable.
  info("ndof = %d.", Space::get_num_dofs(space));
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
