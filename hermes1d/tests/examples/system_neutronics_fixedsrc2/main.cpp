#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

// This test makes sure that example "system_neutronics_fixedsrc2" works correctly.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

// Problem specification (core geometry, material properties, initial FE space).
#include "neutronics_problem_def.cpp"

// Common functions for neutronics problems (requires variable declarations from
// "neutronics_problem_def.cpp").
#include "neutronics_common.cpp"

// Weak forms for the problem (requires the rhs construction routine from 
// "neutronics_common.cpp").
#include "forms.cpp"

// General input (external source problem).
bool flag = false;											// Flag for debugging purposes.
bool verbose = true;										

int N_SLN = 1;              						// Number of solutions.

// Newton's method.
double NEWTON_TOL = 1e-5;               // Tolerance.
int NEWTON_MAX_ITER = 150;              // Max. number of Newton iterations.

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

int main() 
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Create space.
  // Transform input data to the format used by the "Space" constructor.
  SpaceData *md = new SpaceData();		
  Space* space = new Space(md->N_macroel, md->interfaces, md->poly_orders, md->material_markers, md->subdivisions, N_GRP, N_SLN);  
  delete md;
  
  // Enumerate basis functions, info for user.
  int ndof = Space::get_num_dofs(space);
  info("ndof: %d", ndof);

  // Plot the space.
  space->plot("space.gp");

  for (int g = 0; g < N_GRP; g++)  
  {
    space->set_bc_left_dirichlet(g, flux_left_surf[g]);
    space->set_bc_right_dirichlet(g, flux_right_surf[g]);
  }
  
  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, jacobian_mat1_0_0, NULL, mat1);
  wf.add_matrix_form(0, 0, jacobian_mat2_0_0, NULL, mat2);
  wf.add_matrix_form(0, 0, jacobian_mat3_0_0, NULL, mat3);
  
  wf.add_matrix_form(0, 1, jacobian_mat1_0_1, NULL, mat1);
  wf.add_matrix_form(0, 1, jacobian_mat2_0_1, NULL, mat2);
  wf.add_matrix_form(0, 1, jacobian_mat3_0_1, NULL, mat3);
  
  wf.add_matrix_form(1, 0, jacobian_mat1_1_0, NULL, mat1);    
  wf.add_matrix_form(1, 0, jacobian_mat2_1_0, NULL, mat2);
  wf.add_matrix_form(1, 0, jacobian_mat3_1_0, NULL, mat3);
    
  wf.add_matrix_form(1, 1, jacobian_mat1_1_1, NULL, mat1);
  wf.add_matrix_form(1, 1, jacobian_mat2_1_1, NULL, mat2);
  wf.add_matrix_form(1, 1, jacobian_mat3_1_1, NULL, mat3);
  
  wf.add_vector_form(0, residual_mat1_0, NULL, mat1);  
  wf.add_vector_form(0, residual_mat2_0, NULL, mat2);  
  wf.add_vector_form(0, residual_mat3_0, NULL, mat3);
	    
  wf.add_vector_form(1, residual_mat1_1, NULL, mat1);
  wf.add_vector_form(1, residual_mat2_1, NULL, mat2); 
  wf.add_vector_form(1, residual_mat3_1, NULL, mat3);  

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem *dp = new DiscreteProblem(&wf, space, is_linear);
  
  // Newton's loop.
  // Fill vector coeff_vec using dof and coeffs arrays in elements.
  double *coeff_vec = new double[Space::get_num_dofs(space)];
  get_coeff_vector(space, coeff_vec);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  int it = 1;
  bool success = false;
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
