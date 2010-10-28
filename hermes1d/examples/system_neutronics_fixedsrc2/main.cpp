#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"


// This example solves a 1D fixed source problem for the neutron diffusion eqn.
// in a two-group approximation	: 
//		-(D1.u1')' + Sa1.u1 = nSf1.u1 + nSf2.u2 + Q 
//		-(D2.u2')' + Sa2.u2 = S12.u1  
// The core is composed of a seven 100cm wide assemblies. Zero-flux conditions
// on both the right and left boundaries of the core (homogeneous Dirichlet).
// There is an assembly-wise constant source of fast (group 1) neutrons. 
//	Reference:
// 		HP-Space ADAPTATION FOR 1-D MULTIGROUP NEUTRON DIFFUSION PROBLEMS,
// 		A MSc. Thesis by YAQI WANG, Texas A&M University, 2006,
//		Example 3 (pp. 154)


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

// Newton's method
double NEWTON_TOL = 1e-5;               // Tolerance for the Newton's method.
int NEWTON_MAX_ITER = 150;              // Max. number of Newton iterations.

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.


int main() {
  // Create coarse mesh.
  // Transform input data to the format used by the "Space" constructor.
  SpaceData *md = new SpaceData();		
  Space *space = new Space(md->N_macroel, md->interfaces, md->poly_orders, md->material_markers, md->subdivisions, N_GRP, N_SLN);  
  delete md;
  
  info("N_dof = %d", space->assign_dofs());
  space->plot("space.gp");

  for (int g = 0; g < N_GRP; g++)  {
  	space->set_bc_left_dirichlet(g, flux_left_surf[g]);
  	space->set_bc_right_dirichlet(g, flux_right_surf[g]);
	}
  
  // Initialize the FE problem.
  DiscreteProblem *dp = new DiscreteProblem();
  
  dp->add_matrix_form(0, 0, jacobian_mat1_0_0, mat1);
  dp->add_matrix_form(0, 0, jacobian_mat2_0_0, mat2);
  dp->add_matrix_form(0, 0, jacobian_mat3_0_0, mat3);
  
  dp->add_matrix_form(0, 1, jacobian_mat1_0_1, mat1);
  dp->add_matrix_form(0, 1, jacobian_mat2_0_1, mat2);
  dp->add_matrix_form(0, 1, jacobian_mat3_0_1, mat3);
  
  dp->add_matrix_form(1, 0, jacobian_mat1_1_0, mat1);    
  dp->add_matrix_form(1, 0, jacobian_mat2_1_0, mat2);
  dp->add_matrix_form(1, 0, jacobian_mat3_1_0, mat3);
    
  dp->add_matrix_form(1, 1, jacobian_mat1_1_1, mat1);
	dp->add_matrix_form(1, 1, jacobian_mat2_1_1, mat2);
	dp->add_matrix_form(1, 1, jacobian_mat3_1_1, mat3);
  
  dp->add_vector_form(0, residual_mat1_0, mat1);  
	dp->add_vector_form(0, residual_mat2_0, mat2);  
	dp->add_vector_form(0, residual_mat3_0, mat3);
	    
  dp->add_vector_form(1, residual_mat1_1, mat1);
  dp->add_vector_form(1, residual_mat2_1, mat2); 
  dp->add_vector_form(1, residual_mat3_1, mat3);  
	  	
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
    if(res_norm_squared < NEWTON_TOL*NEWTON_TOL && it > 1) break;

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
	 
  // Plot the resulting neutron flux.
  Linearizer l(space);
  l.plot_solution("solution.gp");
	
  info("Done.");
  return 1;
}
