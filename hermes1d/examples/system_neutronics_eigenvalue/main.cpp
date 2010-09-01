//#define NDEBUG

#include "hermes1d.h"

#include <map>
#include <cassert>

using std::cout;
using std::endl;

// *****************************************************************************

// This example solves a 1D eigenvalue problem for the neutron diffusion equation
// in a two-group approximation	: 
//		-(D1.u1')' + Sa1.u1 = 1/k.(nSf1.u1 + nSf2.u2) 
//		-(D2.u2')' + Sa2.u2 = S12.u1  
// The 80cm core is enclosed by a 30cm reflector on both sides. Only the right 
// half is calculated, with reflective boundary condition on the centerline of
// the core and zero-flux condition on the outermost end of the reflector 
// (homogeneous b.c. of Neumann/Dirichlet type, respectively). Example from
// 	 Yoshikawa and Wakabayashi, Journal of NUCLEAR SCIENCE and TECHNOLOGY (7), 
//	 p. 355-365 (July 1970)
	
/******************************************************************************/
// Problem specification (core geometry, material properties, initial FE mesh)
#include "neutronics_problem_def.cpp"

// Common functions for neutronics problems (requires variable declarations from
// "neutronics_problem_def.cpp")
#include "neutronics_common.cpp"

// Weak forms for the problem (requires variable declarations from
// "neutronics_problem_def.cpp")
#include "forms.cpp"
/******************************************************************************/
// General input (eigenvalue problem):

bool flag = false;					// flag for debugging purposes
bool verbose = true;

// Power method initialization
int Max_SI = 1000;          // Max. number of eigenvalue iterations
int N_SLN = 2;              // Number of solutions

// Newton's method
double NEWTON_TOL = 1e-5;               // tolerance for the Newton's method
int NEWTON_MAXITER = 150;               // max. number of Newton iterations
double TOL_SI = 1e-6;                   // tol. for the source (eigenvalue) iteration

/******************************************************************************/

int main() {		
  // Create coarse mesh
  MeshData *md = new MeshData(verbose);		// transform input data to the format used by the "Mesh" constructor
  Mesh *mesh = new Mesh(md->N_macroel, md->interfaces, md->poly_orders, md->material_markers, md->subdivisions, N_GRP, N_SLN);  
  delete md;
  
	// Enumerate basis functions
  printf("N_dof = %d\n", mesh->assign_dofs());
  mesh->plot("mesh.gp");
  
  double K_EFF = 1.0;         // Initial approximation of the dominant eigenvalue
	double init_val = 1.0;			// Initial approximation of the dominant eigenvector

  for (int g = 0; g < N_GRP; g++)  {
  	set_vertex_dofs_constant(mesh, init_val, g);
  	mesh->set_bc_right_dirichlet(g, flux_right_surf[g]);
	}
  
  // Register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  
  dp->add_matrix_form(0, 0, jacobian_fuel_0_0, fuel);
  dp->add_matrix_form(0, 0, jacobian_water_0_0, water);

 	dp->add_matrix_form(0, 1, jacobian_fuel_0_1, fuel);
  dp->add_matrix_form(0, 1, jacobian_water_0_1, water);  

 	dp->add_matrix_form(1, 0, jacobian_fuel_1_0, fuel);
  dp->add_matrix_form(1, 0, jacobian_water_1_0, water);

  dp->add_matrix_form(1, 1, jacobian_fuel_1_1, fuel);
  dp->add_matrix_form(1, 1, jacobian_water_1_1, water);
    
  dp->add_vector_form(0, residual_fuel_0, fuel);
  dp->add_vector_form(0, residual_water_0, water);  
  
  dp->add_vector_form(1, residual_fuel_1, fuel);
  dp->add_vector_form(1, residual_water_1, water); 

  dp->add_vector_form_surf(0, residual_surf_left_0, BOUNDARY_LEFT);
  dp->add_vector_form_surf(1, residual_surf_left_1, BOUNDARY_LEFT);

	Linearizer l(mesh);
	char solution_file[32];

  // Source iteration

	int i;
  int current_solution = 0, previous_solution = 1;
 	double K_EFF_old;
  for (i = 0; i < Max_SI; i++)
  {	
  	// Plot the critical (i.e. steady-state) flux in the actual iteration
  	sprintf(solution_file, "solution_%d.gp", i);
	  l.plot_solution(solution_file); 		
	  
    // Store the previous solution (used at the right-hand side)
    for (int g = 0; g < N_GRP; g++)
	    copy_dofs(current_solution, previous_solution, mesh, g);

    // Newton's loop		
    newton(dp, mesh, NULL, NEWTON_TOL, NEWTON_MAXITER, verbose);
			
    // Update the eigenvalue
    K_EFF_old = K_EFF;
    K_EFF = calc_total_reaction_rate(mesh, nSf, 0., 40.); 
    
    // Convergence test
    if (fabs(K_EFF - K_EFF_old)/K_EFF < TOL_SI) break;
    
    // Normalize total neutron flux to one fission neutron
    multiply_dofs_with_constant(mesh, 1./K_EFF, current_solution);
    
    if (verbose) printf("K_EFF_%d = %.8f\n", i+1, K_EFF);
  }
  
  // Print the converged eigenvalue
  printf("K_EFF = %.8f, err= %.8f%%\n", K_EFF, 100*(K_EFF-1));

  // Plot the converged critical  neutron flux
  sprintf(solution_file, "solution.gp");
  l.plot_solution(solution_file);

	// Comparison with analytical results (see the reference above)
	double flux[N_GRP], J[N_GRP], R;

	get_solution_at_point(mesh, 0.0, flux, J);
	R = flux[0]/flux[1];
	printf("phi_fast/phi_therm at x=0 : %.4f, err = %.2f%%\n", R, 100*(R-2.5332)/2.5332);
	
	get_solution_at_point(mesh, 40.0, flux, J);
	R = flux[0]/flux[1];
	printf("phi_fast/phi_therm at x=40 : %.4f, err = %.2f%%\n", R, 100*(R-1.5162)/1.5162);
	
  printf("Done.\n");
  return 1;
}
