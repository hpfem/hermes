//#define NDEBUG

#include "hermes1d.h"

#include <map>
#include <cassert>

using std::cout;
using std::endl;

// *****************************************************************************

// This example solves a 1D fixed source problem for the neutron diffusion eqn.
// in a two-group approximation	: 
//		-(D1.u1')' + Sa1.u1 = Q 
//		-(D2.u2')' + Sa2.u2 = S12.u1  
// The core is composed of a single, 80cm wide slab. Reflective boundary 
// condition is prescribed on the left end, zero-flux condition on the right end
// (homogeneous b.c. of Neumann/Dirichlet type, respectively). There is
// a uniform source of 1.5 fast neutrons (group 1) per cm per sec. 
//	Reference:
// 		HP-MESH ADAPTATION FOR 1-D MULTIGROUP NEUTRON DIFFUSION PROBLEMS,
// 		A MSc. Thesis by YAQI WANG, Texas A&M University, 2006,
//		Example 4.A (pp. 168)

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
// General input (external source problem):

bool flag = false;											// Flag for debugging purposes
bool verbose = true;										

int N_SLN = 1;              						// Number of solutions

// Newton's method
double NEWTON_TOL = 1e-5;               // tolerance for the Newton's method
int NEWTON_MAXITER = 150;               // max. number of Newton iterations

/******************************************************************************/

int main() {
 // Create coarse mesh
  MeshData *md = new MeshData();		// transform input data to the format used by the "Mesh" constructor
  Mesh *mesh = new Mesh(md->N_macroel, md->interfaces, md->poly_orders, md->material_markers, md->subdivisions, N_GRP, N_SLN);  
  delete md;
  
  printf("N_dof = %d\n", mesh->assign_dofs());
  mesh->plot("mesh.gp");

  for (int g = 0; g < N_GRP; g++)  {
  	mesh->set_bc_right_dirichlet(g, flux_right_surf[g]);
	}
  
  // Register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  
  dp->add_matrix_form(0, 0, jacobian_fuel_0_0, fuel);
  dp->add_matrix_form(0, 1, jacobian_fuel_0_1, fuel);
  dp->add_matrix_form(1, 0, jacobian_fuel_1_0, fuel);    
  dp->add_matrix_form(1, 1, jacobian_fuel_1_1, fuel);
    
  dp->add_vector_form(0, residual_fuel_0, fuel);  
  dp->add_vector_form(1, residual_fuel_1, fuel); 

  dp->add_vector_form_surf(0, residual_surf_left_0, BOUNDARY_LEFT);
  dp->add_vector_form_surf(1, residual_surf_left_1, BOUNDARY_LEFT);
	  	
  // Newton's loop		
  newton(dp, mesh, NULL, NEWTON_TOL, NEWTON_MAXITER, verbose);
	 
  // Plot the resulting neutron flux
  Linearizer l(mesh);
  l.plot_solution("solution.gp");

	// Calculate flux integral for comparison with the reference value
	double I = calc_integrated_flux(mesh, 1, 60., 80.);
	double Iref = 134.9238787715397;
	printf("I = %.13f, err = %.13f%%\n", I, 100.*(I - Iref)/Iref );
	
  printf("Done.\n");
  return 1;
}
