//#define NDEBUG

#include "hermes1d.h"

#include <map>
#include <cassert>

using std::cout;
using std::endl;

// *****************************************************************************

// This example solves a 1D fixed source problem for the neutron diffusion eqn.
// in a two-group approximation	: 
//		-(D1.u1')' + Sa1.u1 = nSf1.u1 + nSf2.u2 + Q 
//		-(D2.u2')' + Sa2.u2 = S12.u1  
// The core is composed of a seven 100cm wide assemblies. Zero-flux conditions
// on both the right and left boundaries of the core (homogeneous Dirichlet).
// There is an assembly-wise constant source of fast (group 1) neutrons. 
//	Reference:
// 		HP-MESH ADAPTATION FOR 1-D MULTIGROUP NEUTRON DIFFUSION PROBLEMS,
// 		A MSc. Thesis by YAQI WANG, Texas A&M University, 2006,
//		Example 3 (pp. 154)

/******************************************************************************/
// Problem specification (core geometry, material properties, initial FE mesh)
#include "neutronics_problem_def.cpp"

// Common functions for neutronics problems (requires variable declarations from
// "neutronics_problem_def.cpp")
#include "neutronics_common.cpp"

// Weak forms for the problem (requires the rhs construction routine from 
// "neutronics_common.cpp")
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
  	mesh->set_bc_left_dirichlet(g, flux_left_surf[g]);
  	mesh->set_bc_right_dirichlet(g, flux_right_surf[g]);
	}
  
  // Register weak forms
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
	  	
  // Newton's loop		
  newton(dp, mesh, NULL, NEWTON_TOL, NEWTON_MAXITER, verbose);
	 
  // Plot the resulting neutron flux
  Linearizer l(mesh);
  l.plot_solution("solution.gp");
	
  printf("Done.\n");
  return 1;
}
