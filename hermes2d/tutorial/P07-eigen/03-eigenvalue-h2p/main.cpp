#define HERMES_REPORT_ALL
#include <stdio.h>
#include <cmath>
#include "hermes2d.h"

using Teuchos::ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Hermes::EigenSolver;

// This is a solver for the Schroedinger equation for the ground state of 
// the Hydrogen molecular Ion in cylindrical coordinates the nuclei are at 
// (0,1) and (0,-1) in atomic units. To make the potential less singular 
// the wave function is written as a product of a function satisfying
// the cusp condition at the nuclei and a smoother function to be determined 
// via the method of finite elements on the domain [0,16]x[0,16] in cylindrical 
// coordinates rho and z. The boundary conditions are essential for rho=rhomax 
// and z=zmax and natural for rho=0 and z=0.

const int NUMBER_OF_EIGENVALUES = 1;               // Desired number of eigenvalues.
int P_INIT = 4;                                    // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 4;                        // DO NOT CHANGE THIS!!
const int FINAL_REF_NUM = 0;                       // Final global refinement. 
const double TARGET_VALUE = -2.0;                  // JDSYM parameter: Eigenvalues in the vicinity of this 
                                                   // number will be computed. 
const double E_LITERATURE=-2.2052684289898; 

// Value from J. Chem. Phys. 43, 3004 (1965); doi:10.1063/1.1697265.
const int Nnuc=2;                                            // Number of nuclei, here 2.
double C[Nnuc]= { 1.0186573603637741, 1.0186573603637741 };  // Coefficients in cusp factor.
double Rnuc[Nnuc][3]={  {0.0,0.0,-1.0}, {0.0,0.0,1.0} };     // Coordinates of the two nuclei.

// The following functions are needed for the cusp factor version of h2plus.
// This function gives the distances from the point (rho,z) to the i-th nucleus.
double ri(int i,double rho,double z){
  return sqrt(rho*rho+(z-Rnuc[i][2])*(z-Rnuc[i][2]));
}

// This is the factor that satisfies the cusp conditions at both nuclei.
double f(double rho,double z){
  return 1.0+C[0]*exp(-2.0*ri(0,rho,z))+C[1]*exp(-2.0*ri(1,rho,z));
}

// This is Laplacian of f.
double laplacef(double rho,double z){
  return 4.0*C[0]*exp(-2.0*ri(0,rho,z))+4.0*C[1]*exp(-2.0*ri(1,rho,z))\
    -4.0/ri(0,rho,z)*C[0]*exp(-2.0*ri(0,rho,z))-4.0/ri(1,rho,z)*C[1]*exp(-2.0*ri(1,rho,z));
}

// This is the potential with the singularities at the nuclei removed 
// and the 1/r replace by a kink 
double pot(double rho,double z){
  return -2.0/ri(0,rho,z)-2.0/ri(1,rho,z)-laplacef(rho,z)/f(rho,z);
}

double  wfun (double rho,double z){
  return f(rho,z)*f(rho,z);
}// this is the weight function f**2

// Boundary markers.
const int BDY_BOTTOM = 1, BDY_RIGHT = 2, BDY_TOP = 3, BDY_LEFT = 4;

#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
  for (int i=0;i<INIT_REF_NUM;i++) mesh.refine_all_elements();

  //The following code finds the 
  // vertex id at the nuclear coordinate (0.0,1.0).
  int NUM_ELEMS=mesh.get_max_element_id();
  int vertex_id_nucleus=-1;
  for (int id=0;id<NUM_ELEMS;id++) {
    Element* elem=mesh.get_element(id);
    for(int i=0;i<4;i++){
      Node* node= elem->vn[i];
      if ( node->x == 0.0 && node->y == 1.0) 
	      vertex_id_nucleus=node->id;
    }
  }
  if (vertex_id_nucleus == -1) error("Nucleus is not a vertex! Slow convergence would result.");
  info("Vertex_id_nucleus = %d",vertex_id_nucleus);
  mesh.refine_towards_vertex(vertex_id_nucleus,10);
  for (int i=0;i<FINAL_REF_NUM;i++) mesh.refine_all_elements();

  // Enter boundary markers. 
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_RIGHT, BDY_TOP));

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(Hermes::vector<int>(BDY_RIGHT, BDY_TOP));

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation for the left hand side, i.e., H.
  WeakForm wfH, wfU;
  wfH.add_matrix_form(bilinear_form_H,bilinear_form_ord);
  wfU.add_matrix_form(bilinear_form_U,bilinear_form_ord);

  // Initialize matrices.
  RCP<SparseMatrix> Hmat = rcp(new CSCMatrix());
  RCP<SparseMatrix> Umat = rcp(new CSCMatrix());

  // Assemble the matrices.
  info("Assembling matrices.");
  bool is_linear = true;
  DiscreteProblem dpH(&wfH, &space, is_linear);
  dpH.assemble(Hmat.get());
  DiscreteProblem dpU(&wfU, &space, is_linear);
  dpU.assemble(Umat.get());
  cpu_time.tick();
  info("Total running time for assembling matrices : %g s.", cpu_time.accumulated());

  // Initialize eigensolver.
  cpu_time.reset();
  EigenSolver es(Hmat, Umat);
  cpu_time.tick();
  info("Total running time for initializing EigenSolver : %g s.", cpu_time.accumulated());

  // Solve the generalized eigenproblem.
  cpu_time.reset();
  es.solve(NUMBER_OF_EIGENVALUES, TARGET_VALUE);
  es.print_eigenvalues();
  info("Total running time for solving generalized eigenvalue problem (new approach): %g s.", cpu_time.accumulated());

  double* coeff_vec;
  Solution sln;
  ScalarView view("Solution", new WinGeom(0, 0, 1024, 768));

  // Read solution vectors from file and visualize.
  double* eigenval = new double[NUMBER_OF_EIGENVALUES];
  int neig = es.get_n_eigs();
  int n;
  if (neig != NUMBER_OF_EIGENVALUES) error("Mismatched number of eigenvectors in the eigensolver versus requested number.");
  for (int ieig = 0; ieig < neig; ieig++) {
    eigenval[ieig] = es.get_eigenvalue(ieig);

    // Convert coefficient vector into a Solution.
    es.get_eigenvector(ieig, &coeff_vec, &n);
    Solution::vector_to_solution(coeff_vec, &space, &sln);

    // Visualize the solution.
    char title[100];
    sprintf(title, "Solution %d, val = %g", ieig, eigenval[ieig]);
    view.set_title(title);
    view.show(&sln);

    // Wait for keypress.
    View::wait(HERMES_WAIT_KEYPRESS);
  }

  double E = eigenval[0];
  info("E = %.16f   Delta E = %.16e", E, E - E_LITERATURE);
  return 0;
};
