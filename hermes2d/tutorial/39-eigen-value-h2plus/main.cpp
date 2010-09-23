#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#include "python_api.h"
#include "hermes2d.h"
#include <stdio.h>
#include <cmath>
using namespace std;

// This is a solver for the Schroedinger equation of
// the Hydrogen molecular Ion in cylindrical coordinates.
// The nuclei are at (0,1) and (0,-1) in atomic units.
// To make the potential less singular 
// the wave function is written as a product of a function satisfying
// the cusp condition at the nuclei and a smoother function to be determined 
// via the method of finite elements on the domain [0,8]x[0,8] in 
// cylindrical coordinates rho and z. The boundary conditions are 
// essential for rho=8 and z=8 and natural for rho=0 and z=0

int P_INIT = 6;                                   // Uniform polynomial degree of mesh elements.
int zp = 1;                                       // Zparity +1/-1.
const int INIT_REF_NUM = 3;                       // DON'T CHANGE THIS!!
const int FINAL_REF_NUM = 2;                      // Final global refinement. 
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
const int Nnuc = 2;                                              // Number of nuclei, here 2.
double C[Nnuc] = { 1.0186573603637741, 1.0186573603637741  };    // Coefficients in cusp factor.
double Rnuc[Nnuc][3]={  {0.0,0.0,-1.0}, {0.0,0.0,1.0} };         // Coordinates of the two nuclei.

// Boundary condition types.
// Note: "essential" means that solution value is prescribed.
BCType bc_types(int marker)
{
  if (marker == 2 or marker == 3)
    return BC_ESSENTIAL;
  if (marker == 1 and zp == -1) 
    return BC_ESSENTIAL;
  if (marker == 1 and zp == 1)
    return BC_NATURAL;
  if (marker == 4 )
    return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

// The following functions are needed to the cusp factor  version of h2plus.
// This function gives the distances from the point (rho,z) to the i-th nucleus.
double ri(int i,double rho,double z){
  return pow(rho*rho+(z-Rnuc[i][2])*(z-Rnuc[i][2]),0.5);
}

// The factor that satisfies the cusp conditions at both nuclei.
double f(double rho,double z){
  return 1.0+C[0]*exp(-2.0*ri(0,rho,z))+C[1]*exp(-2.0*ri(1,rho,z));
}

// Laplace of f.
double laplacef(double rho,double z){
  return 4.0*C[0]*exp(-2.0*ri(0,rho,z))+4.0*C[1]*exp(-2.0*ri(1,rho,z))\
    -4.0/ri(0,rho,z)*C[0]*exp(-2.0*ri(0,rho,z))-4.0/ri(1,rho,z)*C[1]*exp(-2.0*ri(1,rho,z));
}

// Potential with the singularities at the nuclei removed and the 1/r replace by a kink. 
double pot(double rho,double z){
  return -2.0/ri(0,rho,z)-2.0/ri(1,rho,z)-laplacef(rho,z)/f(rho,z);
}

// Weight function f**2.
double  wfun (double rho,double z){
  return f(rho,z)*f(rho,z);
}

// Python interpreter object.
Python p; 

// Weak forms.
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

  // The Python script below calculates the functions for control purposes to compare them with the C++ implementations.   
  p.exec("\
from gencuspfun import cuspfun \n\
from  numpy import *\n\
ZZ=[1.0,1.0];RNUC=[array([0,0,-1.0]),array([0,0,1.0])]\n\
f,pot=cuspfun(ZZ,RNUC)\n\
import pylab\n\
x=linspace(-4.0,4.0,4000)\n\
px=[pot(0,0,xx) for xx in x]\n\
fx=[f(0,0,xx) for xx in x];pylab.plot(x,px);pylab.plot(x,fx);\n\
pylab.show()\n\
");

  // Perform initial mesh refinements. 
  for (int i=0;i<INIT_REF_NUM;i++) mesh.refine_all_elements();
  int VERTEX_ID_NUCLEUS=66;
  mesh.refine_towards_vertex(VERTEX_ID_NUCLEUS,6);
  /* This following code was used to identify the VERTEX_ID_NUCLEUS! Needs to be used again whenever INIT_REF_NUM is changed! 
  // There actually should be a search function on the mesh for this!
  int NUM_ELEMS=mesh.get_max_element_id();
  for (int id=0;id<NUM_ELEMS;id++) {
    Element* elem=mesh.get_element(id);
    for(int i=0;i<4;i++){
      Node* node= elem->vn[i];
      //printf("element %d: %d, %d,  %f, %f\n",id,node->id,i,node->x,node->y);
	}
  }
  */
  for (int i=0;i<FINAL_REF_NUM;i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Calculate the number of active elements.
  int nelem = mesh.get_num_active_elements();
  info("nelem = %d", nelem);

  // Initialize the weak formulation for the left hand side i.e. H. 
  WeakForm wfH;;
  wfH.add_matrix_form(bilinear_form_H, bilinear_form_ord, H2D_SYM);

  // Initialize the linear problem.
  LinearProblem lpH(&wfH, &space);

  // Select matrix solver.
  Matrix *Hmat,*Vmat,*Umat;
  Vector* eivec; CommonSolver* solver;
  init_matrix_solver(matrix_solver, ndof, Hmat, eivec, solver);

  // Assemble stiffness matrix
  lpH.assemble(Hmat, eivec);

  // Initialize the weak formulation for the potential matrix.
  WeakForm wfPot;
  wfPot.add_matrix_form(bilinear_form_V, bilinear_form_ord, H2D_SYM);

  // Initialize the linear problem.
  LinearProblem lpPot(&wfPot, &space);
  init_matrix_solver(matrix_solver, ndof, Vmat, eivec, solver);

  // Assemble potential matrix.
  lpPot.assemble(Vmat,eivec);

  // Initialize the weak formulation for the right hand side i.e. U. 
  WeakForm wfU;
  wfU.add_matrix_form(bilinear_form_U, bilinear_form_ord, H2D_SYM);

  // Initialize the linear problem.
  LinearProblem lpU(&wfU, &space);
  init_matrix_solver(matrix_solver, ndof, Umat, eivec, solver);

  // Assemble overlap matrix. 
  lpU.assemble(Umat,eivec);
  cpu_time.tick();
  verbose("Total running time for assembling matrices : %g s", cpu_time.accumulated());
  cpu_time.reset();
  CSRMatrix  mat1(Hmat);
  CSRMatrix  mat2(Umat);
  CSRMatrix  mat3(Vmat);
  cpu_time.tick();
  verbose("Time for converting matrices to CSR  : %g s", cpu_time.accumulated());
  cpu_time.reset();
  p.exec("print 'Python initialized'");
  p.push("hmat",c2py_CSRMatrix(&mat1));
  p.push("umat",c2py_CSRMatrix(&mat2));
  p.push("vmat",c2py_CSRMatrix(&mat3));
  p.exec("h_csr = hmat.to_scipy_csr()");
  p.exec("u_csr = umat.to_scipy_csr()");
  p.exec("v_csr = vmat.to_scipy_csr()");
  p.exec("H = h_csr.tocoo()");
  p.exec("U = u_csr.tocoo()");
  p.exec("V = v_csr.tocoo()");
  p.exec("from solvers.eigen import solve_eig_numpy, solve_eig_pysparse");
  p.exec("eigs = solve_eig_pysparse(H+V,U,-2.1,n_eigs=1)");
  p.exec("E, v = eigs[0]");
  double E = py2c_double(p.pull("E"));
  int n;
  double *evec;
  numpy2c_double_inplace(p.pull("v"), &evec, &n);
  printf("E = %.16f\n", E);
  cpu_time.tick();
  verbose("Total running time for solving generalized eigenvalue problem: %g s", cpu_time.accumulated());

  // Convert double array to vector and coefficient vector into a Solution.
  for (int i=0;i<ndof;i++) eivec->set(i, abs(evec[i])); // Ensure that the plotted eigenvector is positive.
  Solution* sln = new Solution(&space, eivec);
  ScalarView view("Solution", new WinGeom(0, 0, 1024, 768));
  double PsiNuc=sln->get_pt_value(0.0, 1.0, H2D_FN_VAL_0);  
  p.push("PsiNuc",c2numpy_double(&PsiNuc,1));
  p.exec("print PsiNuc[0]"); // Value of the wave function at the nuclei.

  // Visualize the solution.
  //view.set_3d_mode();
  view.show(sln);

  // Wait for the view to be closed.
  View::wait(H2DV_WAIT_KEYPRESS);

  return 0; 
};

