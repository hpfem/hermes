#include "hermes2d.h"
#include <stdio.h>
#include <cmath>
using namespace std;

// This test makes sure that example 52-eigenvalue-h2p works correctly.

const int NUMBER_OF_EIGENVALUES = 1;             // Desired number of eigenvalues.
int P_INIT =4;  // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 4; // DON'T CHANGE THIS!!
const int FINAL_REF_NUM = 0; // final global refinement 
const double TARGET_VALUE = -2.0;                  // JDSYM parameter: Eigenvalues in the vicinity of this number will be computed. 
const double E_LITERATURE=-2.2052684289898; 
// value from 
// J. Chem. Phys. 43, 3004 (1965); doi:10.1063/1.1697265.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.
const int Nnuc=2; // number of nuclei, here 2
double C[Nnuc]= { 1.0186573603637741, 1.0186573603637741  }; // coefficients in cusp factor 
double Rnuc[Nnuc][3]={  {0.0,0.0,-1.0}, {0.0,0.0,1.0} }; // coordinates of the two nuclei
// the following functions are needed for the cusp factor  version of h2plus
double ri(int i,double rho,double z){
  return sqrt(rho*rho+(z-Rnuc[i][2])*(z-Rnuc[i][2]));
}// this function gives the distances from the point (rho,z)  to the i-th nucleus

double f(double rho,double z){
  return 1.0+C[0]*exp(-2.0*ri(0,rho,z))+C[1]*exp(-2.0*ri(1,rho,z));
}// this is the factor that satisfies the cusp conditions at both nuclei
double laplacef(double rho,double z){
  return 4.0*C[0]*exp(-2.0*ri(0,rho,z))+4.0*C[1]*exp(-2.0*ri(1,rho,z))\
    -4.0/ri(0,rho,z)*C[0]*exp(-2.0*ri(0,rho,z))-4.0/ri(1,rho,z)*C[1]*exp(-2.0*ri(1,rho,z));
}// this is laplace of f

double pot(double rho,double z){
  return -2.0/ri(0,rho,z)-2.0/ri(1,rho,z)-laplacef(rho,z)/f(rho,z);
}// this is the potential with the singularities at the nuclei removed and the 1/r replace by a kink 

double  wfun (double rho,double z){
  return f(rho,z)*f(rho,z);
}// this is the weight function f**2

// Boundary markers.
const int BDY_BOTTOM = 1, BDY_RIGHT = 2, BDY_TOP = 3, BDY_LEFT = 4;

#include "forms.cpp"
// Write the matrix in Matrix Market format.
void write_matrix_mm(const char* filename, Matrix* mat) 
{
  int ndof = mat->get_size();
  FILE *out = fopen(filename, "w" );
  int nz=0;
  for (int i=0; i < ndof; i++) {
    for (int j=0; j <=i; j++) { 
      double tmp = mat->get(i,j);
      if (fabs(tmp) > 1e-15) nz++;
    }
  } 

  fprintf(out,"%%%%MatrixMarket matrix coordinate real symmetric\n");
  fprintf(out,"%d %d %d\n", ndof, ndof, nz);
  for (int i=0; i < ndof; i++) {
    for (int j=0; j <=i; j++) { 
      double tmp = mat->get(i,j);
      if (fabs(tmp) > 1e-15) fprintf(out, "%d %d %24.15e\n", i+1, j+1, tmp);
    }
  } 
  fclose(out);
}

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
  if (vertex_id_nucleus == -1)
    {
      printf("nucleus is not a vertex! Slow convergence would result\n");
      exit(0);
    }
  printf("vertex_id_nucleus=%d\n",vertex_id_nucleus);
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
  printf("ndof = %d\n", ndof);
  // Initialize the weak formulation for the left hand side, i.e., H.
  WeakForm wfH, wfU;
  wfH.add_matrix_form(bilinear_form_H,bilinear_form_ord);
  wfU.add_matrix_form(bilinear_form_U,bilinear_form_ord);

  // Initialize matrices.
  SparseMatrix* Hmat = create_matrix(matrix_solver);
  SparseMatrix* Umat = create_matrix(matrix_solver);

  // Assemble the matrices.
  bool is_linear = true;
  DiscreteProblem dpH(&wfH, &space, is_linear);
  dpH.assemble(Hmat);
  DiscreteProblem dpU(&wfU, &space, is_linear);
  dpU.assemble(Umat);
  cpu_time.tick();
  printf("Total running time for assembling matrices : %g s\n", cpu_time.accumulated());
  cpu_time.reset();
  write_matrix_mm("mat_left.mtx" ,Hmat);
  write_matrix_mm("mat_right.mtx", Umat);
  cpu_time.tick();
  printf("Total running time for writing matrices to disk : %g s\n", cpu_time.accumulated());
  cpu_time.reset();
  printf("Calling JDSYM...\n");
  char call_cmd[255];
  double TOL=1e-10;
  int MAX_ITER=1000;
  sprintf(call_cmd, "python solveGenEigenFromMtx.py mat_left.mtx mat_right.mtx %g %d %g %d", 
	  TARGET_VALUE, NUMBER_OF_EIGENVALUES, TOL, MAX_ITER);
  system(call_cmd);
  printf("JDSYM finished.\n");
  cpu_time.tick();
  printf("Total running time for solving generalized eigenvalue problem: %g s\n", cpu_time.accumulated());
  double* coeff_vec = new double[ndof];
  Solution sln;
  // Reading solution vectors from file and visualizing.
  double* eigenval = new double[NUMBER_OF_EIGENVALUES];
  FILE *file = fopen("eivecs.dat", "r");
  char line [64];                  // Maximum line size.
  fgets(line, sizeof line, file);  // ndof
  int n = atoi(line);            
  if (n != ndof) error("Mismatched ndof in the eigensolver output file.");  
  fgets(line, sizeof line, file);  // Number of eigenvectors in the file.
  int neig = atoi(line);
  if (neig != NUMBER_OF_EIGENVALUES) error("Mismatched number of eigenvectors in the eigensolver output file.");  
  for (int ieig = 0; ieig < neig; ieig++) {
    // Get next eigenvalue from the file
    fgets(line, sizeof line, file);
    eigenval[ieig] = atof(line);            
    // Get the corresponding eigenvector.
    for (int i = 0; i < ndof; i++) {  
      fgets(line, sizeof line, file);
      coeff_vec[i] = atof(line);
    }

    // Convert coefficient vector into a Solution.
    Solution::vector_to_solution(coeff_vec, &space, &sln);
  }  
  fclose(file);

  delete [] coeff_vec;

  double E=eigenval[0];
  printf("E=%.16f   Delta E=%.16e\n",E ,E-E_LITERATURE);

  ndof = Space::get_num_dofs(&space);
  printf("Coordinate ( 0, 0) value = %lf\n", sln.get_pt_value(0, 0));
  printf("Coordinate ( 1, 1) value = %lf\n", sln.get_pt_value(1.0, 1.0));
  printf("Coordinate ( 3, 3) value = %lf\n", sln.get_pt_value(3.0, 3.0));
  printf("Coordinate ( 5, 5) value = %lf\n", sln.get_pt_value(5.0, 5.0));
  printf("Coordinate ( 7, 7) value = %lf\n", sln.get_pt_value(7.0, 7.0));
  printf("Coordinate (10,10) value = %lf\n", sln.get_pt_value(10.0, 10.0));

  double coor_x_y[6] = {0.0, 1.0, 3.0, 5.0, 7.0, 10.0};
  double value[6] = { 0.246677, 0.140036, 0.003992, 0.000072, 0.000001, 0.000000};
  bool success = true;
  for (int i = 0; i < 6; i++)
  {
    if (abs(value[i] - sln.get_pt_value(coor_x_y[i], coor_x_y[i])) > 1E-6) success = false;
  }
  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

