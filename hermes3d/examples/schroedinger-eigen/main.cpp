#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#include "config.h"
#include "hermes3d.h"
#include <stdio.h>
#define HERMES__REPORT_INFO
//  This example solves the eigenproblem for the time-independent Schroedinger 
//  equation in a cube with zero boundary conditions. Python and Pysparse must
//  be installed. 
//
//  PDE: -Laplace u + V(x,y,z) u = lambda_k u,
//  where lambda_0, lambda_1, ... are the eigenvalues.
//
//  Domain: Cube (-pi/2, pi/2)^3.
//
//  BC:  Homogeneous Dirichlet.
//
//  The eigenvalues are of the form 
//  i^2+j^2+k^2, where i,j,k are natural numbers 
//  The following parameters can be changed:

int NUMBER_OF_EIGENVALUES = 1;                    // Desired number of eigenvalues.

int P_INIT_X = 4;                                   // Uniform polynomial degree of mesh elements.
int P_INIT_Y = 4;                                   // Uniform polynomial degree of mesh elements.
int P_INIT_Z = 4;                                   // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial mesh refinements.
double TARGET_VALUE = 3.0;
// PySparse parameter: Eigenvalues in the vicinity of this number will be computed. 
double TOL = 1e-10;                               // Pysparse parameter: Error tolerance.
int MAX_ITER = 1000;                              // PySparse parameter: Maximum number of iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK

// Boundary condition types.
// Note: "essential" means that solution value is prescribed.
BCType bc_types(int marker)
{
  if (marker > 0) return BC_ESSENTIAL;
    else
      return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}

// Weak forms.
double V( double x, double y, double z){
  return 0.0; //-1./sqrt(x*x + y*y + z*z);
}
#include "forms.cpp"

// Write the matrix in Matrix Market format.
void write_matrix_mm(const char* filename, Matrix* mat) 
{
  // Get matrix size.
  int ndof = mat->get_size();
  FILE *out = fopen(filename, "w" );
  if (out == NULL) error("failed to open file for writing.");

  // Calculate the number of nonzeros.
  int nz = 0;
  for (int i = 0; i < ndof; i++) {
    for (int j = 0; j <= i; j++) { 
      double tmp = mat->get(i, j);
      if (fabs(tmp) > 1e-15) nz++;
    }
  } 

  // Write the matrix in MatrixMarket format
  fprintf(out,"%%%%MatrixMarket matrix coordinate real symmetric\n");
  fprintf(out,"%d %d %d\n", ndof, ndof, nz);
  for (int i = 0; i < ndof; i++) {
    for (int j = 0; j <= i; j++) { 
      double tmp = mat->get(i, j);
      if (fabs(tmp) > 1e-15) fprintf(out, "%d %d %24.15e\n", i + 1, j + 1, tmp);
    }
  } 

  fclose(out);
}

int main(int argc, char* argv[])
{
  info("Desired number of eigenvalues: %d.", NUMBER_OF_EIGENVALUES);
  // Time measurement.
  TimePeriod cpu_time;
  // Load the mesh.
  info("Loading mesh...");
  Mesh mesh;
  H3DReader mloader;
  cpu_time.reset();
  mloader.load("hexahedron.mesh3d", &mesh);
  cpu_time.tick();
  printf("Total running time for loading mesh : %g s\n", cpu_time.accumulated());
  cpu_time.reset();


  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
  int ndof = Space::get_num_dofs(&space);
  //info("ndof: %d.", ndof);
  printf("ndof: %d.\n", ndof);
  // Initialize the weak formulation for the left hand side, i.e., H.
  info("Initializing weak form...");
  WeakForm wf_left, wf_right;
  wf_left.add_matrix_form(bilinear_form_left, bilinear_form_left_ord, HERMES_SYM, HERMES_ANY );
  wf_right.add_matrix_form(callback(bilinear_form_right), HERMES_SYM, HERMES_ANY );

  // Initialize matrices and matrix solver.
  SparseMatrix* matrix_left = create_matrix(matrix_solver);
  SparseMatrix* matrix_right = create_matrix(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix_left);

  // Assemble the matrices.
  info("Assembling matrices...");
  bool is_linear = true;
  DiscreteProblem dp_left(&wf_left, &space, is_linear);
  dp_left.assemble(matrix_left);
  DiscreteProblem dp_right(&wf_right, &space, is_linear);
  dp_right.assemble(matrix_right);
  cpu_time.tick();
  printf("Total running time for assembling matrices : %g s\n", cpu_time.accumulated());
  cpu_time.reset();
  // Write matrix_left in MatrixMarket format.
  write_matrix_mm("mat_left.mtx", matrix_left);

  // Write matrix_right in MatrixMarket format.
  write_matrix_mm("mat_right.mtx", matrix_right);

  cpu_time.tick();
  printf("Total running time for writing matrices to disk : %g s\n", cpu_time.accumulated());
  cpu_time.reset();

  // Calling Python eigensolver. Solution will be written to "eivecs.dat".
  info("Using eigensolver...");
  char call_cmd[255];
  sprintf(call_cmd, "python solveGenEigenFromMtx.py mat_left.mtx mat_right.mtx %g %d %g %d", 
	  TARGET_VALUE, NUMBER_OF_EIGENVALUES, TOL, MAX_ITER);
  system(call_cmd);
  cpu_time.tick();
  printf("Total running time for solving generalized eigenvalue problem: %g s\n", cpu_time.accumulated());
  // Initializing solution vector, solution and ScalarView.
  info("Initializing solution vector...");
  double* coeff_vec = new double[ndof];
  Solution sln(space.get_mesh());


  // Reading solution vectors from file and visualizing.
  info("Reading solution vectors from file and saving as solutions in paraview format...");
  FILE *file = fopen("eivecs.dat", "r");
  char line [64];                  // Maximum line size.
  fgets(line, sizeof line, file);  // ndof
  int n = atoi(line);            
  if (n != ndof) error("Mismatched ndof in the eigensolver output file.");  
  fgets(line, sizeof line, file);  // Number of eigenvectors in the file.
  int neig = atoi(line);
  if (neig != NUMBER_OF_EIGENVALUES) error("Mismatched number of eigenvectors in the eigensolver output file.");  
  for (int ieig = 0; ieig < neig; ieig++) {
    // Get next eigenvector from the file.
    for (int i = 0; i < ndof; i++) {  
      fgets(line, sizeof line, file);
      coeff_vec[i] = atof(line);
    }

    // Convert coefficient vector into a Solution.
    Solution::vector_to_solution(coeff_vec, &space, &sln);

    out_fn_vtk(&sln, "sln", ieig );
  }  
  fclose(file);

  delete [] coeff_vec;

  return 0; 
};

