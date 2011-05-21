#define HERMES_REPORT_INFO
#define HERMES_REPORT_ERROR
#include "config.h"
#include "hermes3d.h"
#include <stdio.h>

// describe the purpose i.e. Hartree Fock for Helium 

using Teuchos::ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Hermes::EigenSolver;

int NUMBER_OF_EIGENVALUES=1;
const int INIT_REF_NUM = 0;                       // Number of initial mesh refinements.
double TOL = 1e-10;                               // Pysparse parameter: Error tolerance.
int MAX_ITER = 1000;                              // PySparse parameter: Maximum number of iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const int MAX_SCF_ITER = 8;                                          

// Boundary condition types.
// Note: "essential" means that solution value is prescribed.
BCType bc_types(int marker)
{
  if (marker > 0) return H3D_BC_ESSENTIAL;
  else return H3D_BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values_eigen(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}
scalar essential_bc_values_poisson(int ess_bdy_marker, double x, double y, double z)
{
  return 2.0 / sqrt(x*x + y*y + z*z);
}

// Potential and weight function for respective molecule.
#include "mol.cpp"
double TARGET_VALUE = E0;

// Weak forms.
#include "definitions.cpp"

int P_INIT = P;
int P_INIT_X = P_INIT;                                   // Uniform polynomial degree of mesh elements (x-direction).
int P_INIT_Y = P_INIT;                                   // Uniform polynomial degree of mesh elements (y-direction).
int P_INIT_Z = P_INIT;                                   // Uniform polynomial degree of mesh elements (z-direction).

double expectation_value (double* vec,  SparseMatrix* A , int ndof)
{
  double* tmpvec = new double[ndof];
  A->multiply_with_vector(vec, tmpvec);
  return vec_dot(vec,tmpvec,ndof);
}

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;

  // Load the mesh.
  info("Loading mesh...");
  Mesh mesh;
  H3DReader mloader;
  cpu_time.reset();
  mloader.load("mol.mesh3d", &mesh);
  cpu_time.tick();
  info("Time taken for loading mesh : %g s", cpu_time.accumulated());
  cpu_time.reset();

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
  int NUM_ELEMS = mesh.get_max_element_id();
  info("NUM_ELEMS = %d", NUM_ELEMS);
  cpu_time.tick();
  info("Time for refining mesh: %g s", cpu_time.accumulated());
  cpu_time.reset();

  // Setting up space for eigen value calculation with zero boundary conditions.
  H1Space space(&mesh, bc_types, essential_bc_values_eigen, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));

  // Setting up space for solution of  Poisson equation with (approximate) boundary conditions 2*N/sqrt(x*x+y*y+z*z).
  H1Space space_poisson(&mesh, bc_types, essential_bc_values_poisson, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
  bool is_linear = true;

  // Setting up Laplace matrix for solving the Poisson equation. 
  WeakForm wf_poisson;
  wf_poisson.add_matrix_form(bilinear_form_laplace, bilinear_form_ord1, HERMES_SYM, HERMES_ANY_INT);
  RCP<SparseMatrix> matrix_Laplace = rcp(new CSCMatrix());
  Solver* solver = create_linear_solver(matrix_solver, matrix_Laplace.get());
  DiscreteProblem dp_poisson(&wf_poisson, &space, is_linear);
  dp_poisson.assemble(matrix_Laplace.get());
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);
  ExactSolution pot_exact(&mesh, pot);
  ExactSolution wfun_exact(&mesh, wfun);
  Solution coul_pot(space.get_mesh());
  coul_pot.set_zero();                 // coul_pot zero at the beginning.

  // Initialize the weak formulation for the left hand side, i.e., H.
  info("Initializing weak form...");
  WeakForm wf_left, wf_right;
  wf_left.add_matrix_form(bilinear_form_left, bilinear_form_ord, HERMES_SYM, HERMES_ANY_INT,
                          Hermes::vector<MeshFunction*>(&pot_exact, &wfun_exact,&coul_pot ));
  wf_right.add_matrix_form(bilinear_form_right, bilinear_form_ord, HERMES_SYM, HERMES_ANY_INT, &wfun_exact);   
  DiscreteProblem dp(&wf_left, &space, is_linear);

  // Initialize matrices and matrix solver.
  RCP<SparseMatrix> matrix_left = rcp(new CSCMatrix());
  RCP<SparseMatrix> matrix_right = rcp(new CSCMatrix());
  cpu_time.reset();
  info("Assembling RHS matrix....");
  DiscreteProblem dp_left(&wf_left, &space, is_linear);
  DiscreteProblem dp_right(&wf_right, &space, is_linear);
  dp_right.assemble(matrix_right.get());
  cpu_time.tick();
  info("time taken to assemble RHS matrix: %g s", cpu_time.accumulated());
  WeakForm wf_coulomb;
  wf_coulomb.add_matrix_form(bilinear_form_coul_pot, bilinear_form_ord1, HERMES_SYM, HERMES_ANY_INT,
                             Hermes::vector<MeshFunction*>(&wfun_exact,&coul_pot ));
  RCP<SparseMatrix> matrix_coulomb = rcp(new CSCMatrix());
  DiscreteProblem dp_coulomb(&wf_coulomb, &space, is_linear);
  Solution sln(space.get_mesh());
  RCP<SparseMatrix> matrix_U = rcp(new CSCMatrix());
  bool DONE=false;
  int iter=0;
  info("SELF CONSISTENT LOOP BEGINS:");
  fflush(stdout);

  // Now follows the self consistent loop:
  while (!DONE){
    WeakForm  wf;
    Vector* rhs=create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix_left.get());
    dp_left.assemble(matrix_left.get());
    cpu_time.tick();
    info("time taken for assembling LHS matrix : %g", cpu_time.accumulated());
    cpu_time.reset();
    dp_coulomb.assemble(matrix_coulomb.get());
    cpu_time.tick();
    info("time for assembling  matrix_coulomb: %g  " , cpu_time.accumulated());

    // Initialize eigensolver.
    cpu_time.reset();
    EigenSolver es(matrix_left, matrix_right);
    cpu_time.tick();
    info("Total running time for preparing generalized eigenvalue problem: %g s\n", cpu_time.accumulated());

    cpu_time.reset();
    info("Using eigensolver...");
    es.solve(NUMBER_OF_EIGENVALUES, TARGET_VALUE, TOL, MAX_ITER);
    info("Total running time for solving generalized eigenvalue problem: %g s", cpu_time.accumulated());
    double* coeff_vec;
    double coulomb_energy;
    int neig = es.get_n_eigs();
    double *eival = new double[neig]; 
    if (neig != NUMBER_OF_EIGENVALUES) error("Mismatched number of eigenvectors in eigensolver");  
    for (int ieig = 0; ieig < neig; ieig++) {
      int n;
      es.get_eigenvector(ieig, &coeff_vec, &n);

      // Convert coefficient vector into a Solution.
      Solution::vector_to_solution(coeff_vec, &space, &sln);
      double norm2=expectation_value(coeff_vec, matrix_right.get(),ndof);
      coulomb_energy=expectation_value(coeff_vec,matrix_coulomb.get(),ndof);
      eival[ieig]=es.get_eigenvalue(ieig);
      info("eigenvector %d : norm=%25.16f \n eigenvalue=%25.16f \n coulomb_energy=%25.16f ",
           ieig, pow(norm2, 0.5), eival[ieig], coulomb_energy);
      wf.add_vector_form(linear_form_poisson, linear_form_poisson_ord, HERMES_ANY_INT,
                         Hermes::vector<MeshFunction*>(&sln, &wfun_exact)); 
      out_fn_vtk(&sln, "phi", ieig);
    }  
    double HFenergy = 2*eival[0] - coulomb_energy;
    info("HF energy for two electrons=%25.16f", HFenergy);
    DiscreteProblem dp_density(&wf, &space, is_linear);
    dp_density.assemble(matrix_U.get(), rhs);
    Solver* solver_poisson = create_linear_solver(matrix_solver, matrix_Laplace.get(),rhs);
    info("Solving the matrix problem.");
    fflush(stdout);
    if(solver_poisson->solve())
	Solution::vector_to_solution(solver_poisson->get_solution(), &space_poisson, &coul_pot);
    else
      error ("Matrix solver failed.\n");
    out_fn_vtk(&coul_pot, "coul_pot", 0);
    iter++; 
    if (iter > MAX_SCF_ITER){ 
      delete [] coeff_vec;
      DONE=true;    
    } 
  }
  return 0; 
};

