#define HERMES_REPORT_INFO
#define HERMES_REPORT_ERROR
#include "hermes3d.h"
#include <stdio.h>

// test of refinement towards the origin for hermes3d
// test problem is harmonic oscillator with
// V(x,y,z)=x*x+y*y+z*z with lowest eigenvalue 3 and
// eigenfunction exp(-1/2*(x*x+y*y+z*z))/pi**(3/4)
// INIT_REF_NUM = number of mesh refinements done for all elements
// REF_ORIGIN = number of mesh refinements towards the origin m
using Teuchos::ptr;
using Teuchos::RCP;
using Teuchos::rcp;
using Hermes::EigenSolver;

int NUMBER_OF_EIGENVALUES = 1;
//const int INIT_REF_NUM = 0;                       // Number of initial mesh refinements.
//const int REF_ORIGIN = 5;
double TOL = 1e-10;                               // Pysparse parameter: Error tolerance.
int MAX_ITER = 1000;                              // PySparse parameter: Maximum number of iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Boundary condition types.
// Note: "essential" means that solution value is prescribed.
BCType bc_types(int marker)
{
  if (marker > 0) return H3D_BC_ESSENTIAL;
    else
      return H3D_BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}

// potential and weight function for respective molecule
double TARGET_VALUE = 3.0;

// Weak forms.
#include "../definitions.cpp"

int P_INIT = P;
int P_INIT_X = P_INIT;                                   // Uniform polynomial degree of mesh elements.
int P_INIT_Y = P_INIT;                                   // Uniform polynomial degree of mesh elements.
int P_INIT_Z = P_INIT;                                   // Uniform polynomial degree of mesh elements.
scalar  hopot(double x,double y, double z, scalar &dx, scalar &dy, scalar &dz){
    return x*x+y*y+z*z;
}



int main(int argc, char* argv[])
{
  int INIT_REF_NUM, REF_ORIGIN;
  if (argc <2)
    error("Not enough parameters, provide a test number!");

  int test_type=atoi(argv[1]);
  if ( test_type == 1){
    INIT_REF_NUM = 0;
    REF_ORIGIN = 5;
  }
  else if ( test_type == 2){
    INIT_REF_NUM = 2;
    REF_ORIGIN = 0;
  }
  else
    error("Invalid test number"); 
  TimePeriod cpu_time;

  // Load the mesh.
  info("Loading mesh...");
  Mesh mesh;
  H3DReader mloader;
  cpu_time.reset();
  mloader.load("../box.mesh3d", &mesh);
  cpu_time.tick();
  info("time taken for loading mesh : %g s", cpu_time.accumulated());
  cpu_time.reset();

  // Perform initial mesh refinements
  cpu_time.tick();
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
  int NUM_ELEMS=mesh.get_num_active_elements();
  info("NUM_ELEMS=%d",NUM_ELEMS);

  // iterate over all elements to find those including the origin!
  // There are at most 8 active elements that have the origin 
  // as a vertex  at every step of the refinement
  int list_of_ids[8];
  int io,ir=0;
  while(ir<REF_ORIGIN) {
    io=0;
    for(std::map<unsigned int, Element*>::const_iterator it=mesh.elements.begin(); it != mesh.elements.end(); it++) {
      Element *e=it->second;
      if (e->active) {
	info("element id= %d",it->first);
	for(int i=0;i<8;i++){
	  unsigned int vtx = e->get_vertex(i);    
	  double xx=mesh.vertices[vtx]->x; 
	  double yy=mesh.vertices[vtx]->y; 
	  double zz=mesh.vertices[vtx]->z;
	  if ( xx == 0.0 and yy == 0.0 and zz== 0.0 ){
	    list_of_ids[io]=it->first;
	    io++;
	  }
	}
      }
    }
    for (int i=0;i<8;i++) {
      info("%d %d\n",i,list_of_ids[i]);
      mesh.refine_element(list_of_ids[i], H3D_H3D_H3D_REFT_HEX_XYZ);
    }
    ir++;
    int NUM_ELEMS=mesh.get_max_element_id();
    info("NUM_ELEMS=%d",NUM_ELEMS);
  }
  cpu_time.reset();
  FILE *out = fopen("mesh.gmsh", "w" );
  GmshOutputEngine *gout = new GmshOutputEngine(out);
  gout->out( &mesh);
  fclose(out);

  // setting up space for eigen value calculation with zero boundary conditions
  H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
  bool is_linear = true;
  int ndof = Space::get_num_dofs(&space);
  info("ndof=%d",ndof);
  ExactSolution pot_exact(&mesh,hopot);

  // Initialize the weak formulation for left hand side and right hand side , i.e., H and U.
  info("Initializing weak forms...");
  WeakForm wf_left, wf_right;
  wf_left.add_matrix_form(bilinear_form_left, bilinear_form_ord, HERMES_SYM, HERMES_ANY_INT, Hermes::vector<MeshFunction*>(&pot_exact));
  wf_right.add_matrix_form(bilinear_form_right,bilinear_form_ord, HERMES_SYM, HERMES_ANY_INT);   

  // Initialize matrices and matrix solver.
  RCP<CSCMatrix>  matrix_right = rcp(new CSCMatrix());
  info("Assembling RHS matrix....");
  cpu_time.reset();
  DiscreteProblem dp_right(&wf_right, &space, is_linear);
  dp_right.assemble(matrix_right.get());
  cpu_time.tick();
  info("time taken to assemble RHS matrix: %g s", cpu_time.accumulated());
  Solution sln(&mesh);
  fflush(stdout);
  DiscreteProblem dp_left(&wf_left, &space, is_linear);
  RCP<CSCMatrix> matrix_left = rcp(new CSCMatrix());
  Solver* solver = create_linear_solver(matrix_solver, matrix_left.get());
  cpu_time.reset();
  dp_left.assemble(matrix_left.get());
  cpu_time.tick();
  info("time taken for assembling LHS matrix : %g", cpu_time.accumulated());

  // Initialize eigensolver
  cpu_time.reset();
  EigenSolver es(matrix_left, matrix_right);
  cpu_time.tick();
  info("Total running time for preparing generalized eigenvalue problem: %g s\n", cpu_time.accumulated());
  info("Using eigensolver...");
  cpu_time.reset();
  es.solve(NUMBER_OF_EIGENVALUES, TARGET_VALUE, TOL, MAX_ITER);
  cpu_time.tick();
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
    eival[ieig]=es.get_eigenvalue(ieig);
    info("eival[%d]=%24.15E",ieig,eival[ieig]);
    out_fn_vtk(&sln, "phi", ieig );
    }

  return 0; 
};
