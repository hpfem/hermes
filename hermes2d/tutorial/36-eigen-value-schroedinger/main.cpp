#define H2D_REPORT_INFO
#include "hermes2d.h"
#include <stdio.h>
// This is a solver for the 2D Schroedinger equation

int P_INIT = 6;                                   // Uniform polynomial degree of mesh elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Boundary condition types.
// Note: "essential" means that solution value is prescribed.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

double pot(double x, double y) {
  return x*x+y*y;
}

// Weak forms.

#include "forms.cpp"
int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements (optional).
  for (int i=0;i<2;i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation for the left hand side i.e. H 
  WeakForm wfH;
  wfH.add_matrix_form(callback(bilinear_form_H));
  //wfH.add_vector_form(callback(linear_form));
  // Initialize the linear problem.
  LinearProblem lpH(&wfH, &space);

  // Select matrix solver.
  Matrix *Hmat,*Vmat,*Umat; Vector* eivec; CommonSolver* solver;
  init_matrix_solver(matrix_solver, ndof, Hmat, eivec, solver);

  // Assemble stiffness matrix
  lpH.assemble(Hmat, eivec);


  // Initialize the weak formulation for the potential matrix 
  WeakForm wfPot;
  wfPot.add_matrix_form(bilinear_form_V, bilinear_form_V_ord, H2D_SYM);
  // Initialize the linear problem.
  LinearProblem lpPot(&wfPot, &space);
  init_matrix_solver(matrix_solver, ndof, Vmat, eivec, solver);
  // Assemble potential matrix
  lpPot.assemble(Vmat,eivec);

  


  // Initialize the weak formulation for the right hand side i.e. U 
  WeakForm wfU;
  wfU.add_matrix_form(callback(bilinear_form_U));
  // Initialize the linear problem.
  LinearProblem lpU(&wfU, &space);
  init_matrix_solver(matrix_solver, ndof, Umat, eivec, solver);
  // Assemble overlap matrix 
  lpU.assemble(Umat,eivec);
  // print H+V as MatrixMarket
  FILE *out = fopen( "hmat.mtx", "w" );
  int nz=0;
  for (int i=0; i < ndof; i++) 
    for (int j=0; j <=i; j++)
      { 
	double tmp=Hmat->get(i,j)+Vmat->get(i,j);
	if (tmp != 0) nz++;
      } 
  fprintf(out,"%%%%MatrixMarket matrix coordinate real symmetric\n");
  fprintf(out,"%d %d %d\n",ndof,ndof,nz);
  for (int i=0; i < ndof; i++) 
    for (int j=0; j <=i; j++)
      { 
	double tmp=Hmat->get(i,j)+Vmat->get(i,j);
	if (tmp != 0) fprintf(out, "%d %d %24.15e\n",i+1,j+1,tmp);
      } 
  fclose(out);
  // print U as MatrixMarket
  out = fopen( "umat.mtx", "w" );
  nz=0;
  for (int i=0; i < ndof; i++) 
    for (int j=0; j <=i; j++)
      { 
	double tmp=Umat->get(i,j);
	if (tmp != 0) nz++;
      } 
  fprintf(out,"%%%%MatrixMarket matrix coordinate real symmetric\n");
  fprintf(out,"%d %d %d\n",ndof,ndof,nz);
  for (int i=0; i < ndof; i++) 
    for (int j=0; j <=i; j++)
      { 
	double tmp=Umat->get(i,j);
	if (tmp != 0) fprintf( out,"%d %d %24.15e\n",i+1,j+1,tmp);
      } 
  fclose(out);
  system("python solveGenEigenFromMtx.py hmat.mtx umat.mtx 1.0 1");
  FILE *file=fopen("eivecs.dat","r");
  char line [64]; /* or other suitable maximum line size */
  fgets ( line, sizeof line, file );
  int N=atoi(line);
  fgets ( line, sizeof line, file );
  int neig=atoi(line);
  for (int ieig=0;ieig<neig;ieig++)
    {
      for (int i=0;i<N;i++){  
	fgets ( line, sizeof line, file );
	eivec->set(i,atof(line));
	}
      // Convert coefficient vector into a Solution.
      Solution* sln = new Solution(&space, eivec);
      printf("value at x=0,y=0 is %24.15e\n",sln->get_pt_value(0.0, 0.0, H2D_FN_VAL_0));
      ScalarView view("Solution", new WinGeom(0, 0, 512, 384));
      // Visualize the solution.
      view.set_3d_mode();
      view.show(sln);
      // Wait for the view to be closed.
      View::wait(H2DV_WAIT_KEYPRESS);
    }  
  fclose(file);
  return 0; 
};

