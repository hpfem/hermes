#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <getopt.h>
#include <hermes3d.h>

// Solving a simple heat equation to demonstrate how to use CUBIT with Hermes3D.
//
// Use mesh file 'cylinder2.e. Material IDs corresponds to elements markers, 
// sideset IDs correspond to face (BC) markers.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 1;		          // Number of initial uniform mesh refinements.
const int P_INIT_X = 1, 
          P_INIT_Y = 1, 
          P_INIT_Z = 1;                           // Initial polynomial degree of all mesh elements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.
bool do_output = true;				  // Generate output files (if true).

// Usage info.
void usage() {
  printf("Usage:\n");
  printf("\n");
  printf("  heat-cubit [options] <mesh-file>\n");
  printf("\n");
  printf("Options:\n");
  printf("  --no-output         - do not generate output files\n");
  printf("\n");
}


// Boundary condition types. 
BCType bc_types(int marker)
{
  if (marker == 1) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values. 
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 10;
}

// Weak forms.
#include "forms.cpp"

// Solution output.
void out_fn(MeshFunction *x, const char *name)
{
  char of_name[1024];
  FILE *ofile;
  sprintf(of_name, "%s.vtk", name);
  ofile = fopen(of_name, "w");
  if (ofile != NULL) {
    VtkOutputEngine output(ofile);
    output.out(x, name);
    fclose(ofile);
  }
  else warning("Can not open '%s' for writing.", of_name);
}

// Boundary conditions output.
void out_bc(Mesh *mesh, const char *name)
{
  char of_name[1024];
  FILE *ofile;
  sprintf(of_name, "%s.vtk", name);
  ofile = fopen(of_name, "w");
  if (ofile != NULL) {
    VtkOutputEngine output(ofile);
    output.out_bc(mesh, name);
    fclose(ofile);
  }
  else warning("Can not open '%s' for writing.", of_name);
}

int main(int argc, char **argv)
{
  // Load the inital mesh.
  Mesh mesh;
  ExodusIIReader mesh_loader;
  if (!mesh_loader.load("cylinder2.e", &mesh))
    error("Loading mesh file '%s'\n", "cylinder2.e");

  // Perform initial mesh refinements.
  printf("Performing %d initial mesh refinements.\n", INIT_REF_NUM);
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
  Word_t (nelem) = mesh.get_num_elements();
  printf("New number of elements is %d.\n", (int) nelem);

  // Create H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
  printf("Number of DOF: %d.\n", space.get_num_dofs());

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form1<double, scalar>, bilinear_form1<Ord, Ord>, HERMES_SYM, 1);
  wf.add_matrix_form(bilinear_form2<double, scalar>, bilinear_form2<Ord, Ord>, HERMES_SYM, 2);
  wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY);

  // Initialize the coarse mesh problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Assemble stiffness matrix and rhs.
  printf("  - Assembling... "); fflush(stdout);
  Timer tmr_assemble;
  tmr_assemble.start();
  dp.assemble(matrix, rhs);
  tmr_assemble.stop();
  printf("done in %s (%lf secs).\n", tmr_assemble.get_human_time(), tmr_assemble.get_seconds());
	
  // Solve the matrix problem.
  printf("  - Solving... "); fflush(stdout);
  Timer tmr_solve;
  tmr_solve.start();
  bool solved = solver->solve();
  tmr_solve.stop();
  if (solved) printf("done in %s (%lf secs).\n", tmr_solve.get_human_time(), tmr_solve.get_seconds());
  else error("Failed.\n");

  // Construct a solution.
  Solution sln(&mesh);
  sln.set_coeff_vector(&space, solver->get_solution());

  // Output solution and boundary conditions.
  if (do_output) 
  {
    printf("Solution and BC output.\n");
    out_fn(&sln, "sln");
    out_bc(&mesh, "bc");
  }

  printf("Done.\n");
  return 1;
}
