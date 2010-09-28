#include "config.h"
#include <getopt.h>
#include <hermes3d.h>

// With large K, this is a singularly perturbed problem that exhibits an extremely
// thin and steep boundary layer. Singularly perturbed problems are considered to
// be very difficult, but you'll see that Hermes can solve them easily even for large
// values of K.
//
// PDE: -Laplace u + K*K*u = K*K.
//
// Domain: cube (0, 1) x (0, 1) x (0, 1), see the file singpert-aniso.mesh3d.
//
// BC:  Homogeneous Dirichlet.
//
// The following parameters can be changed:

const int INIT_REF_NUM = 3;		          // Number of initial uniform mesh refinements.
const int P_INIT_X = 1, 
          P_INIT_Y = 1, 
          P_INIT_Z = 1;                           // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;		          // Error threshold for element refinement of the adapt(...) function 
						  // (default) STRATEGY = 0 ... refine elements elements until sqrt(THRESHOLD) 
						  // times total error is processed. If more elements have similar errors, 
						  // refine all to keep the mesh symmetric.
						  // STRATEGY = 1 ... refine all elements whose error is larger
						  // than THRESHOLD times maximum element error.
const double ERR_STOP = 1.0;			  // Stopping criterion for adaptivity (rel. error tolerance between the
						  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;			  // Adaptivity process stops when the number of degrees of freedom grows
						  // over this limit. This is to prevent h-adaptivity to go on forever.
bool do_output = true;				  // Generate output files (if true).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Problem constants.
const double K_squared = 1e4;		 	  // Equation parameter.

// Boundary condition types. 
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values. 
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}

// Weak forms. 
#include "forms.cpp"

// Mesh output.
void out_orders(Space *space, const char *name, int iter)
{
  char fname[1024];
  sprintf(fname, "iter-%s-%d.vtk", name, iter);
  FILE *f = fopen(fname, "w");
  if (f != NULL) {
    VtkOutputEngine vtk(f);
    vtk.out_orders(space, name);
    fclose(f);
  }
  else
    warning("Could not open file '%s' for writing.", fname);
}

// Solution output.
void out_fn(MeshFunction *fn, const char *name, int iter)
{
  char fname[1024];
  sprintf(fname, "iter-%s-%d.vtk", name, iter);
  FILE *f = fopen(fname, "w");
  if (f != NULL) {
    VtkOutputEngine vtk(f);
    vtk.out(fn, name);
    fclose(f);
  }
  else warning("Could not open file '%s' for writing.", fname);
}

int main(int argc, char **args) {

  // Load the inital mesh.
  Mesh mesh;
  Mesh3DReader mesh_loader;
  if(!mesh_loader.load("singpert-aniso.mesh3d", &mesh))
    error("Loading mesh file '%s'\n", "cylinder2.e");

  // Perform initial mesh refinements.
  printf("Performing %d initial mesh refinements.\n", INIT_REF_NUM);
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
  Word_t (nelem) = mesh.get_num_elements();
  printf("New number of elements is %d.\n", (int) nelem);

  // Graphs of DOF convergence.
  GnuplotGraph graph;
  graph.set_captions("", "Degrees of freedom", "Error [%]");
  graph.set_log_y();
  graph.add_row("Total error", "k", "-", "O");

  // Create H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(biform<double, double>, biform<Ord, Ord>, HERMES_SYM, HERMES_ANY);
  wf.add_vector_form(liform<double, double>, liform<Ord, Ord>, HERMES_ANY);

  // Initialize the coarse mesh problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Adaptivity loop.
  int as = 0; bool done = false;
  do {
    printf("\n---- Adaptivity step %d:\n", as);
    printf("\nSolving on coarse mesh:\n");

    // Procedures for coarse mesh problem.
    // Assign DOF.
    int ndof = space.assign_dofs();
    printf("  - Number of DOF: %d\n", ndof);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_solver(matrix_solver, matrix, rhs);

    // Assemble stiffness matrix and rhs.
    printf("  - Assembling...\n"); fflush(stdout);
    dp.assemble(matrix, rhs);
    
    // Solve the matrix problem.
    printf("  - Solving... "); fflush(stdout);
    bool solved = solver->solve();
    if (solved) printf("done in %lf secs\n", solver->get_time());
    else error("Failed\n");
    
    // Construct a solution.
    Solution sln(&mesh);
    sln.set_coeff_vector(&space, solver->get_solution());

    // Output solution and mesh.
    if (do_output) 
    {
      out_fn(&sln, "sln", as);
      out_orders(&space, "order", as);
    }

    // Solving the fine mesh problem.
    printf("Solving on fine mesh:\n");

    // Construct the refined mesh for reference(refined) solution. 
    Mesh rmesh;
    rmesh.copy(mesh);
    rmesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

    // Setup space for the reference solution.
    Space *rspace = space.dup(&rmesh);
    rspace->copy_orders(space, 1);

    // Initialize the mesh problem for reference solution.
    DiscreteProblem rdp(&wf, rspace, is_linear);

    // Assign DOF.
    int rndof = rspace->assign_dofs();
    printf("  - Number of DOF: %d\n", rndof);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* rmatrix = create_matrix(matrix_solver);
    Vector* rrhs = create_vector(matrix_solver);
    Solver* rsolver = create_solver(matrix_solver, rmatrix, rrhs);

    // Assemble stiffness matric and rhs.
    printf("  - Assembling...\n"); fflush(stdout);
    rdp.assemble(rmatrix, rrhs);

    // Solve the system.
    bool rsolved = rsolver->solve();
    if (rsolved) printf("done in %lf secs\n", rsolver->get_time());
    else printf("failed\n");

    // Construct the reference solution.
    Solution rsln(&rmesh);
    rsln.set_coeff_vector(rspace, rsolver->get_solution());

    // Calculate error estimates for hp-adaptivity.
    printf("Adaptivity\n");
    printf("  - calculating error: "); fflush(stdout);
    H1Adapt hp(&space);
    hp.set_aniso(true);							// anisotropic adaptivity.
    double err_est_rel = hp.calc_error(&sln, &rsln) * 100;
    printf("% lf %%\n", err_est_rel);

    // Save it to the graph.
    graph.add_value(0, ndof, err_est_rel);
    if (do_output) graph.save("conv.gp");


    // If error is too large, adapt the mesh. 
    if (err_est_rel < ERR_STOP) 
    {
      printf("\nDone\n");
      break;
    }
    printf("  - adapting... "); fflush(stdout);
    hp.adapt(THRESHOLD);						// run the adaptivity algorithm
    printf("done in %lf secs (refined %d element(s))\n", hp.get_adapt_time(), hp.get_num_refined_elements());

    if (rndof >= NDOF_STOP)
    {
      printf("\nDone.\n");
      break;
    }

    // Clean up.
    delete rspace;
    delete matrix;
    delete rhs;
    delete rmatrix;
    delete rrhs;

    // Next iteration.
    as++;

  } while (!done);

  return 1;
}
