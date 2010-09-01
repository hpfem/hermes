#include "config.h"
#ifdef USE_PETSC
#include <petsc.h>
#endif
#ifdef USE_UMFPACK
#include <umfpack.h>
#endif
#include <getopt.h>
#include <hermes3d.h>

//  This benchmark solves the Poisson equation with an exact solution from fichera corner. 
//  The problem geometry is a cube with missing corner. The exact solution has a singular 
//  gradient (an analogy of infinite stress) at the center. The knowledge of the exact 
//  solution makes it possible to calculate the approximation error exactly.  You can 
//  compare h- and hp-adaptivity from the point of view of both CPU time requirements and 
//  discrete problem size. Also look at the quality of the a-posteriori error estimator 
//  used by Hermes. 
//
//  PDE: -Laplace u = f.
//
//  Known exact solution, see functions fn() and fndd().
//
//  Domain: unit cube (0, 0, 1)x(0, 1, 0)x(1, 0, 0), with a missing corner, see the file "fichera-corner.mesh3d".
//
//  BC:  Essential Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 2;          // Number of initial uniform mesh refinements.
const int P_INIT = 2;                // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;        // Error threshold for element refinement of the adapt(...) function 
                                     // (default) STRATEGY = 0 ... refine elements elements until sqrt(THRESHOLD) 
                                     // times total error is processed. If more elements have similar errors, 
                                     // refine all to keep the mesh symmetric.
                                     // STRATEGY = 1 ... refine all elements whose error is larger
                                     // than THRESHOLD times maximum element error.
const double ERR_STOP  = 1;          // Stopping criterion for adaptivity (rel. error tolerance between the
                                     // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;        // Adaptivity process stops when the number of degrees of freedom grows
                                     // over this limit. This is to prevent h-adaptivity to go on forever.
bool do_output = true;               // Generate output files (if true).


// Exact solution. 
#include "exact_solution.cpp"

// Boundary condition types.
BCType bc_types(int marker) 
{
  return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z) 
{
  return fn(x, y, z);
}

// Weak forms.
#include "forms.cpp"

// Output element orders.
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

// Output solutions. 
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
  else
    warning("Could not open file '%s' for writing.", fname);
}

/***********************************************************************************
 * main program                                                                    *
 ***********************************************************************************/
int main(int argc, char **args) 
{
#ifdef WITH_PETSC
  PetscInitialize(NULL, NULL, PETSC_NULL, PETSC_NULL);
  PetscPushErrorHandler(PetscIgnoreErrorHandler, PETSC_NULL);		// disable PETSc error handler
#endif

  // Load the inital mesh.
  Mesh mesh;
  Mesh3DReader mesh_loader;
  mesh_loader.load("fichera-corner.mesh3d", &mesh);

  // Initial uniform mesh refinements.
  printf("Performing %d initial mesh refinements.\n", INIT_REF_NUM);
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
  Word_t (nelem) = mesh.get_num_elements();
  printf("New number of elements is %d.\n", (int) nelem);

  // Initialize the shapset and the cache. 
  H1ShapesetLobattoHex shapeset;

  // Matrix solver.
#if defined WITH_UMFPACK
  UMFPackMatrix mat;
  UMFPackVector rhs;
  UMFPackLinearSolver solver(&mat, &rhs);
#elif defined WITH_PETSC
  PetscMatrix mat;
  PetscVector rhs;
  PetscLinearSolver solver(&mat, &rhs);
#elif defined WITH_MUMPS
  MumpsMatrix mat;
  MumpsVector rhs;
  MumpsSolver solver(&mat, &rhs);
#endif

  // Graphs of DOF convergence.
  GnuplotGraph graph;
  graph.set_captions("", "Degrees of Freedom", "Error [%]");
  graph.set_log_y();
  graph.add_row("Total error", "k", "-", "O");

  // Create H1 space to setup the problem. 
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_essential_bc_values(essential_bc_values);
  space.set_uniform_order(order3_t(P_INIT, P_INIT, P_INIT)); 

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, double>, bilinear_form<ord_t, ord_t>, SYM, ANY);
  wf.add_vector_form(linear_form<double, double>, linear_form<ord_t, ord_t>, ANY);

  // Initialize the coarse mesh problem. 
  LinProblem lp(&wf);
  lp.set_space(&space);

  // Adaptivity loop. 
  int as = 0; bool done = false;
  do 
  {
    printf("\n---- Adaptivity step #%d:\n", as);

    printf("\nSolving on coarse mesh:\n");

    // Procedure for coarse mesh problem. 
    // Assign DOF. 
    int ndof = space.assign_dofs();
    printf("  - Number of DOF: %d\n", ndof);

    // Initialize timer.
    Timer t_assemble;
    Timer t_solver;

    // Assemble stiffness matrix and rhs.
    printf("  - Assembling... "); fflush(stdout);
    t_assemble.start();
    bool assembled = lp.assemble(&mat, &rhs);
    t_assemble.stop();
    if (assembled)
      printf("done in %s (%lf secs)\n", t_assemble.get_human_time(), t_assemble.get_seconds());
    else
      error("failed!");

    // Solve the system. 
    printf("  - Solving... "); fflush(stdout);
    t_solver.start();
    bool solved = solver.solve();
    t_solver.stop();
    if (solved)
      printf("done in %s (%lf secs)\n", t_solver.get_human_time(), t_solver.get_seconds());
    else {
      printf("Failed\n");
      break;
    }

    // Construct a solution. 
    Solution sln(&mesh);
    sln.set_fe_solution(&space, solver.get_solution());

    // Output the orders and the solution.
    if (do_output) {
      out_orders(&space, "order", as);
      out_fn(&sln, "sln", as);
    }

    // Reference solution.
    printf("Reference solution\n");

#if defined WITH_UMFPACK
    UMFPackLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_PETSC
    PetscLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_MUMPS
    MumpsSolver rsolver(&mat, &rhs);
#endif

    // Construct the mesh for reference solution.
    Mesh rmesh;
    rmesh.copy(mesh);
    rmesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

    // Setup space for the reference (globally refined)  solution.
    Space *rspace = space.dup(&rmesh);
    rspace->copy_orders(space, 1);

    // Initialize the mesh problem for reference solution. 
    LinProblem rlp(&wf);
    rlp.set_space(rspace);

    // Assign DOF. 
    int rndof = rspace->assign_dofs();
    printf("  - Number of DOF: %d\n", rndof);

    // Assemble stiffness matric and rhs. 
    printf("  - Assembling... "); fflush(stdout);
    t_assemble.start();
    assembled = rlp.assemble(&mat, &rhs);
    t_assemble.stop();
    if (assembled)
      printf("done in %s (%lf secs)\n", t_assemble.get_human_time(), t_assemble.get_seconds());
    else
      error("failed!");

    // Solve the system.
    printf("  - Solving... "); fflush(stdout);
    t_solver.start();
    bool rsolved = rsolver.solve();
    t_solver.stop();
    if (rsolved)
      printf("done in %s (%lf secs)\n", t_solver.get_human_time(), t_solver.get_seconds());
    else {
      printf("failed\n");
      break;
    }

    // Construct the reference solution.
    Solution rsln(&rmesh);
    rsln.set_fe_solution(rspace, rsolver.get_solution());

    // Calculate the error estimate.
    double err = h1_error(&sln, &rsln);
    printf("  - H1 error: % lf\n", err * 100);

    // Save it to the graph. 
    graph.add_value(0, ndof, err * 100);
    if (do_output)
    graph.save("conv.gp");

    // Calculate error estimates for adaptivity. 
    printf("Adaptivity\n");
    printf("  - calculating error: "); fflush(stdout);
    H1Adapt hp(&space);
    Timer t_err;
    t_err.start();
    double err_est = hp.calc_error(&sln, &rsln) * 100;		// Calc error estimates on elements.
    t_err.stop();
    printf("% lf\n", err_est);

    // If error is small, then we can stop. 
    if (err_est < ERR_STOP) {
      printf("\nDone\n");
      break;
    }

    // If error is too large, do the hp-adaptivity.
    printf("  - adapting... "); fflush(stdout);
    Timer t_adapt;
    t_adapt.start();
    hp.adapt(THRESHOLD);					// run the adaptivity algorithm
    t_adapt.stop();
    printf("done in %lf secs (refined %d element(s))\n", t_adapt.get_seconds(), hp.get_num_refined_elements());

    if (rndof >= NDOF_STOP)
    {
      printf("\nDone.\n");
      break;
    }

    // Clean up.
    delete rspace;
    mat.free();
    rhs.free();

    // Next iteration.
    as++;
  } while (!done);

#ifdef WITH_PETSC
  PetscFinalize();
#endif

  return 1;
}
