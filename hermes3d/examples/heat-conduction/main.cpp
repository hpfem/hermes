#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#ifdef WITH_UMFPACK
#include <umfpack.h>
#endif
#include <getopt.h>
#include <hermes3d.h>
//  This example shows how to solve a time-dependent PDE discretized
//  in time via the implicit Euler method on a fixed mesh. 
//  You will also learn how to use the
//  solution calculated in the previous time step.
//
//  PDE: stationary heat transfer equation
//  dT/dt - Laplace T = f.
//
//  Domain: (0, 0, 1)x(0, 1, 0)x(1, 0, 0), see the file hexahedron.mesh3d. 
//
//  Known exact solution, see unctions fn() and fndd(). 
//
//  IC:  T = 0.
//  BC:  T = 0. ... Dirichlet,
//
//  Time-stepping: implicit Euler.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 2;					// Number of initial uniform mesh refinements.
const int P_INIT = 4;						// Initial polynomial degree of all mesh elements.
const double TAU = 0.05;					// Time step in seconds. 
bool do_output = true;						// Generate output files (if true).

// Problem parameters. 
const double FINAL_TIME = 2 * M_PI;				// Length of time inveral in seconds. 

// Global time variable. 
double TIME = TAU;

// Exact solution. 
#include "exact_solution.cpp"

// Boundary condition types.
BCType bc_types(int marker) {
  return BC_ESSENTIAL;
}

// Weak forms. 
#include "forms.cpp"

// Output the solutions. 
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
    warning("Can not open '%s' for writing.", fname);
}

/***********************************************************************************
 * main program                                                                    *
 ***********************************************************************************/
int main(int argc, char **argv) {
#ifdef WITH_PETSC
  PetscInitialize(NULL, NULL, (char *) PETSC_NULL, PETSC_NULL);
  PetscPushErrorHandler(PetscIgnoreErrorHandler, PETSC_NULL);			// Disable PETSc error handler.
#endif

  // Load the initial mesh. 
  Mesh mesh;
  Mesh3DReader mesh_loader;
  mesh_loader.load("hexahedron.mesh3d", &mesh);

  // Initial uniform mesh refinements. 
  printf("Performing %d initial mesh refinements.\n", INIT_REF_NUM);
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
  Word_t (nelem) = mesh.get_num_elements();
  printf("New number of elements is %d.\n", (int) nelem);

  // Initialize the shapeset and the cache. 
  H1ShapesetLobattoHex shapeset;

#if defined WITH_UMFPACK
  UMFPackMatrix mat;
  UMFPackVector rhs;
  UMFPackLinearSolver solver(&mat, &rhs);
#elif defined WITH_PARDISO
  PardisoMatrix mat;
  PardisoVector rhs;
  PardisoLinearSolver solver(&mat, &rhs);
#elif defined WITH_PETSC
  PetscMatrix mat;
  PetscVector rhs;
  PetscLinearSolver solver(&mat, &rhs);
#elif defined WITH_MUMPS
  MumpsMatrix mat;
  MumpsVector rhs;
  MumpsSolver solver(&mat, &rhs);
#endif

  // Create H1 space to setup the problem. 
  H1Space space(&mesh, &shapeset);
  space.set_bc_types(bc_types);
  space.set_uniform_order(order3_t(P_INIT, P_INIT, P_INIT));

  // Assign DOF.
  int ndof = space.assign_dofs();
  printf("  - Number of DOFs: %d\n", ndof);

  // Construct initial solution and set zero.
  Solution sln_prev(&mesh);
  sln_prev.set_zero();

  // Initialize the weak formulation. 
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM);
  wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>, ANY, &sln_prev);

  // Initialize the coarse mesh problem. 
  LinearProblem lp(&wf);
  lp.set_spaces(&space);

  // Time stepping. 
  int nsteps = (int) (FINAL_TIME/TAU + 0.5);
  for (int ts = 0; ts < nsteps;  ts++)
  {
    printf("\n---- Time step %d, time %3.5f\n", ts, TIME);

    // Assemble stiffness matrix and rhs. 
    printf("  - Assembling... "); fflush(stdout);
    Timer tmr_assemble;
    tmr_assemble.start();
    bool assembled = lp.assemble(&mat, &rhs);
    tmr_assemble.stop();
    if (assembled)
      printf("done in %s (%lf secs)\n", tmr_assemble.get_human_time(), tmr_assemble.get_seconds());
    else
      error("failed!");

    // Solve the stiffness matrix. 
    printf("  - Solving... "); fflush(stdout);
    Timer tmr_solve;
    tmr_solve.start();
    bool solved = solver.solve();
    tmr_solve.stop();

    if (solved) 
      printf("done in %s (%lf secs)\n", tmr_solve.get_human_time(), tmr_solve.get_seconds());
    else {
      error("failed\n");
      break;
    }

    // Construct a solution. 
    Solution sln(&mesh);
    sln.set_fe_solution(&space, solver.get_solution());

    // Output the solution. 
    if (do_output)
    {
      out_fn(&sln, "sln", ts);
    }

    // Calculate error wrt. exact solution. 
    printf("  - Calculating exact error ...\n"); fflush(stdout);
    ExactSolution esln(&mesh, fndd);
    double err_exact = h1_error(&sln, &esln) * 100; 
    printf("  - err. exact: %le\n", err_exact);

    // next step
    sln_prev = sln;
    TIME += TAU;
  }

#ifdef WITH_PETSC
  mat.free();
  rhs.free();
  PetscFinalize();
#endif

  return 1;
}
