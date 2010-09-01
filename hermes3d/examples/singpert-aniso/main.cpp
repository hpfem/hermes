#include "config.h"
#ifdef USE_PETSC
#include <petsc.h>
#endif
#ifdef USE_UMFPACK
#include <umfpack.h>
#endif
#include <getopt.h>
#include <hermes3d.h>

//  With large K, this is a singularly perturbed problem that exhibits an extremely
//  thin and steep boundary layer. Singularly perturbed problems are considered to
//  be very difficult, but you'll see that Hermes can solve them easily even for large
//  values of K.
//
//  PDE: -Laplace u + K*K*u = K*K.
//
//  Domain: cube (0, 1) x (0, 1) x (0, 1), see the file singpert-aniso.mesh3d.
//
//  BC:  Homogeneous Dirichlet.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 2;			// Number of initial uniform mesh refinements.
const int P_INIT = 1;				// Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;			// Error threshold for element refinement of the adapt(...) function 
						// (default) STRATEGY = 0 ... refine elements elements until sqrt(THRESHOLD) 
						// times total error is processed. If more elements have similar errors, 
						// refine all to keep the mesh symmetric.
						// STRATEGY = 1 ... refine all elements whose error is larger
						// than THRESHOLD times maximum element error.
const double ERR_STOP  = 1;			// Stopping criterion for adaptivity (rel. error tolerance between the
						// fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;			// Adaptivity process stops when the number of degrees of freedom grows
						// over this limit. This is to prevent h-adaptivity to go on forever.
bool do_output = true;				// generate output files (if true)

// Problem constants
const double K_squared = 1e4;			// Equation parameter.

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
template<typename real, typename scalar>
scalar biform(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, fn_t<real> *v, geom_t<real> *e, user_data_t<scalar> *data)
{
  return int_grad_u_grad_v<real, scalar>(n, wt, u, v, e) + K_squared * int_u_v<real, scalar>(n, wt, u, v, e);
}

template<typename real, typename scalar>
scalar liform(int n, double *wt, fn_t<scalar> *u_ext[], fn_t<real> *u, geom_t<real> *e, user_data_t<scalar> *data)
{
  return K_squared * int_u<real, scalar>(n, wt, u, e);
}

// Output the element order.
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
    warning("Could not open file '%s' for writing.", fname);
}

/***********************************************************************************
 * main program                                                                    *
 ***********************************************************************************/
int main(int argc, char **args) {
#ifdef WITH_PETSC
  PetscInitialize(NULL, NULL, PETSC_NULL, PETSC_NULL);
  PetscPushErrorHandler(PetscIgnoreErrorHandler, PETSC_NULL);		// disable PETSc error handler
#endif

  // Load the inital mesh.
  Mesh mesh;
  Mesh3DReader mesh_loader;
  mesh_loader.load("singpert-aniso.mesh3d", &mesh);

  // Perform initial mesh refinements.
  printf("Performing %d initial mesh refinements.\n", INIT_REF_NUM);
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
  Word_t (nelem) = mesh.get_num_elements();
  printf("New number of elements is %d.\n", (int) nelem);

  // Initialize the shapset and the cache.
  H1ShapesetLobattoHex shapeset;

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
  wf.add_matrix_form(biform<double, double>, biform<ord_t, ord_t>, SYM, ANY);
  wf.add_vector_form(liform<double, double>, liform<ord_t, ord_t>, ANY);

  // Initialize the coarse mesh problem.
  LinearProblem lp(&wf);
  lp.set_space(&space);

  // Adaptivity loop.
  int as = 0; bool done = false;
  do {
    printf("\n---- Adaptivity step %d:\n", as);
    printf("\nSolving on coarse mesh:\n");

    // Procedures for coarse mesh problem.
    // Assign DOF.
    int ndof = space.assign_dofs();
    printf("  - Number of DOF: %d\n", ndof);


    // Assemble stiffness matrix and rhs.
    printf("  - Assembling... "); fflush(stdout);
    if (lp.assemble(&mat, &rhs))
      printf("done in %lf secs\n", lp.get_time());
    else
      error("failed!");

    // Solve the system.
    printf("  - Solving... "); fflush(stdout);
    bool solved = solver.solve();
    if (solved)
      printf("done in %lf secs\n", solver.get_time());
    else {
      printf("Failed\n");
      break;
    }

    // Construct a solution.
    Solution sln(&mesh);
    sln.set_fe_solution(&space, solver.get_solution());

    // Output the orders and the solution.
    if (do_output) 
    {
      out_orders(&space, "order", as);
      out_fn(&sln, "sln", as);
    }

    // Solving the fine mesh problem.
    printf("Solving on fine mesh:\n");

    // Matrix solver. 
#if defined WITH_UMFPACK
    UMFPackLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_PETSC
    PetscLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_MUMPS
    MumpsSolver rsolver(&mat, &rhs);
#endif

    // Construct the refined mesh for reference(refined) solution. 
    Mesh rmesh;
    rmesh.copy(mesh);
    rmesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

    // Setup space for the reference solution.
    Space *rspace = space.dup(&rmesh);
    rspace->copy_orders(space, 1);

    // Initialize the mesh problem for reference solution.
    LinearProblem rlp(&wf);
    rlp.set_space(rspace);

    // Assign DOF.
    int rndof = rspace->assign_dofs();
    printf("  - Number of DOF: %d\n", rndof);

    // Assemble stiffness matric and rhs.
    printf("  - Assembling... "); fflush(stdout);
    if (rlp.assemble(&mat, &rhs))
      printf("done in %lf secs\n", rlp.get_time());
    else
      error("failed!");

    // Solve the system.
    printf("  - Solving... "); fflush(stdout);
    bool rsolved = rsolver.solve();
    if (rsolved)
      printf("done in %lf secs\n", rsolver.get_time());
    else {
      printf("failed\n");
      break;
    }

    // Construct the reference solution.
    Solution rsln(&rmesh);
    rsln.set_fe_solution(rspace, rsolver.get_solution());

    // Compare coarse and fine mesh.
    // Calculate the error estimate wrt. refined mesh solution.
    double err = h1_error(&sln, &rsln);
    printf("  - H1 error: % lf\n", err * 100);

    // Save it to the graph.
    graph.add_value(0, ndof, err * 100);
    if (do_output)
      graph.save("conv.gp");

    // Calculate error estimates for hp-adaptivity.
    printf("Adaptivity\n");
    printf("  - calculating error: "); fflush(stdout);
    H1Adapt hp(&space);
    hp.set_aniso(true);							// anisotropic adaptivity.
    double err_est = hp.calc_error(&sln, &rsln) * 100;
    printf("% lf %%\n", err_est);

    // If error is too large, adapt the mesh. 
    if (err_est < ERR_STOP) 
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

    // Next iteration.
    as++;

    mat.free();
    rhs.free();
  } while (!done);

#ifdef WITH_PETSC
  PetscFinalize();
#endif

  return 1;
}
