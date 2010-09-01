#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>
#include <float.h>
#include <getopt.h>

//  This example comes with an exact solution, and it describes the diffraction
//  of an electromagnetic wave from a re-entrant corner. Convergence graphs saved
//  (both exact error and error estimate, and both wrt. dof number and cpu time).
//
//  PDE: time-harmonic Maxwell's equations
//
//  Known exact solution, see functions exact_sol_val(), exact_sol(), exact()
//
//  Domain: L-shape 3D domain
//
//  Meshes: "lshape_hex.mesh3d" (hexahedron mesh) See the mesh.load(...) command below.
//
//  BC: perfect conductor on boundary markers 1 and 6 (essential BC)
//      impedance boundary condition on rest of boundary (natural BC)
//
//  The following parameters can be changed:
const int P_INIT = 1;				// Initial polynomial degree. NOTE: The meaning is different from
						// standard continuous elements in the space H1. Here, P_INIT refers
						// to the maximum poly order of the tangential component, and polynomials
						// of degree P_INIT + 1 are present in element interiors. P_INIT = 0
						// is for Whitney elements.
bool do_output = true;				// generate output files (if true)

// Problem parameters.
const double mu_r   = 1.0;
const double kappa  = 1.0;
const double lambda = 1.0;

// Bessel functions, exact solution, and weak forms.
#include "forms.cpp"

// Boundary condition types. 
BCType bc_types(int marker)
{
  if (marker == 1 || marker == 6)
    return BC_ESSENTIAL; // perfect conductor
  else
    return BC_NATURAL; // impedance
}

// Helper function to output Mesh. 
void out_mesh(Mesh *mesh, const char *name)
{
  char fname[1024];
  sprintf(fname, "%s.vtk", name);
  FILE *f = fopen(fname, "w");
  if (f != NULL) {
    VtkOutputEngine vtk(f);
    vtk.out(mesh);
    fclose(f);
  }
  else
    warning("Could not open file '%s' for writing.", fname);
}

// Helper function to output solution. 
void out_fn(MeshFunction *fn, const char *name)
{
  char of_name[1024];
  FILE *ofile;
  // mesh out
  sprintf(of_name, "%s.vtk", name);
  ofile = fopen(of_name, "w");
  if (ofile != NULL) {
    VtkOutputEngine output(ofile);
    output.out(fn, name);
    fclose(ofile);
  }
  else 
    warning("Can not not open '%s' for writing.", of_name);
}

/***********************************************************************************
 * main program                                                                    *
 ***********************************************************************************/
int main(int argc, char **argv)
{

#ifdef WITH_PETSC
  PetscInitialize(NULL, NULL, (char *) PETSC_NULL, PETSC_NULL);
  PetscPushErrorHandler(PetscIgnoreErrorHandler, PETSC_NULL);           // disable PETSc error handler
#endif

  // Load the mesh. 
  Mesh mesh;
  Mesh3DReader mloader;
  mloader.load("lshape_hex.mesh3d", &mesh);			// hexahedron

  // Perform initial mesh refinements. 
  mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
  mesh.refine_all_elements(H3D_H3D_REFT_HEX_XY);
  mesh.refine_all_elements(H3D_H3D_REFT_HEX_XY);
  mesh.refine_all_elements(H3D_H3D_REFT_HEX_XY);

  // Matrix solver. 
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

  // Create an Hcurl space with detault shapeset.
  HcurlShapesetLobattoHex shapeset;

  // Create Hcurl space to setup the problem.
  HcurlSpace sp(&mesh, &shapeset);
  sp.set_bc_types(bc_types);
  sp.set_uniform_order(order3_t(P_INIT, P_INIT, P_INIT));

  // Assign DOF. 
  int ndof = sp.assign_dofs();
  printf("  - Number of DOFs: %d\n", ndof);

  //  Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(biform<double, scalar>, biform<ord_t, ord_t>, SYM);
  wf.add_matrix_form_surf(biform_surf, biform_surf_ord);
  wf.add_vector_form_surf(liform_surf, liform_surf_ord);

  // Initialize the coarse mesh problem. 
  LinProblem lp(&wf);
  lp.set_space(&sp);

  // Assemble stiffness matrix.
  printf("  - assembling... "); fflush(stdout);
  Timer tmr_assemble;
  tmr_assemble.start();
  lp.assemble(&mat, &rhs);
  tmr_assemble.stop();
  printf("done in %s (%lf secs)\n", tmr_assemble.get_human_time(), tmr_assemble.get_seconds());

  // Solve the stiffness matrix.
  printf("  - solving... "); fflush(stdout);
  Timer tmr_solve;
  tmr_solve.start();
  bool solved = solver.solve();
  tmr_solve.stop();

  // Print solving time. 
  if (solved) {
    printf("done in %s (%lf secs)\n", tmr_solve.get_human_time(), tmr_solve.get_seconds());
  }
  else {
    printf("failed\n");
    return -1;
  }

  // Construct the solution.
  std::complex<double> *s = solver.get_solution();
  Solution sln(&mesh);
  sln.set_fe_solution(&sp, s);

  // Output the solution.
  if (do_output) {
    printf("  - output... "); fflush(stdout);
    out_fn(&sln, "solution");
    out_mesh(&mesh, "mesh");
    printf("done\n");
  }

#ifdef WITH_PETSC
  mat.free();
  rhs.free();
  PetscFinalize();
#endif

  return 0;
}
