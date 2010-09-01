#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

// This example shows how to solve linear elasticity problem on a 
// L-shape beam. Since we are dealing with only the isotropic 
// homogeneous case, displacement components in x-, y- and z-direction 
// have similar qualitative behavior. 
//
// PDE: Lame equations of linear elasticity.
//
// BC: u_x = u_y = 0 
//     du_z/dn = f on top of the L-beam
//     du_x/dn = du_y/dn = du_z/dn = 0 elsewhere
//
// The following parameters can be changed:
const int INIT_REF_NUM = 3;					// Number of initial uniform mesh refinements.
const int P_INIT = 4;						// Initial polynomial degree of all mesh elements.

// Problem parameters. 
const double E  = 200e9; 		// Young modulus for steel: 200 GPa.
const double nu = 0.3;			// Poisson ratio. 
const double f  = 1e4;   		// Load force (N). 
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));

// Boundary condition types. 
BCType bc_types_x(int marker) 
{
  return BC_NATURAL;
}

BCType bc_types_y(int marker) 
{
  return BC_NATURAL;
}

BCType bc_types_z(int marker) 
{
  return (marker == 3) ? BC_ESSENTIAL : BC_NATURAL;
}

// Output the solutions.
void out_fn(MeshFunction *x, MeshFunction *y, MeshFunction *z, const char *name) 
{
  char fname[1024];
  sprintf(fname, "%s.vtk", name);
  FILE *ofile = fopen(fname, "w");
  if (ofile != NULL) {
    VtkOutputEngine output(ofile);
    output.out(x, y, z, name);
    fclose(ofile);
  }
  else 
    warning("Can not open '%s' for writing.", fname);
}

#include "forms.cpp"

/***********************************************************************************
 * main program                                                                    *
 ***********************************************************************************/
int main(int argc, char **argv) {
#ifdef WITH_PETSC
  PetscInitialize(&argc, &argv, (char *) PETSC_NULL, PETSC_NULL);
  PetscPushErrorHandler(PetscIgnoreErrorHandler, PETSC_NULL);                   // Disable PETSc error handler.
#endif

  // Load the initial mesh. 
  Mesh mesh;
  Mesh3DReader mloader;
  mloader.load("l-beam.mesh3d", &mesh);

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

  // Create H1 spaces x-displacement component. 
  H1Space xdisp(&mesh, &shapeset);
  xdisp.set_bc_types(bc_types_x);
  xdisp.set_uniform_order(order3_t(P_INIT, P_INIT, P_INIT));

  // Create H1 spaces y-displacement component. 
  H1Space ydisp(&mesh, &shapeset);
  ydisp.set_bc_types(bc_types_y);
  ydisp.set_uniform_order(order3_t(P_INIT, P_INIT, P_INIT));

  // Create H1 spaces z-displacement component. 
  H1Space zdisp(&mesh, &shapeset);
  zdisp.set_bc_types(bc_types_z);
  zdisp.set_uniform_order(order3_t(P_INIT, P_INIT, P_INIT));

  // Assign DOF. 
  int ndofs = 0;
  ndofs += xdisp.assign_dofs(ndofs);
  ndofs += ydisp.assign_dofs(ndofs);
  ndofs += zdisp.assign_dofs(ndofs);
  printf("  - Number of DOFs: %d\n", ndofs);

  // Initialized the Weak formulation.
  WeakForm wf(3);
  wf.add_matrix_form(0, 0, bilinear_form_0_0<double, scalar>, bilinear_form_0_0<ord_t, ord_t>, SYM);
  wf.add_matrix_form(0, 1, bilinear_form_0_1<double, scalar>, bilinear_form_0_1<ord_t, ord_t>, SYM);
  wf.add_matrix_form(0, 2, bilinear_form_0_2<double, scalar>, bilinear_form_0_2<ord_t, ord_t>, SYM);
  wf.add_vector_form_surf(0, surf_linear_form_0<double, scalar>, surf_linear_form_0<ord_t, ord_t>);

  wf.add_matrix_form(1, 1, bilinear_form_1_1<double, scalar>, bilinear_form_1_1<ord_t, ord_t>, SYM);
  wf.add_matrix_form(1, 2, bilinear_form_1_2<double, scalar>, bilinear_form_1_2<ord_t, ord_t>, SYM);
  wf.add_vector_form_surf(1, surf_linear_form_1<double, scalar>, surf_linear_form_1<ord_t, ord_t>);

  wf.add_matrix_form(2, 2, bilinear_form_2_2<double, scalar>, bilinear_form_2_2<ord_t, ord_t>, SYM);
  wf.add_vector_form_surf(2, surf_linear_form_2<double, scalar>, surf_linear_form_2<ord_t, ord_t>, 5);

  // Initialize the mesh problem.
  LinearProblem lp(&wf);
  lp.set_spaces(Tuple<Space *>(&xdisp, &ydisp, &zdisp));

  // Assemble stiffness matrix
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
    printf("failed\n");
  }

  // Construct a solution. 
  double *s = solver.get_solution();
  Solution xsln(&mesh), ysln(&mesh), zsln(&mesh);
  xsln.set_fe_solution(&xdisp, s);
  ysln.set_fe_solution(&ydisp, s);
  zsln.set_fe_solution(&zdisp, s);

  // Output the solutions. 
  printf("  - Output... "); fflush(stdout);
  out_fn(&xsln, &ysln, &zsln, "disp");
  printf("done\n");

#ifdef WITH_PETSC
  mat.free();
  rhs.free();
  PetscFinalize();
#endif

  return 1;
}
