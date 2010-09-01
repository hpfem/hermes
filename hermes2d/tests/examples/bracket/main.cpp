#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example "bracket" works correctly.

const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const bool MULTI = true;                          // MULTI = true  ... use multi-mesh,
                                                  // MULTI = false ... use single-mesh.
                                                  // Note: In the single mesh option, the meshes are
                                                  // forced to be geometrically the same but the
                                                  // polynomial degrees can still vary.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
const double E  = 200e9;                          // Young modulus for steel: 200 GPa.
const double nu = 0.3;                            // Poisson ratio.
const double f  = 1e3;                            // Load force: 10^3 N.
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));

// Boundary markers.
const int BDY_LEFT = 1;
const int BDY_TOP = 2;

// Boundary condition types.
BCType bc_types(int marker)
  { return (marker == BDY_LEFT) ? BC_ESSENTIAL : BC_NATURAL; }

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh u_mesh, v_mesh;
  H2DReader mloader;
  mloader.load("bracket.mesh", &u_mesh);

  // Initial mesh refinements.
  u_mesh.refine_element(1);
  u_mesh.refine_element(4);

  // Create initial mesh for the vertical displacement component.
  // This also initializes the multimesh hp-FEM.
  v_mesh.copy(&u_mesh);

  // Create H1 spaces with default shapesets.
  H1Space u_space(&u_mesh, bc_types, essential_bc_values, P_INIT);
  H1Space v_space(MULTI ? &v_mesh : &u_mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0), H2D_SYM);  // note that only one symmetric part is
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1), H2D_SYM);  // added in the case of symmetric bilinear
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1), H2D_SYM);  // forms
  wf.add_vector_form_surf(1, callback(linear_form_surf_1), BDY_TOP);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize adaptivity parameters.
  double to_be_processed = 0;
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY, 
                          MESH_REGULARITY, to_be_processed, H2D_TOTAL_ERROR_REL, H2D_ELEMENT_ERROR_REL);
  apt.set_error_form(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
  apt.set_error_form(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
  apt.set_error_form(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
  apt.set_error_form(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);

  // Adaptivity loop.
  Solution *u_sln = new Solution();
  Solution *v_sln = new Solution();
  Solution *ref_u_sln = new Solution();
  Solution *ref_v_sln = new Solution();
  bool verbose = true;  // Print info during adaptivity.
  solve_linear_adapt(Tuple<Space *>(&u_space, &v_space), &wf, NULL, matrix_solver,
                     Tuple<int>(H2D_H1_NORM, H2D_H1_NORM),
                     Tuple<Solution *>(u_sln, v_sln), Tuple<Solution *>(ref_u_sln, ref_v_sln),
                     Tuple<WinGeom *>(), Tuple<WinGeom *>(),
                     Tuple<RefinementSelectors::Selector *> (&selector, &selector),
                     &apt, verbose);

  int ndof = get_num_dofs(Tuple<Space *>(&u_space, &v_space));

#define ERROR_SUCCESS       0
#define ERROR_FAILURE      -1
  int ndof_allowed = 1450;
  printf("ndof actual = %d\n", ndof);
  printf("ndof allowed = %d\n", ndof_allowed);
  if (ndof <= ndof_allowed) {      // ndofs was 1414 at the time this test was created
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
