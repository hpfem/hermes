#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"


using namespace RefinementSelectors;

// This test makes sure that example "thermoelasticity" works correctly. 

const bool SOLVE_ON_COARSE_MESH = false; // If true, coarse mesh FE problem is solved in every adaptivity step.
                                         // If false, projection of the fine mesh solution on the coarse mesh is used. 
const int P_INIT_TEMP = 1;               // Initial polynomial degrees in temperature mesh.
const int P_INIT_DISP = 1;               // Initial polynomial degrees for displacement meshes.
const bool MULTI = true;                 // MULTI = true  ... use multi-mesh,
                                         // MULTI = false ... use single-mesh.
                                         // Note: In the single mesh option, the meshes are
                                         // forced to be geometrically the same but the
                                         // polynomial degrees can still vary.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 3.0;            // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
double HEAT_SRC = 10000.0;               // Heat source in the material (caused by induction heating).
double TEMP_INNER = 50;
double HEAT_FLUX_OUTER = -50;
const double E = 2e11;                   // Steel: E=200 GPa.
const double nu = 0.3;
const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));
const double l2m = lambda + 2*mu;
const double rho = 8000;
const double g = 9.81;
const double alpha = 13e-6;              // See http://hyperphysics.phy-astr.gsu.edu/hbase/tables/thexp.html.

//  Boundary markers.
const int marker_bottom = 1;
const int marker_sides = 2;
const int marker_top = 3;
const int marker_holes = 4;

// Boundary condition types.
BCType bc_types_x(int marker)
  { return (marker == marker_bottom) ? BC_ESSENTIAL : BC_NATURAL; }

BCType bc_types_y(int marker)
  { return (marker == marker_bottom) ? BC_ESSENTIAL : BC_NATURAL; }

BCType bc_types_t(int marker)
  { return (marker == marker_holes) ? BC_ESSENTIAL : BC_NATURAL; }

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values_temp(int ess_bdy_marker, double x, double y)
  { return TEMP_INNER; }

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh xmesh, ymesh, tmesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &xmesh); // Master mesh.

  // Initialize multimesh hp-FEM.
  ymesh.copy(&xmesh);                  // Ydisp will share master mesh with xdisp.
  tmesh.copy(&xmesh);                  // Temp will share master mesh with xdisp.

  // Create H1 spaces with default shapesets.
  H1Space xdisp(&xmesh, bc_types_x, NULL, P_INIT_DISP);
  H1Space ydisp(MULTI ? &ymesh : &xmesh, bc_types_y, NULL, P_INIT_DISP);
  H1Space temp(MULTI ? &tmesh : &xmesh, bc_types_t, essential_bc_values_temp, P_INIT_TEMP);

  // Initialize the weak formulation.
  WeakForm wf(3);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1), H2D_SYM);
  wf.add_matrix_form(0, 2, callback(bilinear_form_0_2));
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1));
  wf.add_matrix_form(1, 2, callback(bilinear_form_1_2));
  wf.add_matrix_form(2, 2, callback(bilinear_form_2_2));
  wf.add_vector_form(1, callback(linear_form_1));
  wf.add_vector_form(2, callback(linear_form_2));
  wf.add_vector_form_surf(2, callback(linear_form_surf_2));

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize adaptivity parameters.
  double to_be_processed = 0;
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY,
                          MESH_REGULARITY, to_be_processed, H2D_TOTAL_ERROR_REL, H2D_ELEMENT_ERROR_REL);

  apt.set_error_form(0, 0, bilinear_form_0_0<scalar, scalar>, bilinear_form_0_0<Ord, Ord>);
  apt.set_error_form(0, 1, bilinear_form_0_1<scalar, scalar>, bilinear_form_0_1<Ord, Ord>);
  apt.set_error_form(0, 2, bilinear_form_0_2<scalar, scalar>, bilinear_form_0_2<Ord, Ord>);
  apt.set_error_form(1, 0, bilinear_form_1_0<scalar, scalar>, bilinear_form_1_0<Ord, Ord>);
  apt.set_error_form(1, 1, bilinear_form_1_1<scalar, scalar>, bilinear_form_1_1<Ord, Ord>);
  apt.set_error_form(1, 2, bilinear_form_1_2<scalar, scalar>, bilinear_form_1_2<Ord, Ord>);
  apt.set_error_form(2, 2, bilinear_form_2_2<scalar, scalar>, bilinear_form_2_2<Ord, Ord>);

  // Adaptivity loop.
  Solution *xdisp_sln = new Solution();
  Solution *ydisp_sln = new Solution();
  Solution *temp_sln = new Solution();
  Solution *ref_xdisp_sln = new Solution();
  Solution *ref_ydisp_sln = new Solution();
  Solution *ref_temp_sln = new Solution();
  bool verbose = true;  // Print info during adaptivity.
  solve_linear_adapt(Tuple<Space *>(&xdisp, &ydisp, &temp), &wf, NULL, matrix_solver,
                     Tuple<int>(H2D_H1_NORM, H2D_H1_NORM, H2D_H1_NORM),
                     Tuple<Solution *>(xdisp_sln, ydisp_sln, temp_sln),
                     Tuple<Solution *>(ref_xdisp_sln, ref_ydisp_sln, ref_temp_sln),
                     Tuple<WinGeom *> (),
                     Tuple<WinGeom *> (),
                     Tuple<RefinementSelectors::Selector *> (&selector, &selector, &selector),
                     &apt, verbose);

  int ndof = get_num_dofs(Tuple<Space *>(&xdisp, &ydisp, &temp));

#define ERROR_SUCCESS       0
#define ERROR_FAILURE      -1
  int ndof_allowed = 2000;
  printf("ndof actual = %d\n", ndof);
  printf("ndof allowed = %d\n", ndof_allowed);
  if (ndof <= ndof_allowed) {      // ndofs was 1903 at the time this test was created
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
};

