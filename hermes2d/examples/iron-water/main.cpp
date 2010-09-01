#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example is a standard nuclear engineering benchmark describing an external-force-driven
//  configuration without fissile materials present, using one-group neutron diffusion approximation.
//  It is very similar to example "saphir", the main difference being that the mesh is loaded in
//  the ExodusII format (created for example by Cubit).
//
//  PDE: -div(D(x,y)grad\Phi) + \Sigma_a(x,y)\Phi = Q_{ext}(x,y)
//  where D(x, y) is the diffusion coefficient, \Sigma_a(x,y) the absorption cross-section,
//  and Q_{ext}(x,y) external sources.
//
//  Domain: square (0, L)x(0, L) where L = 30c (see mesh file domain.mesh).
//
//  BC:  Zero Dirichlet for the right and top edges ("vacuum boundary").
//       Zero Neumann for the left and bottom edges ("reflection boundary").
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.6;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
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
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
double L = 30;                           // Edge of square.
double L0 = 0.75*0.5*L;                  // End of first water layer.
double L1 = 0.5*L;                       // End of second water layer.
double L2 = 0.75*L;                      // End of iron layer.
double Q_EXT = 1.0;                      // Neutron source (nonzero in domain 1 only).
double SIGMA_T_WATER = 3.33;             // Total cross-section.
double SIGMA_T_IRON = 1.33;
double C_WATER = 0.994;                  // Scattering ratio.
double C_IRON = 0.831;
double D_WATER = 1./(3.*SIGMA_T_WATER);  // Diffusion coefficient.
double D_IRON = 1./(3.*SIGMA_T_IRON);
double SIGMA_A_WATER = SIGMA_T_WATER - C_WATER*SIGMA_T_WATER;  // Absorbing cross-section.
double SIGMA_A_IRON = SIGMA_T_IRON - C_IRON*SIGMA_T_IRON;

// Materials.
const int WATER_1 = 1;
const int WATER_2 = 2;
const int IRON = 3;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == 1) return BC_NATURAL;
  else return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0.0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  ExodusIIReader mloader;
  if (!mloader.load("iron-water.e", &mesh)) error("ExodusII mesh load failed.");

  // Perform initial uniform mesh refinement.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation
  WeakForm wf;
  wf.add_matrix_form(bilinear_form_water, bilinear_form_ord, H2D_SYM, WATER_1);
  wf.add_matrix_form(bilinear_form_water, bilinear_form_ord, H2D_SYM, WATER_2);
  wf.add_matrix_form(bilinear_form_iron, bilinear_form_ord, H2D_SYM, IRON);
  wf.add_vector_form(linear_form_source, linear_form_ord, WATER_1);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize adaptivity parameters.
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY, 
                          MESH_REGULARITY);

  // Adaptivity loop.
  Solution *sln = new Solution();
  Solution *ref_sln = new Solution();
  WinGeom* sln_win_geom = new WinGeom(0, 0, 440, 350);
  WinGeom* mesh_win_geom = new WinGeom(450, 0, 400, 350);
  bool verbose = true;     // Print info during adaptivity.
  // The NULL pointer means that we do not want the resulting coefficient vector.
  solve_linear_adapt(&space, &wf, NULL, matrix_solver, H2D_H1_NORM, sln, ref_sln, 
                     sln_win_geom, mesh_win_geom, &selector, &apt, verbose);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
