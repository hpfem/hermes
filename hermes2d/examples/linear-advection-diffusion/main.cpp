#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example solves a linear advection diffusion problem using optional
//  variational multiscale stabilization. To use the stabilization, you must
//  uncomment the definition of H2D_SECOND_DERIVATIVES_ENABLED in common.h
//  and rebuild this example. Note that in our experience, the stabilization
//  does only work for linear elements. Nevertheless, we were able to solve the
//  problem without stabilization using adaptive hp-FEM.
//
//  PDE: div(bu - \epsilon \nabla u) = 0 where b = (b1, b2) is a constant vector.
//
//  Domain: Square (0, 1)x(0, 1).
//
//  BC:  Dirichlet, see the function scalar essential_bc_values() below.
//
//  The following parameters can be changed:

const int P_INIT = 1;                    // Initial polynomial degree of all mesh elements.
const bool STABILIZATION_ON = false;     // Stabilization on/off (assumes that H2D_SECOND_DERIVATIVES_ENABLED is defined).
const bool SHOCK_CAPTURING_ON = false;   // Shock capturing on/off.
const int INIT_REF_NUM = 2;              // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 1;          // Number of initial uniform mesh refinements in the boundary layer region.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                  // Adaptive strategy:
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
const double ERR_STOP = 5.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
const double EPSILON = 0.01;             // Diffusivity.
const double B1 = 1., B2 = 1.;           // Advection direction, div(B) = 0.

// Boundary markers.
const int BDY_LAYER = 2;

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essemtial (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
    if (ess_bdy_marker == 1) return 1;
    else return 2 - pow(x, 0.1) - pow(y, 0.1);
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
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh);     // quadrilaterals
  // mloader.load("square_tri.mesh", &mesh);   // triangles

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_LAYER, INIT_REF_NUM_BDY);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  if (STABILIZATION_ON == true) {
    wf.add_matrix_form(callback(bilinear_form_stabilization));
  }
  if (SHOCK_CAPTURING_ON == true) {
    wf.add_matrix_form(callback(bilinear_form_shock_capturing));
  }

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
  solve_linear_adapt(&space, &wf, NULL, matrix_solver, H2D_H1_NORM, sln, ref_sln, 
                     sln_win_geom, mesh_win_geom, &selector, &apt, verbose);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
