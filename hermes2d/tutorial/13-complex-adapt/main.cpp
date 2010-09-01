#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This problem describes the distribution of the vector potential in
//  a 2D domain comprising a wire carrying electrical current, air, and
//  an iron which is not under voltage.
//
//  PDE: -Laplace A + ii*omega*gamma*mu*A = mu *J_ext.
//
//  Domain: Rectangle of height 0.003 and width 0.004. Different
//  materials for the wire, air, and iron (see mesh file domain2.mesh).
//
//  BC: Zero Dirichlet on the top and right edges, zero Neumann
//  elsewhere.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 1;              // Number of initial uniform mesh refinements.
const int P_INIT = 1;                    // Initial polynomial degree of all mesh elements.
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
const double ERR_STOP = 1.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.


// Problem parameters.
double mu_0 = 4.0*3.141592654E-7;
double J_wire = 5000000.0;
double freq = 5E3;
double omega = 2*3.141592654*freq;
double gamma_iron = 6E6;
double mu_iron = 1000*mu_0;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker==1) {return BC_NATURAL;}
  if (marker==2) {return BC_ESSENTIAL;}
  if (marker==3) {return BC_ESSENTIAL;}
  if (marker==4) {return BC_ESSENTIAL;}
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return cplx(0.0,0.0);
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
  mloader.load("domain2.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form_iron), H2D_SYM, 3);
  wf.add_matrix_form(callback(bilinear_form_wire), H2D_SYM, 2);
  wf.add_matrix_form(callback(bilinear_form_air), H2D_SYM, 1);
  wf.add_vector_form(callback(linear_form_wire), 2);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize adaptivity parameters.
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY, 
                          MESH_REGULARITY);

  // Adaptivity loop.
  Solution *sln = new Solution();
  Solution *ref_sln = new Solution();
  WinGeom* sln_win_geom = new WinGeom(0, 0, 600, 350);
  WinGeom* mesh_win_geom = new WinGeom(610, 0, 520, 350);
  bool verbose = true;     // Print info during adaptivity.
  bool is_complex = true;
  // The NULL pointer means that we do not want the resulting coefficient vector.
  solve_linear_adapt(&space, &wf, NULL, matrix_solver, H2D_H1_NORM, sln, ref_sln, 
                     sln_win_geom, mesh_win_geom, &selector, &apt, verbose,
                     Tuple<ExactSolution *>(), is_complex); // Do not use NULL to pass an empty Tuple.

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
