#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This is the nineth in the series of NIST benchmarks with known exact solutions. This benchmark
//  has four different versions, use the global variable PROB_PARAM below to switch among them.
//
//  Reference: W. Mitchell, A Collection of 2D Elliptic Problems for Testing Adaptive Algorithms, 
//                          NIST Report 7668, February 2010.
//
//  PDE: -Laplace u = f
//
//  Known exact solution; atan(ALPHA * (sqrt(pow(x - X_LOC, 2) + pow(y - Y_LOC, 2)) - R_ZERO));
//  See functions fn() and fndd() in "exact_solution.cpp".
//
//  Domain: unit square (0, 1)x(0, 1), see the file square.mesh.
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

int PROB_PARAM = 3;    // PROB_PARAM determines which parameter values you wish to use for the steepness and location of the wave front. 
                       //    name		ALPHA	X_LOC	Y_LOC	R_ZERO
                       // 0: mild		20	-0.05	-0.05	0.7
                       // 1: steep		1000	-0.05	-0.05	0.7
                       // 2: asymmetric         1000	 1.5	 0.25	0.92
                       // 3: well		50	 0.5	 0.5	0.25

const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
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
const CandList CAND_LIST = H2D_HP_ANISO_H;        // Predefined list of element refinement candidates. Possible values are
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
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
double ALPHA;          // (X_LOC, Y_LOC) is the center of the circular wave front, R_ZERO is the distance from the 
double X_LOC;          // wave front to the center of the circle, and ALPHA gives the steepness of the wave front.
double Y_LOC;
double R_ZERO;

// Exact solution.
#include "exact_solution.cpp"

// Boundary condition types.
BCType bc_types(int marker) { return BC_ESSENTIAL;}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y) { return fn(x, y);}

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

  // Perform initial mesh refinements.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form), H2D_SYM);
  wf.add_vector_form(callback(linear_form));

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize adaptivity parameters.
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY, 
                          MESH_REGULARITY);

  // Adaptivity loop.
  Solution *sln = new Solution();
  Solution *ref_sln = new Solution();
  ExactSolution exact(&mesh, fndd);
  WinGeom* sln_win_geom = new WinGeom(0, 0, 440, 350);
  WinGeom* mesh_win_geom = new WinGeom(450, 0, 400, 350);
  bool verbose = true;     // Print info during adaptivity.
  solve_linear_adapt(&space, &wf, NULL, matrix_solver, H2D_H1_NORM, sln, ref_sln, 
                     sln_win_geom, mesh_win_geom, &selector, &apt, verbose, &exact);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
