#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example 22-newton-timedep-heat-adapt works correctly.

const int INIT_REF_NUM = 2;                // Number of initial uniform mesh refinements.
const int P_INIT = 2;                      // Initial polynomial degree of all mesh elements.
const int TIME_DISCR = 2;                  // 1 for implicit Euler, 2 for Crank-Nicolson.
const double TAU = 0.5;                    // Time step.
const double T_FINAL = 5.0;                // Time interval length.

// Adaptivity
const int UNREF_FREQ = 1;                  // Every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_H;  // Predefined list of element refinement candidates. Possible values are
                                           // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                           // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                           // See the User Documentation for details.
const int MESH_REGULARITY = -1;            // Maximum allowed level of hanging nodes:
                                           // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                           // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                           // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                           // Note that regular meshes are not supported, this is due to
                                           // their notoriously bad performance.
const double CONV_EXP = 1.0;               // Default value is 1.0. This parameter influences the selection of
                                           // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;               // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;     // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;       // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 100;            // Maximum allowed number of Newton iterations.

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u)
{
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) {
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions.
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// Initial condition.
scalar init_cond(double x, double y, double& dx, double& dy)
{
  return dir_lift(x, y, dx, dy);
}

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Heat sources (can be a general function of 'x' and 'y').
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Weak forms.
# include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Solutions.
  Solution sln, ref_sln, u_prev_time;

  // Initialize the weak formulation.
  WeakForm wf;
  if(TIME_DISCR == 1) {
    wf.add_matrix_form(callback(J_euler), H2D_UNSYM, H2D_ANY);
    wf.add_vector_form(callback(F_euler), H2D_ANY, &u_prev_time);
  }
  else {
    wf.add_matrix_form(callback(J_cranic), H2D_UNSYM, H2D_ANY);
    wf.add_vector_form(callback(F_cranic), H2D_ANY, &u_prev_time);
  }

  // Initialize adaptivity parameters.
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY, MESH_REGULARITY);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Assign initial condition to mesh.
  u_prev_time.set_exact(&mesh, init_cond);
  Vector *coeff_vec = new AVector();

  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    info("---- Time step %d:", ts);

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
    }

    // Update the coefficient vector and u_prev_time.
    info("Projecting to obtain coefficient vector on coarse mesh.");
    project_global(&space, H2D_H1_NORM, &u_prev_time, &u_prev_time, coeff_vec);

    // Adaptivity loop (in space):
    bool verbose = true;     // Print info during adaptivity.
    info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
    // The NULL pointers mean that we are not interested in visualization during the Newton's loop.
    solve_newton_adapt(&space, &wf, coeff_vec, matrix_solver, H2D_H1_NORM, &sln, &ref_sln, 
                       NULL, NULL, &selector, &apt, 
                       NEWTON_TOL_COARSE, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose);

    // Copy new time level reference solution into u_prev_time.
    u_prev_time.set_fe_solution(&space, coeff_vec);
  }

  int ndof = get_num_dofs(&space);

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  printf("ndof allowed = %d\n", 1100);
  printf("ndof actual = %d\n", ndof);
  if (ndof < 1100) {      // ndofs was 1038 at the time this test was created
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
