#define H2D_REPORT_INFO
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

//  This example shows how to solve a time-dependent PDE discretized
//  in time via the implicit Euler method. The St. Vitus Cathedral
//  in Prague (http://en.wikipedia.org/wiki/St._Vitus_Cathedral)
//  responds to changes in the surrounding air temperature
//  during one 24-hour cycle. You will also learn how to use the
//  solution calculated in the previous time step.
//
//  PDE: non-stationary heat transfer equation
//  HEATCAP * RHO * dT/dt - LAMBDA * Laplace T = 0.
//
//  Domain: St. Vitus cathedral (cathedral.mesh).
//
//  IC:  T = T_INIT.
//  BC:  T = T_INIT on the bottom edge ... Dirichlet,
//       dT/dn = ALPHA*(t_exterior(time) - T) ... Newton, time-dependent.
//
//  Time-stepping: implicit Euler.
//
//  The following parameters can be changed:

const int P_INIT = 4;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 1;                   // Number of initial uniform mesh refinements towards the boundary.
const double TAU = 300.0;                         // Time step in seconds.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
const double T_INIT = 10;        // Temperature of the ground (also initial temperature).
const double ALPHA = 10;         // Heat flux coefficient for Newton's boundary condition.
const double LAMBDA = 1e5;       // Thermal conductivity of the material.
const double HEATCAP = 1e6;      // Heat capacity.
const double RHO = 3000;         // Material density.
const double FINAL_TIME = 86400; // Length of time interval (24 hours) in seconds.

// Global time variable.
double TIME = 0;

// Time-dependent exterior temperature.
double temp_ext(double t) {
  return T_INIT + 10. * sin(2*M_PI*t/FINAL_TIME);
}

// Boundary markers.
int bdy_ground = 1;
int bdy_air = 2;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == bdy_ground) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return T_INIT;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("cathedral.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(bdy_air, INIT_REF_NUM_BDY);

  // Initialize an H1 space with default shepeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Set initial condition.
  Solution tsln;
  tsln.set_const(&mesh, T_INIT);

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, double>, bilinear_form<Ord, Ord>);
  wf.add_matrix_form_surf(bilinear_form_surf<double, double>, bilinear_form_surf<Ord, Ord>, bdy_air);
  wf.add_vector_form(linear_form<double, double>, linear_form<Ord, Ord>, H2D_ANY, &tsln);
  wf.add_vector_form_surf(linear_form_surf<double, double>, linear_form_surf<Ord, Ord>, bdy_air);

  // Initialize the linear problem.
  LinearProblem lp(&wf, &space);

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, ndof, mat, rhs, solver);

  // Initialize views.
  ScalarView Tview("Temperature", new WinGeom(0, 0, 450, 600));
  char title[100];
  sprintf(title, "Time %3.5f, exterior temperature %3.5f", TIME, temp_ext(TIME));
  Tview.set_min_max_range(0,20);
  Tview.set_title(title);
  Tview.fix_scale_width(3);

  // Time stepping:
  int nsteps = (int)(FINAL_TIME/TAU + 0.5);
  bool rhsonly = false;
  for(int ts = 1; ts <= nsteps; ts++)
  {
    info("---- Time step %d, time %3.5f, ext_temp %g", ts, TIME, temp_ext(TIME));

    // Assemble stiffness matrix and rhs.
    lp.assemble(mat, rhs, rhsonly);
    rhsonly = true;

    // Solve the matrix problem.
    if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

    // Update tsln.
    tsln.set_fe_solution(&space, rhs);

    // Update the time variable.
    TIME += TAU;

    // Visualize the solution.
    sprintf(title, "Time %3.2f, exterior temperature %3.5f", TIME, temp_ext(TIME));
    Tview.set_title(title);
    Tview.show(&tsln);
  }

  // Wait for the view to be closed.
  View::wait();
  return 0;
}
