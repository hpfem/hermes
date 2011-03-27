#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#define DEBUG_ORDER
#include "hermes2d.h"

//  This is a simple version of the quantum billiard problem, described with 
//  a complex-valued linear time-dependent wave equation. The equation is converted 
//  into a system of two PDE which are first-order in time. For time-discretization 
//  one can use either the first-order implicit Euler method or the second-order 
//  Crank-Nicolson method.
//
//  PDE: 
//
//  \frac{\partial^2}{\partial t^2}\Psi - \frac{\partial^2}{\partial x^2}\Psi
//  - \frac{\partial^2}{\partial y^2}\Psi = 0.
//  
//  Domain: square (-1, 1)^2.
//
//  BC: homogeneous Dirichlet everywhere on the boundary.
//
//  IC: Gaussian distribution \Psi(0, x, y) = 

const int INIT_REF_NUM = 5;                       // Number of initial uniform refinements.
const int P_INIT = 1;                             // Initial polynomial degree.
const double TAU = 0.05;                          // Time step.
const double T_FINAL = 1000;                      // Time interval length.
const int TIME_DISCR = 2;                         // 1 for implicit Euler, 2 for Crank-Nicolson.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem constants
cplx C = cplx(1./(30 * sqrt((double)3)), 0.0);
cplx C2 = cplx(200., 0.);

// Imaginary unit.
scalar ii = cplx(0.0, 1.0);

// Boundary markers.
const int BDY_BOTTOM = 1;
const int BDY_RIGHT = 2;
const int BDY_TOP = 3;
const int BDY_LEFT = 4;

// Initial condition for Psi.
scalar init_cond_psi(double x, double y, scalar& dx, scalar& dy)
{
  scalar val = exp(-(x*x + y*y)/(2.*C*C)) * exp(C2 * ii * x);
  dx = (-x/(C*C)+ii*C2)*val;
  dy = (-y/(C*C))*val;
  return val;
}

// Initial condition for Phi.
scalar init_cond_phi(double x, double y, scalar& dx, scalar& dy)
{
  scalar val = ii * C2 * exp(-(x*x + y*y)/(2.*C*C)) * exp(C2 * ii * x);
  dx = (-x/(C*C)+ii*C2)*val;
  dy = (-y/(C*C))*val;
  return val;
}

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_zero(Hermes::vector<int>(BDY_BOTTOM, BDY_RIGHT, BDY_TOP, BDY_LEFT));

  // Create an H1 space.
  H1Space* phi_space = new H1Space(&mesh, &bc_types, &bc_values, P_INIT);
  H1Space* psi_space = new H1Space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(Hermes::vector<Space *>(phi_space, psi_space));
  info("ndof = %d.", ndof);

  // Initialize previous time level solutions.
  Solution phi_prev_time, psi_prev_time;
  phi_prev_time.set_exact(&mesh, init_cond_phi);
  psi_prev_time.set_exact(&mesh, init_cond_psi);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(biform_euler_0_0));
  wf.add_matrix_form(0, 1, callback(biform_euler_0_1));
  wf.add_matrix_form(1, 0, callback(biform_euler_1_0));
  wf.add_matrix_form(1, 1, callback(biform_euler_1_1));
  wf.add_vector_form(0, callback(liform_euler_0), HERMES_ANY, &phi_prev_time);
  wf.add_vector_form(1, callback(liform_euler_1), HERMES_ANY, &psi_prev_time);

  // Initialize views.
  ScalarView view("Psi", new WinGeom(0, 0, 600, 500));
  view.fix_scale_width(80);

  // Time stepping loop:
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= nstep; ts++)
  {

    info("Time step %d:", ts);

    info("Solving linear system.");
    // Initialize the FE problem.
    bool is_linear = true;
    DiscreteProblem dp(&wf, Hermes::vector<Space *>(phi_space, psi_space), is_linear);
   
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Assemble the stiffness matrix and right-hand side vector.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the linear system and if successful, obtain the solution.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution::vector_to_solutions(solver->get_solution(), Hermes::vector<Space *>(phi_space, psi_space), Hermes::vector<Solution *>(&phi_prev_time, &psi_prev_time));
    else
      error ("Matrix solver failed.\n");

    // Show the new time level solution.
    char title[100];
    sprintf(title, "Time step %d", ts);
    view.set_title(title);
    view.show(&psi_prev_time);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
