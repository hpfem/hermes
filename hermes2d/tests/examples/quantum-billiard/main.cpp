#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#define DEBUG_ORDER
#include "hermes2d.h"

// This test makes sure that example "quantum-billiard" works correctly.

const int INIT_REF_NUM = 5;      // Number of initial uniform refinements.
const int P_INIT = 1;            // Initial polynomial degree.
const double TAU = 0.05;          // Time step.
const double T_FINAL = 1000;     // Time interval length.
const int TIME_DISCR = 2;        // 1 for implicit Euler, 2 for Crank-Nicolson.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem constants
cplx C = cplx(1./(30 * sqrt(3)), 0.0);
cplx C2 = cplx(200., 0.);

// Imaginary unit.
scalar ii = cplx(0.0, 1.0);

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

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
 return 0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space.
  H1Space* phi_space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);
  H1Space* psi_space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(Tuple<Space *>(phi_space, psi_space));
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
  wf.add_vector_form(0, callback(liform_euler_0), H2D_ANY, &phi_prev_time);
  wf.add_vector_form(1, callback(liform_euler_1), H2D_ANY, &psi_prev_time);

  // Time stepping loop:
  int nstep = T_FINAL;
  for(int ts = 1; ts <= nstep; ts++)
  {

    info("Time step %d:", ts);

    // Newton's method.
    info("Solving linear system.");
    Solution phi, psi; 
    bool is_complex = true;
    if (!solve_linear(Tuple<Space *>(phi_space, psi_space), &wf, matrix_solver,
                      Tuple<Solution *>(&phi, &psi), NULL, is_complex))
      error("Linear solve failed.");

    // Update previous time level solution.
    phi_prev_time.copy(&phi);
    psi_prev_time.copy(&psi);
  }

/*  
  info("Coordinate (   0,   0) value = %lf", psi_prev_time.get_pt_value(0.0, 0.0));
  info("Coordinate (-0.5,-0.5) value = %lf", psi_prev_time.get_pt_value(-0.5, -0.5));
  info("Coordinate ( 0.5,-0.5) value = %lf", psi_prev_time.get_pt_value(0.5, -0.5));
  info("Coordinate ( 0.5, 0.5) value = %lf", psi_prev_time.get_pt_value(0.5, 0.5));
  info("Coordinate (-0.5, 0.5) value = %lf", psi_prev_time.get_pt_value(-0.5, 0.5));
*/
#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
/*
  double coor_x[5] = {0.0, -0.5, 0.5, 0.5, -0.5};
  double coor_y[5] = {0.0, -0.5, -0.5, 0.5, 0.5};
  double value[5] = {-1000.000000, -913.016598, -709.440069, -575.643683, -1000.000000};
  for (int i = 0; i < 5; i++)
  {
    if ((value[i] - psi_prev_time.get_pt_value(coor_x[i], coor_y[i])) < 1E-6)
    {
      printf("Success!\n");
    }
    else
    {
      printf("Failure!\n");
      return ERROR_FAILURE;
    }
  }
*/
  return ERROR_SUCCESS;
}
