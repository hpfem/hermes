#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function.h"

// This test makes sure that example "richards" works correctly.

const int INIT_GLOB_REF_NUM = 3;       // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;        // Number of initial refinements towards boundary.
const int P_INIT = 3;                  // Initial polynomial degree.
const double TAU = 5e-3;               // Time step.
const double T_FINAL = 2*TAU + 1e-4;   // Time interval length.
const int TIME_INTEGRATION = 1;        // 1... implicit Euler, 2... Crank-Nicolson.
const double NEWTON_TOL = 1e-6;        // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;       // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// For the definition of initial condition.
int Y_POWER = 10;

// Problem parameters.
double K_S = 20.464;
double ALPHA = 0.001;
double THETA_R = 0;
double THETA_S = 0.45;
double H_R = -1000;
double A = 100;
double L = 100;
double STORATIVITY = 0.0;
double M = 0.6;
double N = 2.5;

#ifdef CONSTITUTIVE_GENUCHTEN
#include "constitutive_genuchten.cpp"
#else
#include "constitutive_gardner.cpp"
#endif

// Initial condition. It will be projected on the FE mesh 
// to obtain initial coefficient vector for the Newton's method.
double init_cond(double x, double y, double& dx, double& dy) {
  dx = (100 - 2*x)/2.5 * pow(y/100, Y_POWER);
  dy = x*(100 - x)/2.5 * pow(y/100, Y_POWER - 1) * 1./100;
  return x*(100 - x)/2.5 * pow(y/100, Y_POWER) - 1000;
}

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition markers.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  double dx, dy;
  return init_cond(x, y, dx, dy);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Solution for the previous time level.
  Solution u_prev_time;

  // Initialize the weak formulation.
  WeakForm wf;
  if (TIME_INTEGRATION == 1) {
    wf.add_matrix_form(jac_euler, jac_ord, HERMES_UNSYM, HERMES_ANY, &u_prev_time);
    wf.add_vector_form(res_euler, res_ord, HERMES_ANY, &u_prev_time);
  }
  else {
    wf.add_matrix_form(jac_cranic, jac_ord, HERMES_UNSYM, HERMES_ANY, &u_prev_time);
    wf.add_vector_form(res_cranic, res_ord, HERMES_ANY, &u_prev_time);
  }

  // Project the function init_cond() on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[Space::get_num_dofs(&space)];
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  OGProjection::project_global(&space, init_cond, coeff_vec, matrix_solver);
  Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);
  
  // Time stepping loop:
  double current_time = 0.0;
  int t_step = 1;
  do {
    info("---- Time step %d, t = %g s.", t_step, current_time); t_step++;

    // Initialize the FE problem.
    bool is_linear = false;
    DiscreteProblem dp(&wf, &space, is_linear);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Perform Newton's iteration.
    int it = 1;
    while (1)
    {
      // Obtain the number of degrees of freedom.
      int ndof = Space::get_num_dofs(&space);

      // Assemble the Jacobian matrix and residual vector.
      dp.assemble(coeff_vec, matrix, rhs, false);

      // Multiply the residual vector with -1 since the matrix 
      // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
      for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));
      
      // Calculate the l2-norm of residual vector.
      double res_l2_norm = get_l2_norm(rhs);

      // Info for user.
      info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(&space), res_l2_norm);

      // If l2 norm of the residual vector is within tolerance, or the maximum number 
      // of iteration has been reached, then quit.
      if (res_l2_norm < NEWTON_TOL || it > NEWTON_MAX_ITER) break;

      // Solve the linear system.
      if(!solver->solve())
        error ("Matrix solver failed.\n");

      // Add \deltaY^{n+1} to Y^n.
      for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];
      
      if (it >= NEWTON_MAX_ITER)
        error ("Newton method did not converge.");

      it++;
    }

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solution(coeff_vec, &space, &u_prev_time);

    // Cleanup.
    delete matrix;
    delete rhs;
    delete solver;

    // Update time.
    current_time += TAU;

  } while (current_time < T_FINAL);
  
  delete [] coeff_vec;

  info("Coordinate (  0,   0) value = %lf", u_prev_time.get_pt_value(0.0, 0.0));
  info("Coordinate ( 25,  25) value = %lf", u_prev_time.get_pt_value(25.0, 25.0));
  info("Coordinate ( 50,  50) value = %lf", u_prev_time.get_pt_value(50.0, 50.0));
  info("Coordinate ( 75,  75) value = %lf", u_prev_time.get_pt_value(75.0, 75.0));
  info("Coordinate (100, 100) value = %lf", u_prev_time.get_pt_value(100.0, 100.0));

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  double coor_x_y[5] = {0.0, 25.0, 50.0, 75.0, 100.0};
  double value[5] = {-1000.000000, -969.316013, -836.504249, -651.433710, -1000.000000};
  for (int i = 0; i < 5; i++)
  {
    if ((value[i] - u_prev_time.get_pt_value(coor_x_y[i], coor_x_y[i])) < 1E-6)
    {
    }
    else
    {
      printf("Failure!\n");
      return ERROR_FAILURE;
    }
  }
  printf("Success!\n");
  return ERROR_SUCCESS;
}

