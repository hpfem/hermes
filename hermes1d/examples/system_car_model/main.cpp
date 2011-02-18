#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"

//  This example solves a simplified car control model consisting of
//  five first-order equations equipped with the constraints
//  |alpha| <= alpha_max
//  |zeta| <= zeta_max
//  |v| <= v_max
//  |phi| <= phi_max
//  Goal: Calculate all possible trajectories of the car given the 
//  initial condition and intervals for alpha and zeta. The current 
//  algorithm considers 9 control parameters and moves only on the 
//  boundary of the 8-dimensional rectangle. TODO: forget trajectories 
//  where |v| > v_max or |phi| > phi_max.
//
//  PDE: x' - v cos(phi) cos(theta) = 0
//       y' - v cos(phi) sin(theta) = 0
//       v' - alpha                 = 0
//       phi' - zeta                = 0
//       theta' - v sin(phi)        = 0.
//
//  Time interval: (0, T).
//
//  Interval: (A, B).
//
//  BC: Homogenous Dirichlet.
//
//  The following parameters can be changed:
// Print data.
const int PRINT = 0;

//  The following parameters can be changed:
const int NEQ = 5;                      // Number of equations.
const int NELEM = 5;                    // Number of elements.
const double A = 0, B = 10;             // Domain end points.
const int P_INIT = 2;                   // Polynomial degree.

// Newton's method.
double NEWTON_TOL = 1e-5;               // Tolerance.
int NEWTON_MAX_ITER = 150;              // Max. number of Newton iterations.

// Parameters.
const double Alpha_max = 1.;
const double Zeta_max = M_PI/6;
const int Num_rays = 120;

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary conditions.
Hermes::vector<BCSpec *>DIR_BC_LEFT =  Hermes::vector<BCSpec *>(new BCSpec(0,0), new BCSpec(0,0), new BCSpec(0,0), new BCSpec(0,0), new BCSpec(0,0));
Hermes::vector<BCSpec *> DIR_BC_RIGHT = Hermes::vector<BCSpec *>();

// Controls.
const int N_ctrl = 4;
double time_ctrl[N_ctrl] = 
  {A, A + (A+B)/(N_ctrl-1.), A + 2.*(A+B)/(N_ctrl-1.), B};
double alpha_ctrl[N_ctrl] = {0, 0, 0, 0};
double zeta_ctrl[N_ctrl] = {0, 0, 0, 0};


// Include weak forms.
#include "forms.cpp"

void compute_trajectory(Space *space, DiscreteProblem *dp) 
{
  info("alpha = (%g, %g, %g, %g), zeta = (%g, %g, %g, %g)", 
         alpha_ctrl[0], alpha_ctrl[1], 
         alpha_ctrl[2], alpha_ctrl[3], zeta_ctrl[0], 
         zeta_ctrl[1], zeta_ctrl[2], zeta_ctrl[3]); 

  // Newton's loop.
  // Fill vector coeff_vec using dof and coeffs arrays in elements.
  double *coeff_vec = new double[Space::get_num_dofs(space)];
  get_coeff_vector(space, coeff_vec);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  int it = 1;
  while (1) 
  {
    // Obtain the number of degrees of freedom.
    int ndof = Space::get_num_dofs(space);

    // Assemble the Jacobian matrix and residual vector.
    dp->assemble(coeff_vec, matrix, rhs);

    // Calculate the l2-norm of residual vector.
    double res_l2_norm = get_l2_norm(rhs);

    // Info for user.
    info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(space), res_l2_norm);

    // If l2 norm of the residual vector is within tolerance, then quit.
    // NOTE: at least one full iteration forced
    //       here because sometimes the initial
    //       residual on fine mesh is too small.
    if(res_l2_norm < NEWTON_TOL && it > 1) break;

    // Multiply the residual vector with -1 since the matrix 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    for(int i=0; i<ndof; i++) rhs->set(i, -rhs->get(i));

    // Solve the linear system.
    if(!solver->solve())
      error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];

    // If the maximum number of iteration has been reached, then quit.
    if (it >= NEWTON_MAX_ITER) error ("Newton method did not converge.");
    
    // Copy coefficients from vector y to elements.
    set_coeff_vector(coeff_vec, space);

    it++;
  }
  
  // Cleanup.
  delete matrix;
  delete rhs;
  delete solver;
  delete [] coeff_vec;
}

void plot_trajectory(Space *space, int subdivision) 
{
  static int first_traj = 1;
  const char *traj_filename = "trajectory.gp";
  FILE *f;

  if (first_traj) 
  {
    f = fopen(traj_filename, "wb");
    first_traj = 0;
  }
  else 
  f = fopen(traj_filename, "ab");

  if(f == NULL) error("Problem opening trajectory file.");

  Linearizer l_traj(space);
  // Take first solution component as the x-coordinate.
  int x_comp = 0;

  // Take second solution component as the y-coordinate.
  int y_comp = 1;

  l_traj.plot_trajectory(f, x_comp, y_comp, subdivision);
  fclose(f);
}

void plot_trajectory_endpoint(Space *space) 
{
  static int first_traj = 1;
  const char *traj_filename = "reach.gp";
  FILE *f;

  if (first_traj) 
  {
    f = fopen(traj_filename, "wb");
    first_traj = 0;
  }
  else 
  f = fopen(traj_filename, "ab");

  if(f == NULL) error("Problem opening reachability file.");
  Linearizer l_reach(space);
  double x_ref = 1; 
  double x_phys; 
  double val[MAX_EQN_NUM];
  l_reach. eval_approx(space->last_active_element(), x_ref, &x_phys, val);

  // Take first solution component as the x-coordinate.
  int x_comp = 0;

  // Take second solution component as the y-coordinate.
  int y_comp = 1;

  fprintf(f, "%g %g\n", val[x_comp], val[y_comp]);
  fclose(f);
}

void plot_solution(Space *space, int subdivision)
{
  Linearizer l(space);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, subdivision);
}

// Set global controls.
void set_alpha_and_zeta(int component, double ray_angle, double radius) { 
  double alpha = radius * cos(ray_angle);
  double zeta = radius * sin(ray_angle); 
  if (alpha > Alpha_max) 
  {
    double coeff = Alpha_max/alpha;
    alpha = Alpha_max;
    zeta *= coeff;
  }
  if (alpha < -Alpha_max) 
  {
    double coeff = -Alpha_max/alpha;
    alpha = -Alpha_max;
    zeta *= coeff;
  }
  if (zeta > Zeta_max) 
  {
    double coeff = Zeta_max/zeta;
    zeta = Zeta_max;
    alpha *= coeff;
  }
  if (zeta < -Zeta_max) 
  {
    double coeff = -Zeta_max/zeta;
    zeta = -Zeta_max;
    alpha *= coeff;
  }
  alpha_ctrl[component] = alpha;
  zeta_ctrl[component] = zeta;
}


int main() 
{
  // Create space, set Dirichlet BC, enumerate basis functions.
  Space* space = new Space(A, B, NELEM, DIR_BC_LEFT, DIR_BC_RIGHT, P_INIT, NEQ);
  int ndof = Space::get_num_dofs(space);
  info("ndof: %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf(5);
  wf.add_matrix_form(0, 0, jacobian_0_0);
  wf.add_matrix_form(0, 2, jacobian_0_2);
  wf.add_matrix_form(0, 3, jacobian_0_3);
  wf.add_matrix_form(0, 4, jacobian_0_4);
  wf.add_matrix_form(1, 1, jacobian_1_1);
  wf.add_matrix_form(1, 2, jacobian_1_2);
  wf.add_matrix_form(1, 3, jacobian_1_3);
  wf.add_matrix_form(1, 4, jacobian_1_4);
  wf.add_matrix_form(2, 2, jacobian_2_2);
  wf.add_matrix_form(3, 3, jacobian_3_3);
  wf.add_matrix_form(4, 2, jacobian_4_2);
  wf.add_matrix_form(4, 3, jacobian_4_3);
  wf.add_matrix_form(4, 4, jacobian_4_4);
  wf.add_vector_form(0, residual_0);
  wf.add_vector_form(1, residual_1);
  wf.add_vector_form(2, residual_2);
  wf.add_vector_form(3, residual_3);
  wf.add_vector_form(4, residual_4);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem *dp = new DiscreteProblem(&wf, space, is_linear);
  
  // Move on the boundary of the rectangle 
  // (-Alpha_max, Alpha_max) x (-Zeta_max, Zeta_max) in the CCW
  // direction, starting at the point [Alpha_max, 0].
  double radius = sqrt(Alpha_max*Alpha_max + Zeta_max*Zeta_max);
  double angle_increment = 2.*M_PI/Num_rays;
  for (int ray_0 = 0; ray_0 < Num_rays; ray_0++) 
  {
    // Set alpha_ctrl[0], zeta_ctrl[0].
    set_alpha_and_zeta(0, ray_0*angle_increment, radius); 
    for (int ray_1 = 0; ray_1 < Num_rays; ray_1++) 
    {
      // Set alpha_ctrl[1], zeta_ctrl[1]. 
      set_alpha_and_zeta(1, ray_1*angle_increment, radius); 
      for (int ray_2 = 0; ray_2 < Num_rays; ray_2++) 
      {
        // Set alpha_ctrl[2], zeta_ctrl[2]. 
        set_alpha_and_zeta(2, ray_2*angle_increment, radius); 
        for (int ray_3 = 0; ray_3 < Num_rays; ray_3++) 
        {
          // Set alpha_ctrl[3], zeta_ctrl[3]. 
          set_alpha_and_zeta(3, ray_3*angle_increment, radius); 
  
          // Compute new trajectory for the actual set of global control 
          // parameters alpha_ctrl[] and zeta_ctrl[] via the Newton's 
          // method, using the last computed trajectory as the initial 
          // condition.
          compute_trajectory(space, dp);

          // Save trajectory to a file.
          int plotting_subdivision = 10;
          plot_trajectory_endpoint(space); 
        }
      }
    }
  }

  info("Done.");
  return 0;
}
