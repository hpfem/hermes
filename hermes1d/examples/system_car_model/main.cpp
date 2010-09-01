#include "hermes1d.h"

// ********************************************************************
// This example solves a simplified car control model consisting of
// five first-order equations
//
// x' - v cos(phi) cos(theta) = 0
// y' - v cos(phi) sin(theta) = 0
// v' - alpha                 = 0
// phi' - zeta                = 0
// theta' - v sin(phi)        = 0
//
// equipped with the constraints
//
// |alpha| <= alpha_max
// |zeta| <= zeta_max
// |v| <= v_max
// |phi| <= phi_max

// The problem is considered in a time interval (0, T), and Dirichlet 
// conditions are given for all quantities at the beginning.

// Goal: Calculate all possible trajectories of the car given the 
// initial condition and intervals for alpha and zeta. The current 
// algorithm considers 9 control parameters and moves only on the 
// boundary of the 8-dimensional rectangle. TODO: forget trajectories 
// where |v| > v_max or |phi| > phi_max.

// Print data ?
const int PRINT = 0;

// General input:
const int N_eq = 5;
const int N_elem = 5;              // number of elements
const double A = 0, B = 10;         // domain end points
const int P_init = 2;              // initial polynomal degree 

// Parameters
const double Alpha_max = 1.;
const double Zeta_max = M_PI/6;
const int Num_rays = 120;

// Boundary conditions
const double X0_left = 0;
const double Y0_left = 0;
const double Vel_left = 0;
const double Phi_left = 0;
const double Theta_left = 0;

// Controls
const int N_ctrl = 4;
double time_ctrl[N_ctrl] = 
  {A, A + (A+B)/(N_ctrl-1.), A + 2.*(A+B)/(N_ctrl-1.), B};
double alpha_ctrl[N_ctrl] = {0, 0, 0, 0};
double zeta_ctrl[N_ctrl] = {0, 0, 0, 0};

// Newton's method
double NEWTON_TOL = 1e-5;
int NEWTON_MAXITER = 150;

// ********************************************************************

// Include weak forms
#include "forms.cpp"

void compute_trajectory(Mesh *mesh, DiscreteProblem *dp) 
{
  printf("alpha = (%g, %g, %g, %g), zeta = (%g, %g, %g, %g)\n", 
         alpha_ctrl[0], alpha_ctrl[1], 
         alpha_ctrl[2], alpha_ctrl[3], zeta_ctrl[0], 
         zeta_ctrl[1], zeta_ctrl[2], zeta_ctrl[3]); 

  // Newton's loop
  newton(dp, mesh, NULL, NEWTON_TOL, NEWTON_MAXITER);
}

void plot_trajectory(Mesh *mesh, int subdivision) 
{
  static int first_traj = 1;
  const char *traj_filename = "trajectory.gp";
  FILE *f;
  if (first_traj) {
    f = fopen(traj_filename, "wb");
    first_traj = 0;
  }
  else f = fopen(traj_filename, "ab");
  if(f == NULL) error("Problem opening trajectory file.");
  Linearizer l_traj(mesh);
  int x_comp = 0; // take first solution component as the x-coordinate
  int y_comp = 1; // take second solution component as the y-coordinate
  l_traj.plot_trajectory(f, x_comp, y_comp, subdivision);
  fclose(f);
}

void plot_trajectory_endpoint(Mesh *mesh) 
{
  static int first_traj = 1;
  const char *traj_filename = "reach.gp";
  FILE *f;
  if (first_traj) {
    f = fopen(traj_filename, "wb");
    first_traj = 0;
  }
  else f = fopen(traj_filename, "ab");
  if(f == NULL) error("Problem opening reachability file.");
  Linearizer l_reach(mesh);
  double x_ref = 1; 
  double x_phys; 
  double val[MAX_EQN_NUM];
  l_reach. eval_approx(mesh->last_active_element(), x_ref, &x_phys, val);
  int x_comp = 0; // take first solution component as the x-coordinate
  int y_comp = 1; // take second solution component as the y-coordinate
  fprintf(f, "%g %g\n", val[x_comp], val[y_comp]);
  fclose(f);
}

void plot_solution(Mesh *mesh, int subdivision)
{
  Linearizer l(mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename, subdivision);
}

// set global controls
void set_alpha_and_zeta(int component, double ray_angle, double radius) { 
  double alpha = radius * cos(ray_angle);
  double zeta = radius * sin(ray_angle); 
  if (alpha > Alpha_max) {
    double coeff = Alpha_max/alpha;
    alpha = Alpha_max;
    zeta *= coeff;
  }
  if (alpha < -Alpha_max) {
    double coeff = -Alpha_max/alpha;
    alpha = -Alpha_max;
    zeta *= coeff;
  }
  if (zeta > Zeta_max) {
    double coeff = Zeta_max/zeta;
    zeta = Zeta_max;
    alpha *= coeff;
  }
  if (zeta < -Zeta_max) {
    double coeff = -Zeta_max/zeta;
    zeta = -Zeta_max;
    alpha *= coeff;
  }
  alpha_ctrl[component] = alpha;
  zeta_ctrl[component] = zeta;
}

/******************************************************************************/
int main() {
  // create mesh
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  mesh->set_bc_left_dirichlet(0, X0_left);
  mesh->set_bc_left_dirichlet(1, Y0_left);
  mesh->set_bc_left_dirichlet(2, Vel_left);
  mesh->set_bc_left_dirichlet(3, Phi_left);
  mesh->set_bc_left_dirichlet(4, Theta_left);

  // enumerate shape functions and calculate
  // the number of DOF
  printf("N_dof = %d\n", mesh->assign_dofs());

  // register weak forms
  DiscreteProblem *dp = new DiscreteProblem();
  dp->add_matrix_form(0, 0, jacobian_0_0);
  dp->add_matrix_form(0, 2, jacobian_0_2);
  dp->add_matrix_form(0, 3, jacobian_0_3);
  dp->add_matrix_form(0, 4, jacobian_0_4);
  dp->add_matrix_form(1, 1, jacobian_1_1);
  dp->add_matrix_form(1, 2, jacobian_1_2);
  dp->add_matrix_form(1, 3, jacobian_1_3);
  dp->add_matrix_form(1, 4, jacobian_1_4);
  dp->add_matrix_form(2, 2, jacobian_2_2);
  dp->add_matrix_form(3, 3, jacobian_3_3);
  dp->add_matrix_form(4, 2, jacobian_4_2);
  dp->add_matrix_form(4, 3, jacobian_4_3);
  dp->add_matrix_form(4, 4, jacobian_4_4);
  dp->add_vector_form(0, residual_0);
  dp->add_vector_form(1, residual_1);
  dp->add_vector_form(2, residual_2);
  dp->add_vector_form(3, residual_3);
  dp->add_vector_form(4, residual_4);

  // Move on the boundary of the rectangle 
  // (-Alpha_max, Alpha_max) x (-Zeta_max, Zeta_max) in the CCW
  // direction, starting at the point [Alpha_max, 0]
  double radius = sqrt(Alpha_max*Alpha_max + Zeta_max*Zeta_max);
  double angle_increment = 2.*M_PI/Num_rays;
  for (int ray_0 = 0; ray_0 < Num_rays; ray_0++) {
    // set alpha_ctrl[0], zeta_ctrl[0].
    set_alpha_and_zeta(0, ray_0*angle_increment, radius); 
    for (int ray_1 = 0; ray_1 < Num_rays; ray_1++) {
      // set alpha_ctrl[1], zeta_ctrl[1]. 
      set_alpha_and_zeta(1, ray_1*angle_increment, radius); 
      for (int ray_2 = 0; ray_2 < Num_rays; ray_2++) {
        // set alpha_ctrl[2], zeta_ctrl[2]. 
        set_alpha_and_zeta(2, ray_2*angle_increment, radius); 
        for (int ray_3 = 0; ray_3 < Num_rays; ray_3++) {
          // set alpha_ctrl[3], zeta_ctrl[3]. 
          set_alpha_and_zeta(3, ray_3*angle_increment, radius); 
  
          // Compute new trajectory for the actual set of global control 
          // parameters alpha_ctrl[] and zeta_ctrl[] via the Newton's 
          // method, using the last computed trajectory as the initial 
          // condition  
          compute_trajectory(mesh, dp);

          // save trajectory to a file
          int plotting_subdivision = 10;
          plot_trajectory_endpoint(mesh); 

          // save trajectory to a file
          //int plotting_subdivision = 10;
          //plot_trajectory(mesh, plotting_subdivision); 

          // save solution to a file
          //int plotting_subdivision_2 = 10;
          //plot_solution(mesh, plotting_subdivision_2); 
        }
      }
    }
  }

  printf("Done.\n");
  return 1;
}
