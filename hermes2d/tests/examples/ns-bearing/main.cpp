#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "forms.h"

// This test makes sure that example "ns-bearing" works correctly.

const int INIT_REF_NUM = 2;                // Number of initial uniform mesh refinements. 
const int INIT_BDY_REF_NUM_INNER = 2;      // Number of initial mesh refinements towards boundary. 
const int INIT_BDY_REF_NUM_OUTER = 2;      // Number of initial mesh refinements towards boundary. 

//#define STOKES                     // If this is defined, Stokes problem is solved, otherwise N-S.
#define PRESSURE_IN_L2               // If this is defined, the pressure is approximated using
                                     // discontinuous L2 elements (making the velocity discreetely
                                     // divergence-free, more accurate than using a continuous
                                     // pressure approximation). Otherwise the standard continuous
                                     // elements are used. The results are striking - check the
                                     // tutorial for comparisons.
const bool NEWTON = true;            // If NEWTON == true then the Newton's iteration is performed.
                                     // in every time step. Otherwise the convective term is linearized
                                     // using the velocities from the previous time step.
const int P_INIT_VEL = 2;            // Initial polynomial degree for velocity components.
const int P_INIT_PRESSURE = 1;       // Initial polynomial degree for pressure.
                                     // Note: P_INIT_VEL should always be greater than
                                     // P_INIT_PRESSURE because of the inf-sup condition.
const double RE = 5000.0;            // Reynolds number.
const double VEL = 0.1;              // Surface velocity of inner circle.
const double STARTUP_TIME = 1.0;     // During this time, surface velocity of the inner circle increases 
                                     // gradually from 0 to VEL, then it stays constant.
const double TAU = 10.0;             // Time step.
const double T_FINAL = 2*TAU + 1e-4;       // Time interval length.
const double NEWTON_TOL = 1e-5;      // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 10;      // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Boundary markers.
int bdy_inner = 1;
int bdy_outer = 2;

// Current time (used in weak forms).
double TIME = 0;

// Boundary condition types for x-velocity.
BCType xvel_bc_type(int marker) {
  return BC_ESSENTIAL;
}

// Boundary condition types for y-velocity.
BCType yvel_bc_type(int marker) {
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values for x-velocity.
scalar essential_bc_values_xvel(int ess_bdy_marker, double x, double y) {
  if (ess_bdy_marker == bdy_inner) {
    // Time-dependent surface velocity of inner circle.
    double velocity;
    if (TIME <= STARTUP_TIME) velocity = VEL * TIME/STARTUP_TIME;
    else velocity = VEL;
    double alpha = atan2(x, y);
    double xvel = velocity*cos(alpha);
    //printf("%g %g xvel = %g\n", x, y, xvel);
    return xvel; 
  }
  else return 0;
}

// Essential (Dirichlet) boundary condition values for y-velocity.
scalar essential_bc_values_yvel(int ess_bdy_marker, double x, double y) {
  if (ess_bdy_marker == bdy_inner) {
    // Time-dependent surface velocity of inner circle.
    double velocity;
    if (TIME <= STARTUP_TIME) velocity = VEL * TIME/STARTUP_TIME;
    else velocity = VEL;
    double alpha = atan2(x, y);
    double yvel = -velocity*sin(alpha);
    //printf("%g %g yvel = %g\n", x, y, yvel);
    return yvel; 
  }
  else return 0;
}

// Weak forms.
#include "forms.cpp"

// Custom function to calculate drag coefficient.
double integrate_over_wall(MeshFunction* meshfn, int marker)
{
  Quad2D* quad = &g_quad_2d_std;
  meshfn->set_quad_2d(quad);

  double integral = 0.0;
  Element* e;
  Mesh* mesh = meshfn->get_mesh();

  for_all_active_elements(e, mesh)
  {
    for(int edge = 0; edge < e->nvert; edge++)
    {
      if ((e->en[edge]->bnd) && (e->en[edge]->marker == marker))
      {
        update_limit_table(e->get_mode());
        RefMap* ru = meshfn->get_refmap();

        meshfn->set_active_element(e);
        int eo = quad->get_edge_points(edge);
        meshfn->set_quad_order(eo, H2D_FN_VAL);
        scalar *uval = meshfn->get_fn_values();
        double3* pt = quad->get_points(eo);
        double3* tan = ru->get_tangent(edge);
        for (int i = 0; i < quad->get_num_points(eo); i++)
          integral += pt[i][2] * uval[i] * tan[i][2];
      }
    }
  }
  return integral * 0.5;
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain-excentric.mesh", &mesh);
  //mloader.load("domain-concentric.mesh", &mesh);

  // Initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(bdy_inner, INIT_BDY_REF_NUM_INNER, false);  // true for anisotropic refinements
  mesh.refine_towards_boundary(bdy_outer, INIT_BDY_REF_NUM_OUTER, false);  // false for isotropic refinements

  // Create spaces with default shapesets. 
  H1Space* xvel_space = new H1Space(&mesh, xvel_bc_type, essential_bc_values_xvel, P_INIT_VEL);
  H1Space* yvel_space = new H1Space(&mesh, yvel_bc_type, essential_bc_values_yvel, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space* p_space = new L2Space(&mesh, P_INIT_PRESSURE);
#else
  H1Space* p_space = new H1Space(&mesh, NULL, NULL, P_INIT_PRESSURE);
#endif

  // Calculate and report the number of degrees of freedom.
  int ndof = get_num_dofs(Tuple<Space *>(xvel_space, yvel_space, p_space));
  info("ndof = %d.", ndof);

  // Define projection norms.
  int vel_proj_norm = H2D_H1_NORM;
#ifdef PRESSURE_IN_L2
  int p_proj_norm = H2D_L2_NORM;
#else
  int p_proj_norm = H2D_H1_NORM;
#endif

  // Solutions for the Newton's iteration and time stepping.
  info("Setting initial conditions.");
  Solution xvel_prev_time, yvel_prev_time, p_prev_time;
  xvel_prev_time.set_zero(&mesh);
  yvel_prev_time.set_zero(&mesh);
  p_prev_time.set_zero(&mesh);

  // Initialize weak formulation.
  WeakForm wf(3);
  if (NEWTON) {
    wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_matrix_form(0, 0, callback(newton_bilinear_form_unsym_0_0), H2D_UNSYM, H2D_ANY);
    wf.add_matrix_form(0, 1, callback(newton_bilinear_form_unsym_0_1), H2D_UNSYM, H2D_ANY);
    wf.add_matrix_form(0, 2, callback(bilinear_form_unsym_0_2), H2D_ANTISYM);
    wf.add_matrix_form(1, 0, callback(newton_bilinear_form_unsym_1_0), H2D_UNSYM, H2D_ANY);
    wf.add_matrix_form(1, 1, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_matrix_form(1, 1, callback(newton_bilinear_form_unsym_1_1), H2D_UNSYM, H2D_ANY);
    wf.add_matrix_form(1, 2, callback(bilinear_form_unsym_1_2), H2D_ANTISYM);
    wf.add_vector_form(0, callback(newton_F_0), H2D_ANY, 
                       Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_vector_form(1, callback(newton_F_1), H2D_ANY, 
                       Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_vector_form(2, callback(newton_F_2), H2D_ANY);
  }
  else {
    wf.add_matrix_form(0, 0, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_matrix_form(0, 0, callback(simple_bilinear_form_unsym_0_0_1_1), 
                  H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_matrix_form(1, 1, callback(bilinear_form_sym_0_0_1_1), H2D_SYM);
    wf.add_matrix_form(1, 1, callback(simple_bilinear_form_unsym_0_0_1_1), 
                  H2D_UNSYM, H2D_ANY, Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time));
    wf.add_matrix_form(0, 2, callback(bilinear_form_unsym_0_2), H2D_ANTISYM);
    wf.add_matrix_form(1, 2, callback(bilinear_form_unsym_1_2), H2D_ANTISYM);
    wf.add_vector_form(0, callback(simple_linear_form), H2D_ANY, &xvel_prev_time);
    wf.add_vector_form(1, callback(simple_linear_form), H2D_ANY, &yvel_prev_time);
  }

  // Project initial conditions on FE spaces to obtain initial coefficient 
  // vector for the Newton's method.
  Vector* coeff_vec;
  if (NEWTON) {
    info("Projecting initial conditions to obtain initial vector for the Newton's method.");
    coeff_vec = new AVector(ndof);
    project_global(Tuple<Space *>(xvel_space, yvel_space, p_space),
                   Tuple<int>(vel_proj_norm, vel_proj_norm, p_proj_norm),
                   Tuple<MeshFunction*>(&xvel_prev_time, &yvel_prev_time, &p_prev_time),
                   Tuple<Solution*>(&xvel_prev_time, &yvel_prev_time, &p_prev_time),
                   coeff_vec);
  }
  else coeff_vec = NULL;

  // Time-stepping loop:
  char title[100];
  int num_time_steps = T_FINAL / TAU;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    TIME += TAU;
    info("---- Time step %d, time = %g:", ts, TIME);

    // Update time-dependent essential BC are used.
    if (TIME <= STARTUP_TIME) {
      info("Updating time-dependent essential BC.");
      update_essential_bc_values(Tuple<Space *>(xvel_space, yvel_space, p_space));
    }

    if (NEWTON) {
      // Newton's method.
      info("Performing Newton's method.");
      bool verbose = true; // Default is false.
      if (!solve_newton(Tuple<Space *>(xvel_space, yvel_space, p_space), &wf, coeff_vec,
          matrix_solver, NEWTON_TOL, NEWTON_MAX_ITER, verbose))
        error("Newton's method did not converge.");

      // Update previous time level solutions.
      xvel_prev_time.set_fe_solution(xvel_space, coeff_vec);
      yvel_prev_time.set_fe_solution(yvel_space, coeff_vec);
      p_prev_time.set_fe_solution(p_space, coeff_vec);
    }
    else {
      // Linear solve.  
      info("Assembling and solving linear problem.");
      solve_linear(Tuple<Space *>(xvel_space, yvel_space, p_space), &wf, matrix_solver,
                   Tuple<Solution*>(&xvel_prev_time, &yvel_prev_time, &p_prev_time));
    }
 }
  info("Coordinate ( 0.1, 0.0) xvel value = %lf", xvel_prev_time.get_pt_value(0.1, 0.0));
  info("Coordinate ( 0.5, 0.0) xvel value = %lf", xvel_prev_time.get_pt_value(0.5, 0.0));
  info("Coordinate ( 0.9, 0.0) xvel value = %lf", xvel_prev_time.get_pt_value(0.9, 0.0));
  info("Coordinate ( 1.3, 0.0) xvel value = %lf", xvel_prev_time.get_pt_value(1.3, 0.0));
  info("Coordinate ( 1.7, 0.0) xvel value = %lf", xvel_prev_time.get_pt_value(1.7, 0.0));

  info("Coordinate ( 0.1, 0.0) yvel value = %lf", yvel_prev_time.get_pt_value(0.1, 0.0));
  info("Coordinate ( 0.5, 0.0) yvel value = %lf", yvel_prev_time.get_pt_value(0.5, 0.0));
  info("Coordinate ( 0.9, 0.0) yvel value = %lf", yvel_prev_time.get_pt_value(0.9, 0.0));
  info("Coordinate ( 1.3, 0.0) yvel value = %lf", yvel_prev_time.get_pt_value(1.3, 0.0));
  info("Coordinate ( 1.7, 0.0) yvel value = %lf", yvel_prev_time.get_pt_value(1.7, 0.0));

#define ERROR_SUCCESS                                0
#define ERROR_FAILURE                               -1
  int success = 1;
  double eps = 1e-5;
  if (fabs(xvel_prev_time.get_pt_value(0.1, 0.0) - (0.000000)) > eps) {
    printf("Coordinate ( 0.1, 0.0) xvel value = %lf\n", xvel_prev_time.get_pt_value(0.1, 0.0));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(0.5, 0.0) - (-0.000403)) > eps) {
    printf("Coordinate ( 0.5, 0.0) xvel value = %lf\n", xvel_prev_time.get_pt_value(0.5, 0.0));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(0.9, 0.0) - (-0.000057)) > eps) {
    printf("Coordinate ( 0.9, 0.0) xvel value = %lf\n", xvel_prev_time.get_pt_value(0.9, 0.0));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(1.3, 0.0) - (-0.000006)) > eps) {
    printf("Coordinate ( 1.3, 0.0) xvel value = %lf\n", xvel_prev_time.get_pt_value(1.3, 0.0));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(1.7, 0.0) - (-0.000001)) > eps) {
    printf("Coordinate ( 1.7, 0.0) xvel value = %lf\n", xvel_prev_time.get_pt_value(1.7, 0.0));
    success = 0;
  }

  if (fabs(yvel_prev_time.get_pt_value(0.1, 0.0) - (-0.100000)) > eps) {
    printf("Coordinate ( 0.1, 0.0) yvel value = %lf\n", yvel_prev_time.get_pt_value(0.1, 0.0));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(0.5, 0.0) - (0.001486)) > eps) {
    printf("Coordinate ( 0.5, 0.0) yvel value = %lf\n", yvel_prev_time.get_pt_value(0.5, 0.0));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(0.9, 0.0) - (0.000577)) > eps) {
    printf("Coordinate ( 0.9, 0.0) yvel value = %lf\n", yvel_prev_time.get_pt_value(0.9, 0.0));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(1.3, 0.0) - (0.000299)) > eps) {
    printf("Coordinate ( 1.3, 0.0) yvel value = %lf\n", yvel_prev_time.get_pt_value(1.3, 0.0));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(1.7, 0.0) - (0.000170)) > eps) {
    printf("Coordinate ( 1.7, 0.0) yvel value = %lf\n", yvel_prev_time.get_pt_value(1.7, 0.0));
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
