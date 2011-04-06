#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This test makes sure that example "ns-bearing" works correctly.

const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements. 
const int INIT_BDY_REF_NUM_INNER = 2;             // Number of initial mesh refinements towards boundary. 
const int INIT_BDY_REF_NUM_OUTER = 2;             // Number of initial mesh refinements towards boundary. 

const bool STOKES = false;                        // For application of Stokes flow (creeping flow).
#define PRESSURE_IN_L2                            // If this is defined, the pressure is approximated using
                                                  // discontinuous L2 elements (making the velocity discreetely
                                                  // divergence-free, more accurate than using a continuous
                                                  // pressure approximation). Otherwise the standard continuous
                                                  // elements are used. The results are striking - check the
                                                  // tutorial for comparisons.
const bool NEWTON = true;                         // If NEWTON == true then the Newton's iteration is performed.
                                                  // in every time step. Otherwise the convective term is linearized
                                                  // using the velocities from the previous time step.
const int P_INIT_VEL = 2;                         // Initial polynomial degree for velocity components.
const int P_INIT_PRESSURE = 1;                    // Initial polynomial degree for pressure.
                                                  // Note: P_INIT_VEL should always be greater than
                                                  // P_INIT_PRESSURE because of the inf-sup condition.
const double RE = 5000.0;                         // Reynolds number.
const double VEL = 0.1;                           // Surface velocity of inner circle.
const double STARTUP_TIME = 1.0;                  // During this time, surface velocity of the inner circle increases 
                                                  // gradually from 0 to VEL, then it stays constant.
const double TAU = 10.0;                          // Time step.
const double T_FINAL = 2*TAU + 1e-4;              // Time interval length.
const double NEWTON_TOL = 1e-1;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 10;                   // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const std::string BDY_INNER = "Inner";
const std::string BDY_OUTER = "Outer";

// Current time (used in weak forms).
double current_time = 0;

// Weak forms.
#include "../definitions.cpp"

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
  Hermes2D hermes_2D;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain-excentric.mesh", &mesh);
  //mloader.load("../domain-concentric.mesh", &mesh);

  // Initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_INNER, INIT_BDY_REF_NUM_INNER, false);  // true for anisotropic refinements
  mesh.refine_towards_boundary(BDY_OUTER, INIT_BDY_REF_NUM_OUTER, false);  // false for isotropic refinements

  // Initialize boundary conditions.
  EssentialBCNonConstX bc_inner_vel_x(BDY_INNER, VEL, STARTUP_TIME);
  EssentialBCNonConstY bc_inner_vel_y(BDY_INNER, VEL, STARTUP_TIME);
  DefaultEssentialBCConst bc_outer_vel(BDY_OUTER, 0.0);
  EssentialBCs bcs_vel_x(Hermes::vector<EssentialBoundaryCondition *>(&bc_inner_vel_x, &bc_outer_vel));
  EssentialBCs bcs_vel_y(Hermes::vector<EssentialBoundaryCondition *>(&bc_inner_vel_y, &bc_outer_vel));
  EssentialBCs bcs_pressure;

  // Spaces for velocity components and pressure.
  H1Space xvel_space(&mesh, &bcs_vel_x, P_INIT_VEL);
  H1Space yvel_space(&mesh, &bcs_vel_y, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space p_space(&mesh, &bcs_pressure, P_INIT_PRESSURE);
#else
  H1Space p_space(&mesh, &bcs_pressure, P_INIT_PRESSURE);
#endif

  // Calculate and report the number of degrees of freedom.
  int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space));

  // Define projection norms.
  ProjNormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  ProjNormType p_proj_norm = HERMES_L2_NORM;
#else
  ProjNormType p_proj_norm = HERMES_H1_NORM;
#endif

  // Solutions for the Newton's iteration and time stepping.
  Solution xvel_prev_time, yvel_prev_time, p_prev_time;
  xvel_prev_time.set_zero(&mesh);
  yvel_prev_time.set_zero(&mesh);
  p_prev_time.set_zero(&mesh);

  // Initialize weak formulation.
  WeakForm* wf;
  if (NEWTON)
    wf = new WeakFormNSNewton(STOKES, RE, TAU, &xvel_prev_time, &yvel_prev_time);
  else
    wf = new WeakFormNSSimpleLinearization(STOKES, RE, TAU, &xvel_prev_time, &yvel_prev_time);

  // Initialize the FE problem.
  bool is_linear;
  if (NEWTON) is_linear = false;
  else is_linear = true;
  DiscreteProblem dp(wf, Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[Space::get_num_dofs(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space))];
  if (NEWTON) {
    OGProjection::project_global(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space),
                   Hermes::vector<MeshFunction*>(&xvel_prev_time, &yvel_prev_time, &p_prev_time),
                   coeff_vec, matrix_solver, 
                   Hermes::vector<ProjNormType>(vel_proj_norm, vel_proj_norm, p_proj_norm));
  }

  // Time-stepping loop:
  char title[100];
  int num_time_steps = T_FINAL / TAU;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    current_time += TAU;
    
    // Update time-dependent essential BCs.
    Space::update_essential_bc_values(Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), current_time);

    if (NEWTON) 
    {
      // Perform Newton's iteration.
      bool verbose = false;
      if (!hermes_2D.solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
          NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");
      
      // Update previous time level solutions.
      Solution::vector_to_solutions(coeff_vec, Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), 
                                    Hermes::vector<Solution *>(&xvel_prev_time, &yvel_prev_time, &p_prev_time));
    }
    else {
      // Linear solve.  
      dp.assemble(matrix, rhs, false);
      if(solver->solve())
        Solution::vector_to_solutions(solver->get_solution(), 
                  Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), 
                                      Hermes::vector<Solution *>(&xvel_prev_time, &yvel_prev_time, &p_prev_time));
      else
        error ("Matrix solver failed.\n");
    }
  }

  // Clean up.
  delete [] coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;
  
  info("Coordinate ( -0.7, 0.0) xvel value = %lf", xvel_prev_time.get_pt_value(-0.7, 0.1));
  info("Coordinate ( 0.65, 0.0) xvel value = %lf", xvel_prev_time.get_pt_value(0.65, 0.0));
  info("Coordinate ( 0.7, 0.0) xvel value = %lf", xvel_prev_time.get_pt_value(0.7, 0.0));
  info("Coordinate ( 0.9, 0.0) xvel value = %lf", xvel_prev_time.get_pt_value(0.9, 0.0));
  info("Coordinate ( 1.1, 0.0) xvel value = %lf", xvel_prev_time.get_pt_value(1.1, 0.0));

  info("Coordinate ( -0.7, 0.0) yvel value = %lf", yvel_prev_time.get_pt_value(-0.7, 0.1));
  info("Coordinate ( 0.65, 0.0) yvel value = %lf", yvel_prev_time.get_pt_value(0.65, 0.0));
  info("Coordinate ( 0.7, 0.0) yvel value = %lf", yvel_prev_time.get_pt_value(0.7, 0.0));
  info("Coordinate ( 0.9, 0.0) yvel value = %lf", yvel_prev_time.get_pt_value(0.9, 0.0));
  info("Coordinate ( 1.1, 0.0) yvel value = %lf", yvel_prev_time.get_pt_value(1.1, 0.0));

  int success = 1;
  double eps = 1e-5;
  if (fabs(xvel_prev_time.get_pt_value(-0.7, 0.1) - (0.003199)) > eps) {
    printf("Coordinate ( -0.7, 0.1) xvel value = %lf\n", xvel_prev_time.get_pt_value(-0.7, 0.1));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(0.65, 0.0) - (-0.001314)) > eps) {
    printf("Coordinate ( 0.65, 0.0) xvel value = %lf\n", xvel_prev_time.get_pt_value(0.65, 0.0));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(0.7, 0.0) - (-0.001888)) > eps) {
    printf("Coordinate ( 0.7, 0.0) xvel value = %lf\n", xvel_prev_time.get_pt_value(0.7, 0.0));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(0.9, 0.0) - (-0.001262)) > eps) {
    printf("Coordinate ( 0.9, 0.0) xvel value = %lf\n", xvel_prev_time.get_pt_value(0.9, 0.0));
    success = 0;
  }
  if (fabs(xvel_prev_time.get_pt_value(1.1, 0.0) - (-0.000218)) > eps) {
    printf("Coordinate ( 1.1, 0.0) xvel value = %lf\n", xvel_prev_time.get_pt_value(1.1, 0.0));
    success = 0;
  }

  if (fabs(yvel_prev_time.get_pt_value(-0.7, 0.1) - (0.014744)) > eps) {
    printf("Coordinate ( -0.7, 0.1) yvel value = %lf\n", yvel_prev_time.get_pt_value(-0.7, 0.1));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(0.65, 0.0) - (-0.033332)) > eps) {
    printf("Coordinate ( 0.65, 0.0) yvel value = %lf\n", yvel_prev_time.get_pt_value(0.65, 0.0));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(0.7, 0.0) - (-0.012353)) > eps) {
    printf("Coordinate ( 0.7, 0.0) yvel value = %lf\n", yvel_prev_time.get_pt_value(0.7, 0.0));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(0.9, 0.0) - (-0.001691)) > eps) {
    printf("Coordinate ( 0.9, 0.0) yvel value = %lf\n", yvel_prev_time.get_pt_value(0.9, 0.0));
    success = 0;
  }
  if (fabs(yvel_prev_time.get_pt_value(1.1, 0.0) - (-0.000986)) > eps) {
    printf("Coordinate ( 1.1, 0.0) yvel value = %lf\n", yvel_prev_time.get_pt_value(1.1, 0.0));
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
