#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// Flow in between two circles, inner circle is rotating with surface 
// velocity VEL. The time-dependent laminar incompressible Navier-Stokes equations
// are discretized in time via the implicit Euler method. If NEWTON == true,
// the Newton's method is used to solve the nonlinear problem at each time
// step. If NEWTON == false, the convective term is only linearized using the
// velocities from the previous time step. Obviously the latter approach is wrong,
// but people do this frequently because it is faster and simpler to implement.
// Therefore we include this case for comparison purposes. We also show how
// to use discontinuous ($L^2$) elements for pressure and thus make the
// velocity discretely divergence-free. Comparison to approximating the
// pressure with the standard (continuous) Taylor-Hood elements is enabled.
// The Reynolds number Re = 200 which is very low. You can increase it but 
// then you will need to make the mesh finer, and the computation will take 
// more time.
//
// PDE: incompressible Navier-Stokes equations in the form
//     \partial v / \partial t - \Delta v / Re + (v \cdot \nabla) v + \nabla p = 0,
//     div v = 0.
//
// BC: tangential velocity V on Gamma_1 (outer circle)
//
// Geometry: Area in between two concentric circles with radiuses r1 and r2,
//           r1 < r2. These radiuses can be changed in the file "domain.mesh".
//
// The following parameters can be changed:

const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements. 
const int INIT_BDY_REF_NUM_INNER = 2;             // Number of initial mesh refinements towards boundary. 

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
const double T_FINAL = 3600.0;                    // Time interval length.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 10;                   // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Current time (used in weak forms).
double current_time = 0;

// Weak forms.
#include "definitions.cpp"

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
  mloader.load("domain.mesh", &mesh);

  //mloader.load("domain-concentric.mesh", &mesh);

  // Initial mesh refinements.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(Hermes::vector<std::string>("Bdy-1", "Bdy-2", "Bdy-3", "Bdy-4"), 
                               INIT_BDY_REF_NUM_INNER, false);  // true for anisotropic refinement

  // Initialize boundary conditions.
  EssentialBCNonConstX bc_inner_vel_x(Hermes::vector<std::string>("Bdy-1", "Bdy-2", "Bdy-3","Bdy-4"), VEL, STARTUP_TIME);
  EssentialBCNonConstY bc_inner_vel_y(Hermes::vector<std::string>("Bdy-1", "Bdy-2", "Bdy-3","Bdy-4"), VEL, STARTUP_TIME);
  EssentialBCs bcs_vel_x(&bc_inner_vel_x);
  EssentialBCs bcs_vel_y(&bc_inner_vel_y);

  // Spaces for velocity components and pressure.
  H1Space xvel_space(&mesh, &bcs_vel_x, P_INIT_VEL);
  H1Space yvel_space(&mesh, &bcs_vel_y, P_INIT_VEL);
#ifdef PRESSURE_IN_L2
  L2Space p_space(&mesh, P_INIT_PRESSURE);
#else
  H1Space p_space(&mesh, P_INIT_PRESSURE);
#endif
  Hermes::vector<Space *> spaces = Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space);

  // Calculate and report the number of degrees of freedom.
  int ndof = Space::get_num_dofs(spaces);
  info("ndof = %d.", ndof);

  // Define projection norms.
  ProjNormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  ProjNormType p_proj_norm = HERMES_L2_NORM;
#else
  ProjNormType p_proj_norm = HERMES_H1_NORM;
#endif

  // Solutions for the Newton's iteration and time stepping.
  info("Setting initial conditions.");
  Solution xvel_prev_time, yvel_prev_time, p_prev_time; 
  xvel_prev_time.set_zero(&mesh);
  yvel_prev_time.set_zero(&mesh);
  p_prev_time.set_zero(&mesh);
  Hermes::vector<Solution*> slns = Hermes::vector<Solution*>(&xvel_prev_time, &yvel_prev_time, 
                                                             &p_prev_time);

  // Initialize weak formulation.
  WeakForm* wf;
  if (NEWTON)
    wf = new WeakFormNSNewton(STOKES, RE, TAU, &xvel_prev_time, &yvel_prev_time);
  else
    wf = new WeakFormNSSimpleLinearization(STOKES, RE, TAU, &xvel_prev_time, &yvel_prev_time);

  // Initialize the FE problem.
  DiscreteProblem dp(wf, spaces);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize views.
  VectorView vview("velocity [m/s]", new WinGeom(0, 0, 600, 500));
  ScalarView pview("pressure [Pa]", new WinGeom(610, 0, 600, 500));
  //vview.set_min_max_range(0, 1.6);
  vview.fix_scale_width(80);
  //pview.set_min_max_range(-0.9, 1.0);
  pview.fix_scale_width(80);
  pview.show_mesh(true);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  scalar* coeff_vec = new scalar[Space::get_num_dofs(spaces)];
  if (NEWTON) {
    // Newton's vector is set to zero (no OG projection needed).
    memset(coeff_vec, 0, ndof * sizeof(double));
    /*
    // This can be used for more complicated initial conditions.
      info("Projecting initial condition to obtain initial vector for the Newton's method.");
      OGProjection::project_global(spaces, slns, coeff_vec, matrix_solver, 
                                   Hermes::vector<ProjNormType>(vel_proj_norm, vel_proj_norm, p_proj_norm));
    */
  }

  // Time-stepping loop:
  char title[100];
  int num_time_steps = T_FINAL / TAU;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    current_time += TAU;
    info("---- Time step %d, time = %g:", ts, current_time);

    // Update time-dependent essential BCs.
    info("Updating time-dependent essential BC.");
    Space::update_essential_bc_values(Hermes::vector<Space *>(&xvel_space, &yvel_space), current_time);

    if (NEWTON) 
    {
      // Perform Newton's iteration.
      info("Solving nonlinear problem:");
      bool verbose = true;
      bool jacobian_changed = true;
      if (!hermes_2D.solve_newton(coeff_vec, &dp, solver, matrix, rhs, jacobian_changed,
          NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

      // Update previous time level solutions.
      Solution::vector_to_solutions(coeff_vec, spaces, slns);
    }
    else {
      // Linear solve.
      info("Assembling and solving linear problem.");
      dp.assemble(matrix, rhs, false);
      if(solver->solve()) 
        Solution::vector_to_solutions(solver->get_solution(), 
                  Hermes::vector<Space *>(&xvel_space, &yvel_space, &p_space), 
                  Hermes::vector<Solution *>(&xvel_prev_time, &yvel_prev_time, &p_prev_time));
      else 
        error ("Matrix solver failed.\n");
    }

    // Show the solution at the end of time step.
    sprintf(title, "Velocity, time %g", current_time);
    vview.set_title(title);
    vview.show(&xvel_prev_time, &yvel_prev_time);
    sprintf(title, "Pressure, time %g", current_time);
    pview.set_title(title);
    pview.show(&p_prev_time);
 }

  // Clean up.
  delete [] coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
