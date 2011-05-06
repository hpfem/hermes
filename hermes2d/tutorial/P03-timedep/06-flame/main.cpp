#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example is a very simple flame propagation model (laminar flame,
//  zero flow velocity), and its purpose is to show how the Newton's method
//  is applied to a time-dependent two-equation system.
//
//  PDEs:
//
//  dT/dt - laplace T = omega(T,Y),
//  dY/dt - 1/Le * laplace Y = - omega(T,Y).
//
//  Domain: rectangle with cooled rods.
//
//  BC:  T = 1, Y = 0 on the inlet,
//       dT/dn = - kappa T on cooled rods,
//       dT/dn = 0, dY/dn = 0 elsewhere.
//
//  Time-stepping: a second order BDF formula.

const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
const int P_INIT = 1;                             // Initial polynomial degree.

const double T_FINAL = 60.0;                      // Time interval length.
const double NEWTON_TOL = 1e-4;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 50;                   // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
const bool NEWTON = true;                         // If NEWTON == true then the Newton's iteration is performed.
                                                  // in every time step. Otherwise the convective term is linearized
                                                  // using the velocities from the previous time step.
// Problem constants.
const double time_step = 0.5;                     // Time step.
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

// Boundary markers.
const std::string BDY_LEFT = "Bdy_left", BDY_NEUMANN = "Bdy_neumann", BDY_NEWTON_COOLED = "Bdy_newton_cooled";

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

/*
  // Initialize boundary conditions.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_LEFT);
  bc_types.add_bc_neumann(BDY_NEUMANN);
  bc_types.add_bc_newton(BDY_NEWTON_COOLED);

  // Enter Dirichlet boundary values.
  BCValues bc_values_t;
  bc_values_t.add_const(BDY_LEFT, 1.0);

  BCValues bc_values_c;
  bc_values_c.add_zero(BDY_LEFT);
*/

  // Initialize boundary conditions.
  DefaultEssentialBCConst bc_t(BDY_LEFT, 1.0);
  EssentialBCs bcs_t(&bc_t);
  DefaultEssentialBCConst bc_c(BDY_LEFT, 0.0);
  EssentialBCs bcs_c(&bc_c);

  // Create H1 spaces with default shapesets.
  H1Space tspace(&mesh, &bcs_t, P_INIT);
  H1Space cspace(&mesh, &bcs_c, P_INIT);
  int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&tspace, &cspace));
  info("ndof = %d.", ndof);

  // Previous time level solutions and their initial exact setting.
  InitialSolutionTemperature t_prev_time_1(&mesh, x1);
  InitialSolutionTemperature t_prev_time_2(&mesh, x1);
  InitialSolutionTemperature t_prev_newton(&mesh, x1);

  InitialSolutionConcentration c_prev_time_1(&mesh, x1);
  InitialSolutionConcentration c_prev_time_2(&mesh, x1);
  InitialSolutionConcentration c_prev_newton(&mesh, x1);

  // Filters for the reaction rate omega and its derivatives.
  DXDYFilterOmega    omega(Hermes::vector<MeshFunction*>(&t_prev_newton, &c_prev_newton));
  DXDYFilterOmega_dt omega_dt(Hermes::vector<MeshFunction*>(&t_prev_newton, &c_prev_newton));
  DXDYFilterOmega_dc omega_dc(Hermes::vector<MeshFunction*>(&t_prev_newton, &c_prev_newton));

  bool JFNK = false;
  // Initialize the weak formulation.
  CustomWeakForm wf(BDY_NEUMANN, BDY_NEWTON_COOLED, 
                    time_step, Le, kappa, 
                    &omega_dt, &omega_dc, &omega, 
                    &t_prev_time_1, &t_prev_time_2, &t_prev_newton, 
                    &c_prev_time_1, &c_prev_time_2, &c_prev_newton, JFNK);

  // Initialize the FE problem.
  //bool is_linear = false;
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&tspace, &cspace));

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  scalar* coeff_vec = new scalar[ndof];
  OGProjection::project_global(Hermes::vector<Space *>(&tspace, &cspace), 
                               Hermes::vector<MeshFunction *>(&t_prev_newton, &c_prev_newton), 
                               coeff_vec, matrix_solver);

  // Initialize views.
  ScalarView rview("Reaction rate", new WinGeom(0, 0, 800, 230));

  // Time stepping loop:
  double current_time = 0.0; 
  int num_time_steps = T_FINAL / time_step;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    current_time += time_step;
    info("---- Time step %d, time = %g:", ts, current_time);

    // Update time-dependent essential BCs.
    //if (current_time <= STARTUP_TIME) {
    //  info("Updating time-dependent essential BC.");
      Space::update_essential_bc_values(Hermes::vector<Space *>(&tspace, &cspace), current_time);
    //}

    if (NEWTON)
    {
      // Perform Newton's iteration.
      info("Solving nonlinear problem:");
      bool verbose = true;
      bool jacobian_changed = true;
      if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, jacobian_changed,
          NEWTON_TOL, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

      // Update previous time level solutions.
      Solution::vector_to_solutions(coeff_vec, Hermes::vector<Space *>(&tspace, &cspace),
                                    Hermes::vector<Solution *>(&t_prev_time_1, &c_prev_time_1));
    }
    else {
      // Linear solve.
      info("Assembling and solving linear problem.");
      dp.assemble(matrix, rhs, false);
      if(solver->solve())
        Solution::vector_to_solutions(solver->get_solution(),
                  Hermes::vector<Space *>(&tspace, &cspace),
                  Hermes::vector<Solution *>(&t_prev_time_1, &c_prev_time_1));
      else
        error ("Matrix solver failed.\n");
    }


    // Visualization.
    DXDYFilterOmega omega_view(Hermes::vector<MeshFunction*>(&t_prev_newton, &c_prev_newton));
    rview.set_min_max_range(0.0,2.0);
    char title[100];
    sprintf(title, "Reaction rate, t = %g", current_time);
    rview.set_title(title);
    rview.show(&omega_view);

    // Update current time.
    current_time += time_step;

    // Store two time levels of previous solutions.
    t_prev_time_2.copy(&t_prev_time_1);
    c_prev_time_2.copy(&c_prev_time_1);
    Solution::vector_to_solutions(coeff_vec, Hermes::vector<Space *>(&tspace, &cspace), 
                                  Hermes::vector<Solution *>(&t_prev_time_1, &c_prev_time_1));

  } 

  // Cleanup.
  delete [] coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;
  
  // Wait for all views to be closed.
  View::wait();
  return 0;
}
