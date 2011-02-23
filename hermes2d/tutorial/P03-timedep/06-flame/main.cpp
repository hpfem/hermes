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
const double time_step = 0.5;                     // Time step.
const double T_FINAL = 60.0;                      // Time interval length.
const double NEWTON_TOL = 1e-4;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 50;                   // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem constants.
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

// Boundary markers.
const int BDY_LEFT = 1, BDY_NEUMANN = 2, BDY_COOLED = 3;

// Initial conditions.
scalar temp_ic(double x, double y, scalar& dx, scalar& dy)
{ return (x <= x1) ? 1.0 : exp(x1 - x); }

scalar conc_ic(double x, double y, scalar& dx, scalar& dy)
{ return (x <= x1) ? 0.0 : 1.0 - exp(Le*(x1 - x)); }

// Weak forms, definition of reaction rate omega.
# include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_LEFT);
  bc_types.add_bc_neumann(BDY_NEUMANN);
  bc_types.add_bc_newton(BDY_COOLED);

  // Enter Dirichlet boundary values.
  BCValues bc_values_t;
  bc_values_t.add_const(BDY_LEFT, 1.0);

  BCValues bc_values_c;
  bc_values_c.add_zero(BDY_LEFT);

  // Create H1 spaces with default shapesets.
  H1Space tspace(&mesh, &bc_types, &bc_values_t, P_INIT);
  H1Space cspace(&mesh, &bc_types, &bc_values_c, P_INIT);
  int ndof = Space::get_num_dofs(Hermes::vector<Space *>(&tspace, &cspace));
  info("ndof = %d.", ndof);

  // Previous time level solutions.
  Solution t_prev_time_1, c_prev_time_1, t_prev_time_2, c_prev_time_2, 
           t_prev_newton, c_prev_newton;

  // And their initial exact setting.
  t_prev_time_1.set_exact(&mesh, temp_ic); c_prev_time_1.set_exact(&mesh, conc_ic);
  t_prev_time_2.set_exact(&mesh, temp_ic); c_prev_time_2.set_exact(&mesh, conc_ic);
  t_prev_newton.set_exact(&mesh, temp_ic); c_prev_newton.set_exact(&mesh, conc_ic);

  // Filters for the reaction rate omega and its derivatives.
  DXDYFilter omega(omega_fn, Hermes::vector<MeshFunction*>(&t_prev_newton, &c_prev_newton));
  DXDYFilter omega_dt(omega_dt_fn, Hermes::vector<MeshFunction*>(&t_prev_newton, &c_prev_newton));
  DXDYFilter omega_dc(omega_dc_fn, Hermes::vector<MeshFunction*>(&t_prev_newton, &c_prev_newton));

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(newton_bilinear_form_0_0), HERMES_NONSYM, HERMES_ANY, &omega_dt);
  wf.add_matrix_form_surf(0, 0, callback(newton_bilinear_form_0_0_surf), 3);
  wf.add_matrix_form(0, 1, callback(newton_bilinear_form_0_1), HERMES_NONSYM, HERMES_ANY, &omega_dc);
  wf.add_matrix_form(1, 0, callback(newton_bilinear_form_1_0), HERMES_NONSYM, HERMES_ANY, &omega_dt);
  wf.add_matrix_form(1, 1, callback(newton_bilinear_form_1_1), HERMES_NONSYM, HERMES_ANY, &omega_dc);
  wf.add_vector_form(0, callback(newton_linear_form_0), HERMES_ANY, 
                     Hermes::vector<MeshFunction*>(&t_prev_time_1, &t_prev_time_2, &omega));
  wf.add_vector_form_surf(0, callback(newton_linear_form_0_surf), 3);
  wf.add_vector_form(1, callback(newton_linear_form_1), HERMES_ANY, 
                     Hermes::vector<MeshFunction*>(&c_prev_time_1, &c_prev_time_2, &omega));

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&tspace, &cspace), is_linear);

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
  double current_time = 0.0; int ts = 1;
  do 
  {
    info("---- Time step %d, t = %g s.", ts, current_time);

    // Perform Newton's iteration. Note that we are not using the 
    // function solve_newton() now since extra actions are needed
    // inside the loop.
    int it = 1;
    while (1)
    {
      dp.assemble(coeff_vec, matrix, rhs, false);
      
      // Multiply the residual vector with -1 since the matrix 
      // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
      rhs->change_sign();
      
      // Calculate the l2-norm of residual vector.
      double res_l2_norm = get_l2_norm(rhs);

      // Info for user.
      info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, 
           Space::get_num_dofs(Hermes::vector<Space *>(&tspace, &cspace)), res_l2_norm);

      // If l2 norm of the residual vector is within tolerance, or the maximum number 
      // of iteration has been reached, then quit.
      if (res_l2_norm < NEWTON_TOL || it > NEWTON_MAX_ITER) break;

      // Solve the linear system and if successful, obtain the solutions.
      if(!solver->solve())
        error ("Matrix solver failed.\n");

        // Add \deltaY^{n+1} to Y^n.
      for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];
      
      if (it >= NEWTON_MAX_ITER)
        error ("Newton method did not converge.");
     
      // Set current solutions to the latest Newton iterate 
      // and reinitialize filters of these solutions.
      Solution::vector_to_solutions(coeff_vec, Hermes::vector<Space *>(&tspace, &cspace), 
                                    Hermes::vector<Solution *>(&t_prev_newton, &c_prev_newton));
      omega.reinit();
      omega_dt.reinit();
      omega_dc.reinit();

      it++;
    };

    // Visualization.
    DXDYFilter omega_view(omega_fn, Hermes::vector<MeshFunction*>(&t_prev_newton, &c_prev_newton));
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

    ts++;
  } 
  while (current_time <= T_FINAL);

  // Cleanup.
  delete [] coeff_vec;
  delete matrix;
  delete rhs;
  delete solver;
  
  // Wait for all views to be closed.
  View::wait();
  return 0;
}
