#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
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

const int INIT_REF_NUM = 2;            // Number of initial uniform mesh refinements.
const int P_INIT = 1;                  // Initial polynomial degree.
const double TAU = 0.5;                // Time step.
const double T_FINAL = 60.0;           // Time interval length.
const double NEWTON_TOL = 1e-4;        // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 50;        // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem constants.
const double Le    = 1.0;
const double alpha = 0.8;
const double beta  = 10.0;
const double kappa = 0.1;
const double x1    = 9.0;

// Boundary markers.
int bdy_left = 1;

// Boundary conditions.
BCType bc_types(int marker)
  { return (marker == bdy_left) ? BC_ESSENTIAL : BC_NATURAL; }

scalar essential_bc_values_t(int ess_bdy_marker, double x, double y)
  { return (ess_bdy_marker == bdy_left) ? 1.0 : 0; }

scalar essential_bc_values_c(int ess_bdy_marker, double x, double y)
  { return 0; }

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

  // Create H1 spaces with default shapesets.
  H1Space* tspace = new H1Space(&mesh, bc_types, essential_bc_values_t, P_INIT);
  H1Space* cspace = new H1Space(&mesh, bc_types, essential_bc_values_c, P_INIT);
  int ndof = get_num_dofs(Tuple<Space *>(tspace, cspace));
  info("ndof = %d.", ndof);

  // Previous time level solutions.
  Solution t_prev_time_1, c_prev_time_1, t_prev_time_2, c_prev_time_2, 
           t_prev_newton, c_prev_newton;

  // And their initial exact setting.
  t_prev_time_1.set_exact(&mesh, temp_ic); c_prev_time_1.set_exact(&mesh, conc_ic);
  t_prev_time_2.set_exact(&mesh, temp_ic); c_prev_time_2.set_exact(&mesh, conc_ic);
  t_prev_newton.set_exact(&mesh, temp_ic); c_prev_newton.set_exact(&mesh, conc_ic);

  // Projecting initial conditions to obtain initial vector for the Newton's method.
  info("Projecting initial conditions to obtain initial vector for the Newton's method.");
  Vector* coeff_vec = new AVector(); 
  project_global(Tuple<Space *>(tspace, cspace), Tuple<int>(H2D_H1_NORM, H2D_H1_NORM),
  Tuple<MeshFunction*>(&t_prev_newton, &c_prev_newton),Tuple<Solution*>(), coeff_vec);

  // Filters for the reaction rate omega and its derivatives.
  DXDYFilter omega(omega_fn, Tuple<MeshFunction*>(&t_prev_newton, &c_prev_newton));
  DXDYFilter omega_dt(omega_dt_fn, Tuple<MeshFunction*>(&t_prev_newton, &c_prev_newton));
  DXDYFilter omega_dc(omega_dc_fn, Tuple<MeshFunction*>(&t_prev_newton, &c_prev_newton));

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(newton_bilinear_form_0_0), H2D_UNSYM, H2D_ANY, &omega_dt);
  wf.add_matrix_form_surf(0, 0, callback(newton_bilinear_form_0_0_surf), 3);
  wf.add_matrix_form(0, 1, callback(newton_bilinear_form_0_1), H2D_UNSYM, H2D_ANY, &omega_dc);
  wf.add_matrix_form(1, 0, callback(newton_bilinear_form_1_0), H2D_UNSYM, H2D_ANY, &omega_dt);
  wf.add_matrix_form(1, 1, callback(newton_bilinear_form_1_1), H2D_UNSYM, H2D_ANY, &omega_dc);
  wf.add_vector_form(0, callback(newton_linear_form_0), H2D_ANY, 
                     Tuple<MeshFunction*>(&t_prev_time_1, &t_prev_time_2, &omega));
  wf.add_vector_form_surf(0, callback(newton_linear_form_0_surf), 3);
  wf.add_vector_form(1, callback(newton_linear_form_1), H2D_ANY, 
                     Tuple<MeshFunction*>(&c_prev_time_1, &c_prev_time_2, &omega));

  // Initialize view.
  ScalarView rview("Reaction rate", new WinGeom(0, 0, 800, 230));

  // Time stepping loop:
  double current_time = 0.0; int ts = 1;
  do 
  {
    info("---- Time step %d, t = %g s.", ts, current_time);

    // Newton's method.
    info("Performing Newton's method.");
   
    bool verbose = true; // Default is false.

    int ndof = get_num_dofs(Tuple<Space *>(tspace, cspace));
    int n_mesh_fns = 3;

    // Initialize the discrete problem.
    DiscreteProblem dp(&wf, Tuple<Space *>(tspace, cspace));
    //info("ndof = %d", dp.get_num_dofs());

    // Select matrix solver.
    Matrix* mat; Vector* rhs; CommonSolver* solver;
		
    init_matrix_solver(matrix_solver, ndof, mat, rhs, solver, false);
		
    // Newton's loop.
    int it = 1;
    for (; it <= NEWTON_MAX_ITER; it++)
    {
      // Assemble the Jacobian matrix and residual vector.
      bool rhsonly = false;

      // The NULL stands for the dir vector which is not needed here.
      dp.assemble(coeff_vec, mat, NULL, rhs, rhsonly, false);
			
      // Multiply the residual vector with -1 since the matrix 
      // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
      for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));
			
      // Calculate the l2-norm of residual vector.
      double res_l2_norm = 0;
      for (int i = 0; i < ndof; i++) res_l2_norm += rhs->get(i)*rhs->get(i);
      res_l2_norm = sqrt(res_l2_norm);

      if (verbose) info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, ndof, res_l2_norm);

      // If l2 norm of the residual vector is in tolerance, quit.
      if (res_l2_norm < NEWTON_TOL) 
        break;

      // Solve the matrix problem.
      if (!solver->solve(mat, rhs)) 
        error ("Matrix solver failed.\n");

      // Add \deltaY^{n+1} to Y^n.
      for (int i = 0; i < ndof; i++) coeff_vec->add(i, rhs->get(i));
      
      // Set current solutions to the latest Newton iterate and reinitialize filters of these solutions.
      t_prev_newton.set_fe_solution(tspace, coeff_vec);
      c_prev_newton.set_fe_solution(cspace, coeff_vec);
      omega.reinit();
      omega_dt.reinit();
      omega_dc.reinit();
    };

    delete rhs;
    delete mat;
    delete solver;

    // If we have not reached the desired tolerance.
    if(it == NEWTON_MAX_ITER)
       error("Newton's method did not converge.");

    // Visualization.
    DXDYFilter omega_view(omega_fn, Tuple<MeshFunction*>(&t_prev_newton, &c_prev_newton));
    rview.set_min_max_range(0.0,2.0);
    char title[100];
    sprintf(title, "Reaction rate, t = %g", current_time);
    rview.set_title(title);
    rview.show(&omega_view);

    // Update current time.
    current_time += TAU;

    // Store two time levels of previous solutions.
    t_prev_time_2.copy(&t_prev_time_1);
    c_prev_time_2.copy(&c_prev_time_1);
    t_prev_time_1.set_fe_solution(tspace, coeff_vec);
    c_prev_time_1.set_fe_solution(cspace, coeff_vec);

    ts++;
  } 
  while (current_time <= T_FINAL);

  delete coeff_vec;
  
  // Wait for all views to be closed.
  View::wait();
  return 0;
}
