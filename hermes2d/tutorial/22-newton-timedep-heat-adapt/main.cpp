#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

//  This example shows how to combine automatic adaptivity with the Newton's
//  method for a nonlinear time-dependent PDE discretized implicitly in time
//  using implicit Euler or Crank-Nicolson.
//
//  PDE: time-dependent heat transfer equation with nonlinear thermal
//  conductivity, du/dt - div[lambda(u)grad u] = f.
//
//  Domain: square (-10,10)^2.
//
//  BC:  Dirichlet, given by the function dir_lift() below.
//  IC: Same function dir_lift().
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 2;                // Number of initial uniform mesh refinements.
const int P_INIT = 2;                      // Initial polynomial degree of all mesh elements.
const int TIME_DISCR = 2;                  // 1 for implicit Euler, 2 for Crank-Nicolson.
const double TAU = 0.5;                    // Time step.
const double T_FINAL = 5.0;                // Time interval length.

// Adaptivity
const int UNREF_FREQ = 1;                  // Every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_H; // Predefined list of element refinement candidates. Possible values are
                                           // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                           // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                           // See the User Documentation for details.
const int MESH_REGULARITY = -1;            // Maximum allowed level of hanging nodes:
                                           // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                           // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                           // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                           // Note that regular meshes are not supported, this is due to
                                           // their notoriously bad performance.
const double CONV_EXP = 1.0;               // Default value is 1.0. This parameter influences the selection of
                                           // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 3.0;               // Stopping criterion for adaptivity (rel. error tolerance between the
                                           // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Newton's method
const double NEWTON_TOL_COARSE = 0.01;     // Stopping criterion for Newton on fine mesh.
const double NEWTON_TOL_FINE = 0.05;       // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;           // Maximum allowed number of Newton iterations.

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u)
{
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) {
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions.
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// Initial condition.
scalar init_cond(double x, double y, double& dx, double& dy)
{
  return dir_lift(x, y, dx, dy);
}

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Heat sources (can be a general function of 'x' and 'y').
template<typename Real>
Real heat_src(Real x, Real y)
{
  return 1.0;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("square.mesh", &basemesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) basemesh.refine_all_elements();
  mesh.copy(&basemesh);

  // Create an H1 space with default shapeset.
  H1Space* space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(space);

  // Convert initial condition into a Solution.
  Solution sln_prev_time;
  sln_prev_time.set_exact(&mesh, init_cond);

  // Initialize the weak formulation.
  WeakForm wf;
  if(TIME_DISCR == 1) {
    wf.add_matrix_form(callback(J_euler), H2D_UNSYM, H2D_ANY);
    wf.add_vector_form(callback(F_euler), H2D_ANY, &sln_prev_time);
  }
  else {
    wf.add_matrix_form(callback(J_cranic), H2D_UNSYM, H2D_ANY);
    wf.add_vector_form(callback(F_cranic), H2D_ANY, &sln_prev_time);
  }

  // Initialize adaptivity parameters.
  AdaptivityParamType apt(ERR_STOP, NDOF_STOP, THRESHOLD, STRATEGY, MESH_REGULARITY);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Project initial condition to coarse mesh. 
  Vector *coeff_vec = new AVector();
  info("Projecting initial condition to obtain coefficient vector on coarse mesh.");

  project_global(space, H2D_H1_NORM, &sln_prev_time, Tuple<Solution*>(), coeff_vec);

  // Visualize initial condition.
  char title[100];
  ScalarView view("Initial condition", new WinGeom(0, 0, 440, 350));
  OrderView ordview("Initial mesh", new WinGeom(450, 0, 400, 350));
  view.show(&sln_prev_time);
  ordview.show(space);

  // Newton's method on coarse mesh (moving one time step forward).
  info("Solving on coarse mesh.");
  bool verbose = true; // Default is false.
  if (!solve_newton(space, &wf, coeff_vec, matrix_solver, 
                    NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
    error("Newton's method did not converge.");
  
  // Set initial coarse mesh solution, create variable for the reference solution.
  Solution sln, ref_sln;
  sln.set_fe_solution(space, coeff_vec);
  
  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      space->set_uniform_order(P_INIT);

      // Project on globally derefined mesh.
      info("Projecting previous fine mesh solution on derefined mesh.");
      project_global(space, H2D_H1_NORM, &ref_sln, Tuple<Solution *>(), coeff_vec);

      /*
      // Newton's method on derefined mesh (moving one time step forward).
      info("Solving on derefined mesh.");
      bool verbose = true; // Default is false.
      if (!solve_newton(space, &wf, coeff_vec, matrix_solver, 
                        NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose))
        error("Newton's method did not converge.");
      */
      
      sln.set_fe_solution(space, coeff_vec);
    }

    // Adaptivity loop:
    bool done = false; int as = 1;
    double err_est;
    do {
      info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh
      // and setup reference space.
      Space* ref_space;
      Mesh* ref_mesh = new Mesh();
      ref_mesh->copy(space->get_mesh());
      ref_mesh->refine_all_elements();
      ref_space = space->dup(ref_mesh);
      int order_increase = 1;
      ref_space->copy_orders(space, order_increase);

      // Calculate initial coefficient vector for Newton on the fine mesh.
      if (as == 1) {
        info("Projecting coarse mesh solution to obtain coefficient vector on new fine mesh.");
        project_global(ref_space, H2D_H1_NORM, &sln, Tuple<Solution*>(), coeff_vec);
      }
      else {
        info("Projecting previous fine mesh solution to obtain coefficient vector on new fine mesh.");
        project_global(ref_space, H2D_H1_NORM, &ref_sln, Tuple<Solution*>(), coeff_vec);
      }

      // Newton's method on fine mesh
      info("Solving on fine mesh.");
      bool verbose = true; // Default is false.
      if (!solve_newton(ref_space, &wf, coeff_vec, matrix_solver, 
                        NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose))
        error("Newton's method did not converge.");

      // Store the result in ref_sln.
      ref_sln.set_fe_solution(ref_space, coeff_vec);

      // Calculate element errors.
      info("Calculating error (est).");
      Adapt hp(space, H2D_H1_NORM);
      // Pass coarse mesh and reference solutions for error estimation.
      hp.set_solutions(&sln, &ref_sln);
      double err_est_rel_total = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100.;

      // Report results.
      info("ndof: %d, ref_ndof: %d, err_est_rel: %g%%", 
           get_num_dofs(space), get_num_dofs(ref_space), err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
      else {
        if (verbose) info("Adapting the coarse mesh.");
        done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

        if (get_num_dofs(space) >= NDOF_STOP) {
          done = true;
          break;
        }

        info("Projecting fine mesh solution on new coarse mesh.");
        // The NULL pointer means that we do not want the resulting coefficient vector.
        project_global(space, H2D_H1_NORM, &ref_sln, &sln, NULL);
      }

      // Free the reference space and mesh.
      ref_space->free();
      ref_mesh->free();

      as++;
    }
    while (done == false);

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution, time level %d", ts);
    view.set_title(title);
    view.show(&sln);
    sprintf(title, "Mesh, time level %d", ts);
    ordview.set_title(title);
    ordview.show(space);

    // Copy last reference solution into sln_prev_time.
    sln_prev_time.copy(&ref_sln);
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
