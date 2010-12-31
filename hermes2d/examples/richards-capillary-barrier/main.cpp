#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"


using namespace RefinementSelectors;

//  This example uses adaptivity with dynamical meshes to solve
//  the time-dependent Richard's equation. The time discretization 
//  is backward Euler or Crank-Nicolson, and the nonlinear solver 
//  in each time step is either Newton or Picard. 
//
//  PDE: C(h)dh/dt - div(K(h)grad(h)) - (dK/dh)*(dh/dy) = 0
//  where K(h) = K_S*exp(alpha*h)                          for h < 0,
//        K(h) = K_S                                       for h >= 0,
//        C(h) = alpha*(theta_s - theta_r)*exp(alpha*h)    for h < 0,
//        C(h) = alpha*(theta_s - theta_r)                 for h >= 0.
//
//  Picard's linearization: C(h^k)dh^{k+1}/dt - div(K(h^k)grad(h^{k+1})) - (dK/dh(h^k))*(dh^{k+1}/dy) = 0
//                          Note: the index 'k' does not refer to time stepping.
//  Newton's method is more involved, see the file forms.cpp.
//
//  Domain: rectangle (0, 8) x (0, 6.5).
//
//  BC: Dirichlet, given by the initial condition.
//  IC: See the function init_cond().
//
//  The following parameters can be changed:

// Constitutive relations.
#define CONSTITUTIVE_GENUCHTEN                    // Van Genuchten or Gardner.

// Methods.
const int ITERATIVE_METHOD = 1;		          // 1 = Newton, 2 = Picard.
const int TIME_INTEGRATION = 1;                   // 1 = implicit Euler, 2 = Crank-Nicolson.

// Elements orders and initial refinements.
const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 0;                   // Number of initial mesh refinements towards the top edge.

// Adaptivity.
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See the User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.1;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.


// Newton's and Picard's methods.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for Newton on fine mesh.
int NEWTON_MAX_ITER = 20;                         // Maximum allowed number of Newton iterations.
const double PICARD_TOL = 1e-5;                   // Stopping criterion for Picard on fine mesh.
int PICARD_MAX_ITER = 20;                         // Maximum allowed number of Picard iterations.

// Time stepping and such.
double TAU = 5e-1;                                // Time step.
const double STARTUP_TIME = 5e-2;                 // Start-up time for time-dependent Dirichlet boundary condition.
const double T_FINAL = 1000.0;                    // Time interval length.
double TIME = TAU;                                // Global time variable initialized with first time step.


// Problem parameters.
const char* TABLES_FILENAME = "tables.txt";       // Filename for constitutive tables.
double H_INIT = -50.0;                            // Initial pressure head.
double H_ELEVATION = 10.0;
const double K_S_vals[4] = {350.2, 712.8, 1.68, 18.64}; 
const double ALPHA_vals[4] = {0.01, 1.0, 0.01, 0.01};
const double N_vals[4] = {2.5, 2.0, 1.23, 2.5};
const double M_vals[4] = {0.864, 0.626, 0.187, 0.864};

const double THETA_R_vals[4] = {0.064, 0.0, 0.089, 0.064};
const double THETA_S_vals[4] = {0.14, 0.43, 0.43, 0.24};
const double STORATIVITY_vals[4] = {0.1, 0.1, 0.1, 0.1};

// Precalculation of constitutive tables.
const int MATERIAL_COUNT = 4;
double K_TABLE[4][1500000];                       // Four materials, 15000 values.
double dKdh_TABLE[4][1500000];
double ddKdhh_TABLE[4][1500000];
double C_TABLE[4][1500000];
double dCdh_TABLE[4][1500000];
bool CONSTITUTIVE_TABLES_READY = false;
double*** POLYNOMIALS;                            // Polynomial approximation of the K(h) function close to saturation 
                                                  // (this function has singularity in its second derivative).
const double LOW_LIMIT=-1.0;                      // Lower bound of K(h) function approximated by polynomials.
const int NUM_OF_INSIDE_PTS = 0;
bool POLYNOMIALS_READY = false;
bool POLYNOMIALS_ALLOCATED = false;

// Global variables for forms.
double K_S, ALPHA, THETA_R, THETA_S, N, M, STORATIVITY;

// Choose here which constitutive relations should be used.
#ifdef CONSTITUTIVE_GENUCHTEN
#include "constitutive_genuchten.cpp"
#else
#include "constitutive_gardner.cpp"
#endif

// Boundary markers.
const int BDY_1 = 1;
const int BDY_2 = 2;
const int BDY_3 = 3;
const int BDY_4 = 4;

// Initial condition.
double init_cond(double x, double y, double& dx, double& dy) {
  dx = 0;
  dy = 0;
  return H_INIT;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(double x, double y, double time)
{
  if (time < STARTUP_TIME)
    return H_INIT + time/STARTUP_TIME*(H_ELEVATION-H_INIT);
  else 
    return H_ELEVATION;
}

// Weak forms.
#include "forms.cpp"

// Additional functionality.
#include "extras.cpp"

// Main function.
int main(int argc, char* argv[])
{

  // Either use exact constitutive relations (slow) or precalculate 
  // their polynomial approximations (faster).
  CONSTITUTIVE_TABLES_READY = get_constitutive_tables(TABLES_FILENAME, ITERATIVE_METHOD);

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements.
  mesh.copy(&basemesh);
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_1);
  bc_types.add_bc_neumann(Hermes::Tuple<int>(BDY_2, BDY_3, BDY_4));

  // Enter Dirichlet boundary values.
  BCValues bc_values(&TIME);
  bc_values.add_timedep_function(BDY_1, essential_bc_values);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Solutions for the time stepping and the Newton's method.
  Solution sln, ref_sln, sln_prev_time, sln_prev_iter;
  
  // Assign the function f() to the fine mesh.
  sln_prev_time.set_exact(&mesh, init_cond);
  sln_prev_iter.set_exact(&mesh, init_cond);

  // Initialize the weak formulation.
  WeakForm wf;
  if (ITERATIVE_METHOD == 1) {
    if (TIME_INTEGRATION == 1) {
      info("Registering forms for the Newton's method (implicit Euler in time).");
      wf.add_matrix_form(jac_form_vol_euler, jac_form_vol_ord, HERMES_UNSYM, HERMES_ANY, 
	                 &sln_prev_time);
      wf.add_vector_form(res_form_vol_euler, res_form_vol_ord, HERMES_ANY, 
			 &sln_prev_time);
    }
    else {
      info("Registering forms for the Newton's method (Crank-Nicolson in time).");
      wf.add_matrix_form(jac_form_vol_cranic, jac_form_vol_ord, HERMES_UNSYM, HERMES_ANY, 
      		         &sln_prev_time);
      wf.add_vector_form(res_form_vol_cranic, res_form_vol_ord, HERMES_ANY, 
			 &sln_prev_time);
    }
  }
  else {
    if (TIME_INTEGRATION == 1) {
      info("Registering forms for the Picard's method (implicit Euler in time).");
      wf.add_matrix_form(bilinear_form_picard_euler, bilinear_form_picard_euler_ord, HERMES_UNSYM, HERMES_ANY, 
	                 &sln_prev_iter);
      wf.add_vector_form(linear_form_picard_euler, linear_form_picard_euler_ord, HERMES_ANY, 
			 Hermes::Tuple<MeshFunction*>(&sln_prev_iter, &sln_prev_time));
    }
    else {
      info("Registering forms for the Picard's method (Crank-Nicolson in time).");
      error("Not implemented yet.");
      wf.add_matrix_form(bilinear_form_picard_euler, bilinear_form_picard_euler_ord, HERMES_UNSYM, HERMES_ANY, 
	                 &sln_prev_iter);
      wf.add_vector_form(linear_form_picard_euler, linear_form_picard_euler_ord, HERMES_ANY, 
			 Hermes::Tuple<MeshFunction*>(&sln_prev_iter, &sln_prev_time));
    }
  }

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err_est, graph_time_err_exact, graph_time_dof, graph_time_cpu;
 
  // Visualize the projection and mesh.
  ScalarView view("Initial condition", new WinGeom(0, 0, 630, 350));
  view.fix_scale_width(50);
  OrderView ordview("Initial mesh", new WinGeom(640, 0, 600, 350));
  view.show(&sln_prev_time);
  ordview.show(&space);

  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    info("---- Time step %d:", ts);

    // Time measurement.
    cpu_time.tick();

    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) {
      info("Global mesh derefinement.");
      mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
    }

    // Spatial adaptivity loop. Note: sln_prev_time must not be touched during adaptivity.
    bool done = false;
    int as = 1;
    double err_est_rel;
    // Save time step so that it can be restored after the adaptivity loop,
    // if it needs to be reduced during adaptivity.
    double save_tau = TAU;
    do
    {
      info("---- Time step %d, time %g (days), adaptivity step %d:", ts, TIME, as);

      // Construct globally refined reference mesh
      // and setup reference space.
      Space* ref_space = construct_refined_space(&space);
      ndof = Space::get_num_dofs(ref_space);

      // Next we need to calculate the reference solution.
      // Newton's method:
      if(ITERATIVE_METHOD == 1) {
        scalar* coeff_vec = new scalar[Space::get_num_dofs(ref_space)];
     
        // Calculate initial coefficient vector for Newton on the fine mesh.
        if (as == 1 && ts == 1) {
          info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
          OGProjection::project_global(ref_space, &sln_prev_time, coeff_vec, matrix_solver);
        }
        else {
          info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
          OGProjection::project_global(ref_space, &ref_sln, coeff_vec, matrix_solver);
          delete ref_sln.get_mesh();
        }

        // Initialize the FE problem.
        bool is_linear = false;
        DiscreteProblem dp(&wf, ref_space, is_linear);

        // Set up the solver, matrix, and rhs according to the solver selection.
        SparseMatrix* matrix = create_matrix(matrix_solver);
        Vector* rhs = create_vector(matrix_solver);
        Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

        // Perform Newton's iteration on the reference mesh. If necessary, 
        // reduce time step to make it converge, but then restore time step 
        // size to its original value.
        info("Performing Newton's iteration:");
        bool success, verbose = true;
        double* save_coeff_vec = new double[ndof];
        // Save coefficient vector.
        for (int i=0; i < ndof; i++) save_coeff_vec[i] = coeff_vec[i];
        double damping_coeff = 1.0;
        while (!solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
                             NEWTON_TOL, NEWTON_MAX_ITER, verbose, damping_coeff)) {
          static int failures = 0;
          if (failures >= 10) 
            error("Newton's method did not converge, even with substantial reduction of time step.");
          // Restore solution from the beginning of time step.
          for (int i=0; i < ndof; i++) coeff_vec[i] = save_coeff_vec[i];
          // Reducing time step to 50%.
          warn("Reducing time step size from %g to %g for the rest of this time step.", TAU, TAU *= 0.5);
          TAU *= 0.5;
          // Counting failures.
	  failures++;
        }  
        // Delete the saved coefficient vector.
        delete [] save_coeff_vec;

        // Translate the resulting coefficient vector 
        // into the desired reference solution. 
        Solution::vector_to_solution(coeff_vec, ref_space, &ref_sln);

        // Cleanup.
        delete [] coeff_vec;
        delete solver;
        delete matrix;
        delete rhs;
      }
      else {
        // Calculate initial condition for Picard on the fine mesh.
        if (as == 1 && ts == 1) {
          info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
          OGProjection::project_global(ref_space, &sln_prev_time, &sln_prev_iter, matrix_solver);
        }
        else {
          info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
          OGProjection::project_global(ref_space, &ref_sln, &sln_prev_iter, matrix_solver);
          //delete ref_sln.get_mesh();
        }

        // Perform Picard iteration on the reference mesh. If necessary, 
        // reduce time step to make it converge, but then restore time step 
        // size to its original value.
        info("Performing Picard's iteration:");
        bool success, verbose = true;
        while(!solve_picard(&wf, ref_space, &sln_prev_iter, matrix_solver, PICARD_TOL, 
                            PICARD_MAX_ITER, verbose)) {
          static int failures = 0;
          if (failures >= 10) 
            error("Newton's method did not converge, even with substantial reduction of time step.");
          // Restore solution from the beginning of time step.
          sln_prev_iter.copy(&sln_prev_time);
          // Reducing time step to 50%.
          warn("Reducing time step size from %g to %g for the rest of this time step", TAU, TAU *= 0.5);
          TAU *= 0.5;
          // Counting failures.
	  failures++;
        } 

        ref_sln.copy(&sln_prev_iter);
      }

      /*** ADAPTIVITY ***/

      // Project the fine mesh solution on the coarse mesh.
      info("Projecting fine mesh solution on coarse mesh for error calculation.");
      if(space.get_mesh() == NULL) error("it is NULL");
      OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver);

      // Calculate element errors.
      info("Calculating error estimate."); 
      Adapt* adaptivity = new Adapt(&space, HERMES_H1_NORM);
      bool solutions_for_adapt = true;
      
      // Calculate error estimate wrt. fine mesh solution.
      err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, 
                    HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS) * 100;

      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, space_err_est_rel: %g%%", 
        Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel);

      // If space_err_est too large, adapt the mesh.
      if (err_est_rel < ERR_STOP) done = true;
      else {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (Space::get_num_dofs(&space) >= NDOF_STOP) {
          done = true;
          break;
        }
        as++;
      }

      delete adaptivity;
      delete ref_space;
    }
    while (!done);

    // Updating time step. Note that TAU might have been reduced during adaptivity.
    TIME += TAU;

    // Add entries to convergence graphs.
    graph_time_err_est.add_values(TIME, err_est_rel);
    graph_time_err_est.save("time_error_est.dat");
    graph_time_dof.add_values(TIME, Space::get_num_dofs(&space));
    graph_time_dof.save("time_dof.dat");
    graph_time_cpu.add_values(TIME, cpu_time.accumulated());
    graph_time_cpu.save("time_cpu.dat");

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution, time %g (days)", TIME);
    view.set_title(title);
    view.show(&sln);
    sprintf(title, "Mesh, time %g (days)", TIME);
    ordview.set_title(title);
    ordview.show(&space);
    
    // Save complete Solution.
    char* filename = new char[100];
    sprintf(filename, "outputs/tsln_%f.dat", TIME);
    bool compress = false;   // Gzip compression not used as it only works on Linux.
    sln.save(filename, compress);
    info("Complete Solution saved to file %s.", filename);

    // Copy new reference level solution into sln_prev_time.
    // This starts new time step.
    sln_prev_time.copy(&ref_sln);

    // Restore time step if it was reduced during adaptivity.
    if (fabs(save_tau - TAU) > 1e-12) {
      warn("Restoring time step from %g to %g.", TAU, save_tau);
      TAU = save_tau;
    }
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
