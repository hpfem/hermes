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
//  Units: length: cm
//         time: days
//
//
//  BC: Dirichlet, given by the initial condition.
//  IC: See the function init_cond().
//
//  The following parameters can be changed:

// Constitutive relations.
#define CONSTITUTIVE_GENUCHTEN                    // Van Genuchten or Gardner.

// Choose full domain or half domain.
//const char* mesh_file = "domain-full.mesh";
const char* mesh_file = "domain-half.mesh";

// Methods.
const int ITERATIVE_METHOD = 2;		          // 1 = Newton, 2 = Picard.
const int TIME_INTEGRATION = 1;                   // 1 = implicit Euler, 2 = Crank-Nicolson.

// Adaptive time stepping.
double time_step = 0.5;                           // Time step (in days).
double time_step_dec = 0.5;                       // Timestep decrease ratio after unsuccessful nonlinear solve.
double time_step_inc = 1.1;                       // Timestep increase ratio after successful nonlinear solve.
double time_step_min = 1e-8; 			  // Computation will stop if time step drops below this value. 
double time_step_max = 1.0;                       // Maximal time step.

// Elements orders and initial refinements.
const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY_TOP = 0;               // Number of initial mesh refinements towards the top edge.

// Adaptivity.
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is unrefined.
const int UNREF_LEVEL = 1;                        // 1 = one layer of refinements is shaved off and poly degrees
                                                  // of all elements reset to P_INIT; 2 = mesh reset to basemesh.  
                                                  // TODO: Add a third option where one layer will be taken off 
                                                  // and just one polynomial degree subtracted.
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
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.


// Newton's and Picard's methods.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for Newton on fine mesh.
int NEWTON_MAX_ITER = 10;                         // Maximum allowed number of Newton iterations.
const double PICARD_TOL = 1e-2;                   // Stopping criterion for Picard on fine mesh.
int PICARD_MAX_ITER = 23;                         // Maximum allowed number of Picard iterations.


// Times.
const double STARTUP_TIME = 5.0;                  // Start-up time for time-dependent Dirichlet boundary condition.
const double T_FINAL = 1000.0;                    // Time interval length.
const double PULSE_END_TIME = 1000.0;             // Time interval of the top layer infiltration.
double current_time = time_step;                  // Global time variable initialized with first time step.



// Problem parameters.
double H_INIT = -15.0;                            // Initial pressure head.
double H_ELEVATION = 10.0;                        // Top constant pressure head -- an infiltration experiment.
const double K_S_vals[4] = {350.2, 712.8, 1.68, 18.64}; 
const double ALPHA_vals[4] = {0.01, 1.0, 0.01, 0.01};
const double N_vals[4] = {2.5, 2.0, 1.23, 2.5};
const double M_vals[4] = {0.864, 0.626, 0.187, 0.864};

const double THETA_R_vals[4] = {0.064, 0.0, 0.089, 0.064};
const double THETA_S_vals[4] = {0.14, 0.43, 0.43, 0.24};
const double STORATIVITY_vals[4] = {0.1, 0.1, 0.1, 0.1};


// Precalculation of constitutive tables.
const int MATERIAL_COUNT = 4;
const int CONSTITUTIVE_TABLE_METHOD = 2;          // 0 - constitutive functions are evaluated directly (slow performance).
					          // 1 - constitutive functions are linearly approximated on interval <TABLE_LIMIT; LOW_LIMIT> 
						  //	 (very efficient CPU utilization less efficient memory consumption (depending on TABLE_PRECISION)).
						  // 2 - constitutive functions are aproximated by quintic splines.
						  
//!Use only if 	CONSTITUTIVE_TABLE_METHOD == 2 !//					  
const int NUM_OF_INTERVALS = 16;                  // Number of intervals.                      
const double INTERVALS_4_APPROX[16] = {-1.0, -2.0, -3.0, -4.0, -5.0, -8.0, -10.0, -12.0, // Low limits of intervals approximated by quintic splines.
				     -15.0, -20.0, -30.0, -50.0, -75.0, -100.0,-300.0, -1000.0}; 
int* POL_SEARCH_HELP;                             // This array contains for each integer of h function appropriate polynomial ID.
double**** K_POLS;                                // First DIM is the interval ID, second DIM is the material ID, third DIM is 
                                                  // the derivative degree, fourth DIM are the coefficients.
double**** C_POLS;
//!------------------------------------------!//


//!Use only if CONSTITUTIVE_TABLE_METHOD == 1 !//
double TABLE_LIMIT = -1000.0; 		          // Limit of precalculated functions (should be always negative value lower 
						  // then the lowest expect value of the solution (consider DMP!!)
const double TABLE_PRECISION = 0.1;               // Precision of precalculated table use 1.0, 0,1, 0.01, etc.....
double** K_TABLE;                                  
double** dKdh_TABLE;
double** ddKdhh_TABLE;
double** C_TABLE;
double** dCdh_TABLE;
bool CONSTITUTIVE_TABLES_READY = false;
double*** POLYNOMIALS;                            // Polynomial approximation of the K(h) function close to saturation.
                                                  // This function has singularity in its second derivative.
						  // First dimension is material ID
						  // Second dimension is the polynomial derivative.
						  // Third dimension are the polynomial's coefficients.
						  
						  
const double LOW_LIMIT=-1.0;                      // Lower bound of K(h) function approximated by polynomials.
const int NUM_OF_INSIDE_PTS = 0;
//!------------------------------------------!//


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
const int BDY_TOP = 1;
const int BDY_RIGHT = 2;
const int BDY_BOTTOM = 3;
const int BDY_LEFT = 4;

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
  else if (time > PULSE_END_TIME)
    return H_INIT;
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
  
  // Either use exact constitutive relations (slow) (method 0) or precalculate 
  // their linear approximations (faster) (method 1) or
  // precalculate their quintic polynomial approximations (method 2) -- managed by 
  // the following loop "Initializing polynomial approximation".
  if (CONSTITUTIVE_TABLE_METHOD == 1)
    CONSTITUTIVE_TABLES_READY = get_constitutive_tables(ITERATIVE_METHOD);
  // Points to be used for polynomial approximation of K(h).
  double* points = new double[NUM_OF_INSIDE_PTS];

  // The van Genuchten + Mualem K(h) function is approximated by polynomials close 
  // to zero in case of CONSTITUTIVE_TABLE_METHOD==1.
  // In case of CONSTITUTIVE_TABLE_METHOD==2, all constitutive functions are approximated by polynomials.
  info("Initializing polynomial approximations.");
  for (int i=0; i < MATERIAL_COUNT; i++) {
    init_polynomials(6 + NUM_OF_INSIDE_PTS, LOW_LIMIT, points, NUM_OF_INSIDE_PTS, i);
  }
  POLYNOMIALS_READY = true;
  if (CONSTITUTIVE_TABLE_METHOD == 2) {
    CONSTITUTIVE_TABLES_READY = true ;
    //Assign table limit to global definition.
    TABLE_LIMIT = INTERVALS_4_APPROX[NUM_OF_INTERVALS-1];
  }
  

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh, basemesh;
  H2DReader mloader;
  mloader.load(mesh_file, &basemesh);
  
  // Perform initial mesh refinements.
  mesh.copy(&basemesh);
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BDY_TOP, INIT_REF_NUM_BDY_TOP);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_TOP);
  bc_types.add_bc_neumann(Hermes::vector<int>(BDY_RIGHT, BDY_BOTTOM, BDY_LEFT));

  // Enter Dirichlet boundary values.
  BCValues bc_values(&current_time);
  bc_values.add_timedep_function(BDY_TOP, essential_bc_values);

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
      wf.add_matrix_form(jac_form_vol_euler, jac_form_vol_ord, HERMES_NONSYM, HERMES_ANY, 
	                 &sln_prev_time);
      wf.add_vector_form(res_form_vol_euler, res_form_vol_ord, HERMES_ANY, 
			 &sln_prev_time);
    }
    else {
      info("Registering forms for the Newton's method (Crank-Nicolson in time).");
      wf.add_matrix_form(jac_form_vol_cranic, jac_form_vol_ord, HERMES_NONSYM, HERMES_ANY, 
      		         &sln_prev_time);
      wf.add_vector_form(res_form_vol_cranic, res_form_vol_ord, HERMES_ANY, 
			 &sln_prev_time);
    }
  }
  else {
    if (TIME_INTEGRATION == 1) {
      info("Registering forms for the Picard's method (implicit Euler in time).");
      wf.add_matrix_form(bilinear_form_picard_euler, bilinear_form_picard_euler_ord, HERMES_NONSYM, HERMES_ANY, 
	                 &sln_prev_iter);
      wf.add_vector_form(linear_form_picard_euler, linear_form_picard_euler_ord, HERMES_ANY, 
			 Hermes::vector<MeshFunction*>(&sln_prev_iter, &sln_prev_time));
    }
    else {
      info("Registering forms for the Picard's method (Crank-Nicolson in time).");
      error("Not implemented yet.");
      wf.add_matrix_form(bilinear_form_picard_euler, bilinear_form_picard_euler_ord, HERMES_NONSYM, HERMES_ANY, 
	                 &sln_prev_iter);
      wf.add_vector_form(linear_form_picard_euler, linear_form_picard_euler_ord, HERMES_ANY, 
			 Hermes::vector<MeshFunction*>(&sln_prev_iter, &sln_prev_time));
    }
  }

  // Error estimate and discrete problem size as a function of physical time.
  SimpleGraph graph_time_err_est, graph_time_err_exact, 
    graph_time_dof, graph_time_cpu, graph_time_step;
 
  // Visualize the projection and mesh.
  ScalarView view("Initial condition", new WinGeom(0, 0, 630, 350));
  view.fix_scale_width(50);
  OrderView ordview("Initial mesh", new WinGeom(640, 0, 600, 350));
  view.show(&sln_prev_time, HERMES_EPS_VERYHIGH);
  ordview.show(&space);
  //MeshView mview("Mesh", new WinGeom(840, 0, 600, 350));
  //mview.show(&mesh);
  //View::wait();

  // Time stepping loop.
  int num_time_steps = (int)(T_FINAL/time_step + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    info("---- Time step %d:", ts);

    // Time measurement.
    cpu_time.tick();

    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      if (UNREF_LEVEL == 1) mesh.unrefine_all_elements();
      else mesh.copy(&basemesh);
      space.set_uniform_order(P_INIT);
      ndof = Space::get_num_dofs(&space);
    }

    // Spatial adaptivity loop. Note: sln_prev_time must not be touched during adaptivity.
    bool done = false;
    int as = 1;
    double err_est_rel;
    do
    {
      info("---- Time step %d, time step lenght %g, time %g (days), adaptivity step %d:", ts, time_step, current_time, as);

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
        info("Performing Newton's iteration (tau = %g days):", time_step);
        bool success, verbose = true;
        double* save_coeff_vec = new double[ndof];
        // Save coefficient vector.
        for (int i=0; i < ndof; i++) save_coeff_vec[i] = coeff_vec[i];
        double damping_coeff = 1.0;
        while (!solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
                             NEWTON_TOL, NEWTON_MAX_ITER, verbose, damping_coeff)) {
          // Restore solution from the beginning of time step.
          for (int i=0; i < ndof; i++) coeff_vec[i] = save_coeff_vec[i];
          // Reducing time step to 50%.
          info("Reducing time step size from %g to %g days for the rest of this time step.", 
               time_step, time_step * time_step_dec);
          time_step *= time_step_dec;
          // If time_step less than the prescribed minimum, stop.
          if (time_step < time_step_min) error("Time step dropped below prescribed minimum value.");
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
        }

        // Perform Picard iteration on the reference mesh. If necessary, 
        // reduce time step to make it converge, but then restore time step 
        // size to its original value.
        info("Performing Picard's iteration (tau = %g days):", time_step);
        bool success, verbose = true;
        while(!solve_picard(&wf, ref_space, &sln_prev_iter, matrix_solver, PICARD_TOL, 
                            PICARD_MAX_ITER, verbose)) {
          // Restore solution from the beginning of time step.
          sln_prev_iter.copy(&sln_prev_time);
          // Reducing time step to 50%.
          info("Reducing time step size from %g to %g days for the rest of this time step", time_step, time_step * time_step_inc);
          time_step *= time_step_dec;
          // If time_step less than the prescribed minimum, stop.
          if (time_step < time_step_min) error("Time step dropped below prescribed minimum value.");
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
      Adapt* adaptivity = new Adapt(&space);
      
      // Calculate error estimate wrt. fine mesh solution.
      err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

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

    // Add entries to graphs.
    graph_time_err_est.add_values(current_time, err_est_rel);
    graph_time_err_est.save("time_error_est.dat");
    graph_time_dof.add_values(current_time, Space::get_num_dofs(&space));
    graph_time_dof.save("time_dof.dat");
    graph_time_cpu.add_values(current_time, cpu_time.accumulated());
    graph_time_cpu.save("time_cpu.dat");
    graph_time_step.add_values(current_time, time_step);
    graph_time_step.save("time_step_history.dat");

    // Visualize the solution and mesh.
    char title[100];
    sprintf(title, "Solution, time %g days", current_time);
    view.set_title(title);
    view.show(&sln, HERMES_EPS_VERYHIGH);
    sprintf(title, "Mesh, time %g days", current_time);
    ordview.set_title(title);
    ordview.show(&space);
    
    // Save complete Solution.
    char* filename = new char[100];
    sprintf(filename, "outputs/tsln_%f.dat", current_time);
    bool compress = false;   // Gzip compression not used as it only works on Linux.
    sln.save(filename, compress);
    info("Solution at time %g saved to file %s.", current_time, filename);

    // Copy new reference level solution into sln_prev_time.
    // This starts new time step.
    sln_prev_time.copy(&ref_sln);

    // Updating time step. Note that time_step might have been reduced during adaptivity.
    current_time += time_step;

    // Increase time step.
    if (time_step*time_step_inc < time_step_max) {
      info("Increasing time step from %g to %g days.", time_step, time_step * time_step_inc);
      time_step *= time_step_inc;
    }
  }

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
