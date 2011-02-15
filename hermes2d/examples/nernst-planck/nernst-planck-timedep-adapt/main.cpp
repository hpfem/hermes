#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"

#include "timestep_controller.h"

using namespace RefinementSelectors;

/** \addtogroup e_newton_np_timedep_adapt_system Newton Time-dependant System with Adaptivity
 \{
 \brief This example shows how to combine the automatic adaptivity with the Newton's method for a nonlinear time-dependent PDE system.

 This example shows how to combine the automatic adaptivity with the
 Newton's method for a nonlinear time-dependent PDE system.
 The time discretization is done using implicit Euler or
 Crank Nicholson method (see parameter TIME_DISCR).
 The following PDE's are solved:
 Nernst-Planck (describes the diffusion and migration of charged particles):
 \f[dC/dt - D*div[grad(C)] - K*C*div[grad(\phi)]=0,\f]
 where D and K are constants and C is the cation concentration variable,
 phi is the voltage variable in the Poisson equation:
 \f[ - div[grad(\phi)] = L*(C - C_0),\f]
 where \f$C_0\f$, and L are constant (anion concentration). \f$C_0\f$ is constant
 anion concentration in the domain and L is material parameter.
 So, the equation variables are phi and C and the system describes the
 migration/diffusion of charged particles due to applied voltage.
 The simulation domain looks as follows:
 \verbatim
      2
  +----------+
  |          |
 1|          |1
  |          |
  +----------+
      3
 \endverbatim
 For the Nernst-Planck equation, all the boundaries are natural i.e. Neumann.
 Which basically means that the normal derivative is 0:
 \f[ BC: -D*dC/dn - K*C*d\phi/dn = 0 \f]
 For Poisson equation, boundary 1 has a natural boundary condition
 (electric field derivative is 0).
 The voltage is applied to the boundaries 2 and 3 (Dirichlet boundaries)
 It is possible to adjust system paramter VOLT_BOUNDARY to apply
 Neumann boundary condition to 2 (instead of Dirichlet). But by default:
  - BC 2: \f$\phi = VOLTAGE\f$
  - BC 3: \f$\phi = 0\f$
  - BC 1: \f$\frac{d\phi}{dn} = 0\f$
 */

// Parameters to tweak the amount of output to the console.
#define NOSCREENSHOT

/*** Fundamental coefficients ***/
const double D = 10e-11; 	                  // [m^2/s] Diffusion coefficient.
const double R = 8.31; 		                  // [J/mol*K] Gas constant.
const double T = 293; 		                  // [K] Aboslute temperature.
const double F = 96485.3415;	                  // [s * A / mol] Faraday constant.
const double eps = 2.5e-2; 	                  // [F/m] Electric permeability.
const double mu = D / (R * T);                    // Mobility of ions.
const double z = 1;		                  // Charge number.
const double K = z * mu * F;                      // Constant for equation.
const double L =  F / eps;	                  // Constant for equation.
const double VOLTAGE = 1;	                  // [V] Applied voltage.
const scalar C0 = 1200;	                          // [mol/m^3] Anion and counterion concentration.


/* Simulation parameters */
const double T_FINAL = 1;
double INIT_TAU = 0.05;
double *TAU = &INIT_TAU;                          // Size of the time step
const int P_INIT = 2;       	                  // Initial polynomial degree of all mesh elements.
const int REF_INIT = 3;     	                  // Number of initial refinements.
const bool MULTIMESH = true;	                  // Multimesh?
const int TIME_DISCR = 2;                         // 1 for implicit Euler, 2 for Crank-Nicolson.

const double NEWTON_TOL_COARSE = 0.01;            // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;              // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.

const int UNREF_FREQ = 1;                         // every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
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
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int NDOF_STOP = 5000;	                  // To prevent adaptivity from going on forever.
const double ERR_STOP = 0.1;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Weak forms
#include "forms.cpp"


/*** Boundary types and conditions ***/

// Boundary markers.
const int BDY_SIDE = 1;
const int BDY_TOP = 2;
const int BDY_BOT = 3;

scalar voltage_ic(double x, double y, double &dx, double &dy) {
  // y^2 function for the domain.
  //return (y+100e-6) * (y+100e-6) / (40000e-12);
  return 0.0;
}

scalar concentration_ic(double x, double y, double &dx, double &dy) {
  return C0;
}

int main (int argc, char* argv[]) {

  // Load the mesh file.
  Mesh C_mesh, phi_mesh, basemesh;
  H2DReader mloader;
  mloader.load("small.mesh", &basemesh);
  
  // When nonadaptive solution, refine the mesh.
  basemesh.refine_towards_boundary(BDY_TOP, REF_INIT);
  basemesh.refine_towards_boundary(BDY_BOT, REF_INIT - 1);
  basemesh.refine_all_elements(1);
  basemesh.refine_all_elements(1);
  C_mesh.copy(&basemesh);
  phi_mesh.copy(&basemesh);

  // Enter Neumann boundary markers for Nernst-Planck.
  BCTypes C_bc_types;
  C_bc_types.add_bc_neumann(Hermes::vector<int>(BDY_SIDE, BDY_TOP, BDY_BOT));

  // Enter Dirichlet and Neumann boundary markers for Poisson.
  BCTypes phi_bc_types;
  phi_bc_types.add_bc_neumann(BDY_SIDE);
  phi_bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_TOP, BDY_BOT));

  // Enter Dirichlet boundary values.
  BCValues phi_bc_values;
  phi_bc_values.add_const(BDY_TOP, VOLTAGE);
  phi_bc_values.add_zero(BDY_BOT);

  BCValues C_bc_values;
//  C_bc_values.add_zero(Hermes::vector<int>(BDY_SIDE, BDY_TOP, BDY_BOT));

  // Spaces for concentration and the voltage.
  H1Space C_space(&C_mesh, &C_bc_types, &C_bc_values, P_INIT);
  H1Space phi_space(MULTIMESH ? &phi_mesh : &C_mesh, &phi_bc_types, &phi_bc_values, P_INIT);
  int ndof = Space::get_num_dofs(Hermes::vector<Space*>(&C_space, &phi_space));

  Solution C_sln, C_ref_sln;
  Solution phi_sln, phi_ref_sln; 

  // Assign initial condition to mesh.
  Solution C_prev_time(&C_mesh, concentration_ic);
  Solution phi_prev_time(MULTIMESH ? &phi_mesh : &C_mesh, voltage_ic);

  // The weak form for 2 equations.
  WeakForm wf(2);
  // Add the bilinear and linear forms.
  if (TIME_DISCR == 1) {  // Implicit Euler.
  wf.add_matrix_form(0, 0, callback(J_euler_DFcDYc), HERMES_NONSYM);
  wf.add_matrix_form(0, 1, callback(J_euler_DFcDYphi), HERMES_NONSYM);
  wf.add_matrix_form(1, 0, callback(J_euler_DFphiDYc), HERMES_NONSYM);
  wf.add_matrix_form(1, 1, callback(J_euler_DFphiDYphi), HERMES_NONSYM);
  wf.add_vector_form(0, callback(Fc_euler), HERMES_ANY, 
                     Hermes::vector<MeshFunction*>(&C_prev_time, &phi_prev_time));
  wf.add_vector_form(1, callback(Fphi_euler), HERMES_ANY, 
                     Hermes::vector<MeshFunction*>(&C_prev_time, &phi_prev_time));
  } else {
    wf.add_matrix_form(0, 0, callback(J_cranic_DFcDYc), HERMES_NONSYM);
    wf.add_matrix_form(0, 1, callback(J_cranic_DFcDYphi), HERMES_NONSYM);
    wf.add_matrix_form(1, 0, callback(J_cranic_DFphiDYc), HERMES_NONSYM);
    wf.add_matrix_form(1, 1, callback(J_cranic_DFphiDYphi), HERMES_NONSYM);
    wf.add_vector_form(0, callback(Fc_cranic), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&C_prev_time, &phi_prev_time));
    wf.add_vector_form(1, callback(Fphi_cranic), HERMES_ANY);
  }

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  scalar* coeff_vec_coarse = new scalar[ndof];
  OGProjection::project_global(Hermes::vector<Space *>(&C_space, &phi_space), 
                               Hermes::vector<MeshFunction *>(&C_prev_time, &phi_prev_time), 
                               coeff_vec_coarse, matrix_solver);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp_coarse(&wf, Hermes::vector<Space *>(&C_space, &phi_space), is_linear);

  // Set up the solver, matrix, and rhs for the coarse mesh according to the solver selection.
  SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
  Vector* rhs_coarse = create_vector(matrix_solver);
  Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Visualization windows.
  char title[1000];
  ScalarView Cview("Concentration [mol/m3]", new WinGeom(0, 0, 800, 800));
  ScalarView phiview("Voltage [V]", new WinGeom(650, 0, 600, 600));
  OrderView Cordview("C order", new WinGeom(0, 300, 600, 600));
  OrderView phiordview("Phi order", new WinGeom(600, 300, 600, 600));

  Cview.show(&C_prev_time);
  Cordview.show(&C_space);
  phiview.show(&phi_prev_time);
  phiordview.show(&phi_space);

  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh:");
  bool verbose = true;
  if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, 
      NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution::vector_to_solutions(coeff_vec_coarse, Hermes::vector<Space *>(&C_space, &phi_space), 
                                Hermes::vector<Solution *>(&C_sln, &phi_sln));

  Cview.show(&C_sln);
  phiview.show(&phi_sln);

  // Cleanup after the Newton loop on the coarse mesh.
  delete matrix_coarse;
  delete rhs_coarse;
  delete solver_coarse;
  delete[] coeff_vec_coarse;
  
  // Time stepping loop.
  PidTimestepController pid(T_FINAL, true, INIT_TAU);
  TAU = pid.timestep;
  info("Starting time iteration with the step %g", *TAU);

  do {
    pid.begin_step();
    // Periodic global derefinements.
    if (pid.get_timestep_number() > 1 && pid.get_timestep_number() % UNREF_FREQ == 0)
    {
      info("Global mesh derefinement.");
      C_mesh.copy(&basemesh);
      if (MULTIMESH)
      {
        phi_mesh.copy(&basemesh);
      }
      C_space.set_uniform_order(P_INIT);
      phi_space.set_uniform_order(P_INIT);

      // Project on globally derefined mesh.
      //info("Projecting previous fine mesh solution on derefined mesh.");
      //OGProjection::project_global(Hermes::vector<Space *>(&C, &phi), Hermes::vector<Solution *>(&C_ref_sln, &phi_ref_sln), 
       //                            Hermes::vector<Solution *>(&C_sln, &phi_sln));
    }

    // Adaptivity loop. Note: C_prev_time and Phi_prev_time must not be changed during spatial adaptivity.
    bool done = false; int as = 1;
    double err_est;
    do {
      info("Time step %d, adaptivity step %d:", pid.get_timestep_number(), as);

      // Construct globally refined reference mesh
      // and setup reference space.
      Hermes::vector<Space *>* ref_spaces = construct_refined_spaces(Hermes::vector<Space *>(&C_space, &phi_space));

      scalar* coeff_vec = new scalar[Space::get_num_dofs(*ref_spaces)];
      DiscreteProblem* dp = new DiscreteProblem(&wf, *ref_spaces, is_linear);
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Calculate initial coefficient vector for Newton on the fine mesh.
      if (as == 1 && pid.get_timestep_number() == 1) {
        info("Projecting coarse mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces, Hermes::vector<MeshFunction *>(&C_sln, &phi_sln), 
                                     coeff_vec, matrix_solver);
      }
      else {
        info("Projecting previous fine mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces, Hermes::vector<MeshFunction *>(&C_ref_sln, &phi_ref_sln), 
                                     coeff_vec, matrix_solver);
      }
      if (as > 1) {
        // Now deallocate the previous mesh
        info("Delallocating the previous mesh");
        delete C_ref_sln.get_mesh();
        delete phi_ref_sln.get_mesh();
      }

      // Newton's loop on the fine mesh.
      info("Solving on fine mesh:");
      if (!solve_newton(coeff_vec, dp, solver, matrix, rhs, 
	  	      NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");


      // Store the result in ref_sln.
      Solution::vector_to_solutions(coeff_vec, *ref_spaces, 
                                    Hermes::vector<Solution *>(&C_ref_sln, &phi_ref_sln));
      // Projecting reference solution onto the coarse mesh
      info("Projecting fine mesh solution on coarse mesh.");
      OGProjection::project_global(Hermes::vector<Space *>(&C_space, &phi_space), 
                                   Hermes::vector<Solution *>(&C_ref_sln, &phi_ref_sln), 
                                   Hermes::vector<Solution *>(&C_sln, &phi_sln),
                                   matrix_solver);

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&C_space, &phi_space));
      Hermes::vector<double> err_est_rel;
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&C_sln, &phi_sln), 
                                 Hermes::vector<Solution *>(&C_ref_sln, &phi_ref_sln), &err_est_rel) * 100;

      // Report results.
      info("ndof_coarse[0]: %d, ndof_fine[0]: %d",
           C_space.get_num_dofs(), (*ref_spaces)[0]->get_num_dofs());
      info("err_est_rel[0]: %g%%", err_est_rel[0]*100);
      info("ndof_coarse[1]: %d, ndof_fine[1]: %d",
           phi_space.get_num_dofs(), (*ref_spaces)[1]->get_num_dofs());
      info("err_est_rel[1]: %g%%", err_est_rel[1]*100);
      // Report results.
      info("ndof_coarse_total: %d, ndof_fine_total: %d, err_est_rel: %g%%", 
           Space::get_num_dofs(Hermes::vector<Space *>(&C_space, &phi_space)), 
                               Space::get_num_dofs(*ref_spaces), err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
      else 
      {
        info("Adapting the coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&selector, &selector),
          THRESHOLD, STRATEGY, MESH_REGULARITY);
        
        info("Adapted...");

        if (Space::get_num_dofs(Hermes::vector<Space *>(&C_space, &phi_space)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }

      // Visualize the solution and mesh.
      info("Visualization procedures: C");
      char title[100];
      sprintf(title, "Solution[C], time step# %d, step size %g, time %g",
          pid.get_timestep_number(), *TAU, pid.get_time());
      Cview.set_title(title);
      Cview.show(&C_ref_sln);
      sprintf(title, "Mesh[C], time step# %d, step size %g, time %g",
          pid.get_timestep_number(), *TAU, pid.get_time());
      Cordview.set_title(title);
      Cordview.show(&C_space);
      
      info("Visualization procedures: phi");
      sprintf(title, "Solution[phi], time step# %d, step size %g, time %g",
          pid.get_timestep_number(), *TAU, pid.get_time());
      phiview.set_title(title);
      phiview.show(&phi_ref_sln);
      sprintf(title, "Mesh[phi], time step# %d, step size %g, time %g",
          pid.get_timestep_number(), *TAU, pid.get_time());
      phiordview.set_title(title);
      phiordview.show(&phi_space);
      //View::wait(HERMES_WAIT_KEYPRESS);

      // Clean up.
      info("delete solver");
      delete solver;
      info("delete matrix");
      delete matrix;
      info("delete rhs");
      delete rhs;
      info("delete adaptivity");
      delete adaptivity;
      info("delete[] ref_spaces");
      delete ref_spaces;
      info("delete dp");
      delete dp;
      info("delete[] coeff_vec");
      delete[] coeff_vec;
    }
    while (done == false);


    pid.end_step(Hermes::vector<Solution*> (&C_ref_sln, &phi_ref_sln), Hermes::vector<Solution*> (&C_prev_time, &phi_prev_time));
    // TODO! Time step reduction when necessary.

    // Copy last reference solution into sln_prev_time.
    C_prev_time.copy(&C_ref_sln);
    phi_prev_time.copy(&phi_ref_sln);

  } while (pid.has_next());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

