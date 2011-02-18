#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"

#include "timestep_controller.h"

using namespace RefinementSelectors;

// This test makes sure that the example "nernst-planck-ipmc" works correctly.

// Parameters to tweak the amount of output to the console.
#define NOSCREENSHOT

#define TWO_BASE_MESH

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
const double mech_E = 0.5e9;              //[Pa]
const double mech_nu = 0.487;              // Poisson ratio
const double mech_mu = mech_E / (2 * (1 + mech_nu));
const double mech_lambda = mech_E * mech_nu / ((1 + mech_nu) * (1 - 2 * mech_nu));
const double lin_force_coup = 1e5;


/* Simulation parameters */
const double T_FINAL = 0.1;
double INIT_TAU = 0.05;
double *TAU = &INIT_TAU;                          // Size of the time step
const int P_INIT = 2;       	                  // Initial polynomial degree of all mesh elements.
const int REF_INIT = 3;     	                  // Number of initial refinements.
const bool MULTIMESH = true;	                  // Multimesh?
const int TIME_DISCR = 1;                         // 1 for implicit Euler, 2 for Crank-Nicolson.

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
const double ERR_STOP = 5.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Weak forms
#include "forms.cpp"


/*** Boundary types and conditions ***/

// Boundary markers.
const int BDY_SIDE_FIXED = 1;
const int BDY_SIDE_FREE = 2;
const int BDY_TOP = 3;
const int BDY_BOT = 4;

scalar voltage_ic(double x, double y, double &dx, double &dy) {
  // y^2 function for the domain.
  //return (y+100e-6) * (y+100e-6) / (40000e-12);
  return 0.0;
}

scalar concentration_ic(double x, double y, double &dx, double &dy) {
  return C0;
}

scalar u1_ic(double x, double y, double &dx, double &dy) {
	return 0.0;
}

scalar u2_ic(double x, double y, double &dx, double &dy) {
	return 0.0;
}

int main (int argc, char* argv[]) {

  // Load the mesh file.
  Mesh C_mesh, phi_mesh, u1_mesh, u2_mesh, basemesh;
  H2DReader mloader;
  mloader.load("small.mesh", &basemesh);
#ifdef TWO_BASE_MESH

  Mesh basemesh_electrochem;  // Base mesh to hold C and phi
  Mesh basemesh_deformation;  // Base mesh for deformation displacements u1 and u2

  basemesh_electrochem.copy(&basemesh);
  basemesh_deformation.copy(&basemesh);

  basemesh_electrochem.refine_towards_boundary(BDY_BOT, REF_INIT - 1);
  basemesh_electrochem.refine_all_elements(1);
  C_mesh.copy(&basemesh_electrochem);
  phi_mesh.copy(&basemesh_electrochem);

  for (int i = 0; i < REF_INIT - 1; i++) {
    basemesh_deformation.refine_all_elements(1); //horizontal
    basemesh_deformation.refine_all_elements(2); //vertical
  }

  u1_mesh.copy(&basemesh_deformation);
  u2_mesh.copy(&basemesh_deformation);

#else
  basemesh.refine_towards_boundary(BDY_BOT, REF_INIT);
  basemesh.refine_towards_boundary(BDY_SIDE_FIXED, REF_INIT);
  basemesh.refine_towards_boundary(BDY_SIDE_FREE, REF_INIT - 1);
  basemesh.refine_all_elements(1);
  C_mesh.copy(&basemesh);
  phi_mesh.copy(&basemesh);
  u1_mesh.copy(&basemesh);
  u2_mesh.copy(&basemesh);
#endif

  // Enter Neumann boundary markers for Nernst-Planck.
  BCTypes C_bc_types;
  C_bc_types.add_bc_neumann(Hermes::vector<int>(BDY_SIDE_FIXED, BDY_SIDE_FREE, BDY_TOP, BDY_BOT));

  // Enter Dirichlet and Neumann boundary markers for Poisson.
  BCTypes phi_bc_types;
  phi_bc_types.add_bc_neumann(Hermes::vector<int>(BDY_SIDE_FIXED, BDY_SIDE_FREE));
  phi_bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_TOP, BDY_BOT));

  // Dirichlet and Neumann boundary markers for u1.
  BCTypes u1_bc_types;
  u1_bc_types.add_bc_neumann(Hermes::vector<int>(BDY_SIDE_FREE, BDY_TOP, BDY_BOT));
  u1_bc_types.add_bc_dirichlet(BDY_SIDE_FIXED);

  // Dirichlet and Neumann boundary markers for u1.
  BCTypes u2_bc_types;
  u2_bc_types.add_bc_neumann(Hermes::vector<int>(BDY_SIDE_FREE, BDY_TOP, BDY_BOT));
  u2_bc_types.add_bc_dirichlet(BDY_SIDE_FIXED);

  // Enter Dirichlet boundary values.
  BCValues phi_bc_values;
  phi_bc_values.add_const(BDY_TOP, VOLTAGE);
  phi_bc_values.add_zero(BDY_BOT);

  BCValues u1_bc_values, u2_bc_values;
  u1_bc_values.add_const(BDY_SIDE_FIXED, 0.0);
  u2_bc_values.add_const(BDY_SIDE_FIXED, 0.0);

  BCValues C_bc_values;

  // Spaces for concentration and the voltage.
  H1Space C_space(&C_mesh, &C_bc_types, &C_bc_values, P_INIT);
  H1Space phi_space(MULTIMESH ? &phi_mesh : &C_mesh, &phi_bc_types, &phi_bc_values, P_INIT);
  H1Space u1_space(MULTIMESH ? &u1_mesh : &C_mesh, &u1_bc_types, &u1_bc_values, P_INIT);
  H1Space u2_space(MULTIMESH ? &u2_mesh : &C_mesh, &u2_bc_types, &u2_bc_values, P_INIT);

  int ndof = Space::get_num_dofs(Hermes::vector<Space*>(&C_space, &phi_space, &u1_space, &u2_space));

  Solution C_sln, C_ref_sln;
  Solution phi_sln, phi_ref_sln; 
  Solution u1_sln, u1_ref_sln;
  Solution u2_sln, u2_ref_sln;

  // Assign initial condition to mesh.
  Solution C_prev_time(&C_mesh, concentration_ic);
  Solution phi_prev_time(MULTIMESH ? &phi_mesh : &C_mesh, voltage_ic);
  Solution u1_prev_time(MULTIMESH ? &u1_mesh : &C_mesh, u1_ic);
  Solution u2_prev_time(MULTIMESH ? &u2_mesh : &C_mesh, u2_ic);

  // The weak form for 2 equations.
  WeakForm wf(4);
  // Add the bilinear and linear forms.
  if (TIME_DISCR == 1) {  // Implicit Euler.
    // the first row in the 4 x 4 matrix
    wf.add_matrix_form(0, 0, callback(J_euler_DFcDYc), HERMES_NONSYM);
    wf.add_matrix_form(0, 1, callback(J_euler_DFcDYphi), HERMES_NONSYM);
    wf.add_matrix_form(0, 2, callback(J_euler_DFcDYu1), HERMES_NONSYM);
    wf.add_matrix_form(0, 3, callback(J_euler_DFcDYu2), HERMES_NONSYM);
    // the second row
    wf.add_matrix_form(1, 0, callback(J_euler_DFphiDYc), HERMES_NONSYM);
    wf.add_matrix_form(1, 1, callback(J_euler_DFphiDYphi), HERMES_NONSYM);
    wf.add_matrix_form(1, 2, callback(J_euler_DFphiDYu1), HERMES_NONSYM);
    wf.add_matrix_form(1, 3, callback(J_euler_DFphiDYu2), HERMES_NONSYM);
    // the third row
    wf.add_matrix_form(2, 0, callback(J_euler_DFu1DYc), HERMES_NONSYM);
    wf.add_matrix_form(2, 1, callback(J_euler_DFu1DYphi), HERMES_NONSYM);
    wf.add_matrix_form(2, 2, callback(J_euler_DFu1DYu1), HERMES_NONSYM);
    wf.add_matrix_form(2, 3, callback(J_euler_DFu1DYu2), HERMES_NONSYM);
    // the fourth row
    wf.add_matrix_form(3, 0, callback(J_euler_DFu2DYc), HERMES_NONSYM);
    wf.add_matrix_form(3, 1, callback(J_euler_DFu2DYphi), HERMES_NONSYM);
    wf.add_matrix_form(3, 2, callback(J_euler_DFu2DYu1), HERMES_NONSYM);
    wf.add_matrix_form(3, 3, callback(J_euler_DFu2DYu2), HERMES_NONSYM);
    // the vector forms

    wf.add_vector_form(0, callback(Fc_euler), HERMES_ANY,
                       Hermes::vector<MeshFunction*>(&C_prev_time, &phi_prev_time));
    wf.add_vector_form(1, callback(Fphi_euler), HERMES_ANY,
                       Hermes::vector<MeshFunction*>(&C_prev_time, &phi_prev_time));
    wf.add_vector_form(2, callback(Fu1_euler), HERMES_ANY);
    wf.add_vector_form(3, callback(Fu2_euler), HERMES_ANY);
  } else {
	  error("Crank-Nicholson forms are not implemented yet");
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
  OGProjection::project_global(Hermes::vector<Space *>(&C_space, &phi_space, &u1_space, &u2_space),
                               Hermes::vector<MeshFunction *>(&C_prev_time, &phi_prev_time, &u1_prev_time, &u2_prev_time),
                               coeff_vec_coarse, matrix_solver);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp_coarse(&wf, Hermes::vector<Space *>(&C_space, &phi_space, &u1_space, &u2_space), is_linear);

  // Set up the solver, matrix, and rhs for the coarse mesh according to the solver selection.
  SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
  Vector* rhs_coarse = create_vector(matrix_solver);
  Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh:");
  bool verbose = true;
  if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, 
      NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

  // Translate the resulting coefficient vector into the Solution sln.
  Solution::vector_to_solutions(coeff_vec_coarse, Hermes::vector<Space *>(&C_space, &phi_space, &u1_space, &u2_space),
                                Hermes::vector<Solution *>(&C_sln, &phi_sln, &u1_sln, &u2_sln));

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

#ifdef TWO_BASE_MESH
      C_mesh.copy(&basemesh_electrochem);
#else
      C_mesh.copy(&basemesh);
#endif
      if (MULTIMESH)
      {
#ifdef TWO_BASE_MESH
        phi_mesh.copy(&basemesh_electrochem);
        u1_mesh.copy(&basemesh_deformation);
        u2_mesh.copy(&basemesh_deformation);
#else
        phi_mesh.copy(&basemesh);
        u1_mesh.copy(&basemesh);
        u2_mesh.copy(&basemesh);
#endif

      }
      C_space.set_uniform_order(P_INIT);
      phi_space.set_uniform_order(P_INIT);
      u1_space.set_uniform_order(P_INIT);
      u2_space.set_uniform_order(P_INIT);

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
      Hermes::vector<Space *>* ref_spaces = construct_refined_spaces(
          Hermes::vector<Space *>(&C_space, &phi_space, &u1_space, &u2_space));

      scalar* coeff_vec = new scalar[Space::get_num_dofs(*ref_spaces)];
      DiscreteProblem* dp = new DiscreteProblem(&wf, *ref_spaces, is_linear);
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Calculate initial coefficient vector for Newton on the fine mesh.
      if (as == 1 && pid.get_timestep_number() == 1) {
        info("Projecting coarse mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces, Hermes::vector<MeshFunction *>(&C_sln, &phi_sln, &u1_sln, &u2_sln),
                                     coeff_vec, matrix_solver);
      }
      else {
        info("Projecting previous fine mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces,
            Hermes::vector<MeshFunction *>(&C_ref_sln, &phi_ref_sln, &u1_ref_sln, &u2_ref_sln),
            coeff_vec, matrix_solver);
      }
      if (as > 1) {
        // Now deallocate the previous mesh
        info("Delallocating the previous mesh");
        delete C_ref_sln.get_mesh();
        delete phi_ref_sln.get_mesh();
        delete u1_ref_sln.get_mesh();
        delete u2_ref_sln.get_mesh();
      }

      // Newton's loop on the fine mesh.
      info("Solving on fine mesh:");
      if (!solve_newton(coeff_vec, dp, solver, matrix, rhs, 
	  	      NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");


      // Store the result in ref_sln.
      Solution::vector_to_solutions(coeff_vec, *ref_spaces, 
                                    Hermes::vector<Solution *>(&C_ref_sln, &phi_ref_sln, &u1_ref_sln, &u2_ref_sln));
      // Projecting reference solution onto the coarse mesh
      info("Projecting fine mesh solution on coarse mesh.");
      OGProjection::project_global(Hermes::vector<Space *>(&C_space, &phi_space, &u1_space, &u2_space),
                                   Hermes::vector<Solution *>(&C_ref_sln, &phi_ref_sln, &u1_ref_sln, &u2_ref_sln),
                                   Hermes::vector<Solution *>(&C_sln, &phi_sln, &u1_sln, &u2_sln),
                                   matrix_solver);

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&C_space, &phi_space, &u1_space, &u2_space));
      Hermes::vector<double> err_est_rel;
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&C_sln, &phi_sln, &u1_sln, &u2_sln),
                                 Hermes::vector<Solution *>(&C_ref_sln, &phi_ref_sln, &u1_ref_sln, &u2_ref_sln),
                                 &err_est_rel) * 100;

      // Report results.
      info("ndof_coarse[0]: %d, ndof_fine[0]: %d",
           C_space.get_num_dofs(), (*ref_spaces)[0]->get_num_dofs());
      info("err_est_rel[0]: %g%%", err_est_rel[0]*100);
      info("ndof_coarse[1]: %d, ndof_fine[1]: %d",
           phi_space.get_num_dofs(), (*ref_spaces)[1]->get_num_dofs());
      info("err_est_rel[1]: %g%%", err_est_rel[1]*100);
      info("ndof_coarse[2]: %d, ndof_fine[2]: %d",
           u1_space.get_num_dofs(), (*ref_spaces)[2]->get_num_dofs());
      info("err_est_rel[2]: %g%%", err_est_rel[3]*100);
      info("ndof_coarse[3]: %d, ndof_fine[3]: %d",
            u2_space.get_num_dofs(), (*ref_spaces)[3]->get_num_dofs());
      info("err_est_rel[3]: %g%%", err_est_rel[3]*100);

      // Report results.
      info("ndof_coarse_total: %d, ndof_fine_total: %d, err_est_rel: %g%%", 
           Space::get_num_dofs(Hermes::vector<Space *>(&C_space, &phi_space, &u1_space, &u2_space)),
                               Space::get_num_dofs(*ref_spaces), err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
      else 
      {
        info("Adapting the coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&selector, &selector, &selector, &selector),
          THRESHOLD, STRATEGY, MESH_REGULARITY);
        
        info("Adapted...");

        if (Space::get_num_dofs(Hermes::vector<Space *>(&C_space, &phi_space, &u1_space, &u2_space)) >= NDOF_STOP)
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }

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


    pid.end_step(Hermes::vector<Solution*> (&C_ref_sln, &phi_ref_sln, &u1_ref_sln, &u2_ref_sln),
        Hermes::vector<Solution*> (&C_prev_time, &phi_prev_time, &u1_prev_time, &u2_prev_time));
    // TODO! Time step reduction when necessary.

    // Copy last reference solution into sln_prev_time.
    C_prev_time.copy(&C_ref_sln);
    phi_prev_time.copy(&phi_ref_sln);
    u1_prev_time.copy(&u1_ref_sln);
    u2_prev_time.copy(&u2_ref_sln);

  } while (pid.has_next());

  ndof = Space::get_num_dofs(Hermes::vector<Space*>(&C_space, &phi_space, &u1_space, &u2_space));

  printf("ndof allowed = %d\n", 300);
  printf("ndof actual = %d\n", ndof);
  if (ndof < 300) {      // ndofs was 258 at the time this test was created
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

