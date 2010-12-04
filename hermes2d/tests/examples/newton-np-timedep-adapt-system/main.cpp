#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"

#include "hermes2d.h"
#include <string>

using namespace RefinementSelectors;

// This test makes sure that the example "newton-np-timedep-adapt-system" works correctly.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

#define SIDE_MARKER 1
#define TOP_MARKER 2
#define BOT_MARKER 3

// Parameters to tweak the amount of output to the console.
#define NOSCREENSHOT

/*** Fundamental coefficients ***/
const double D = 10e-11; 	            // [m^2/s] Diffusion coefficient.
const double R = 8.31; 		            // [J/mol*K] Gas constant.
const double T = 293; 		            // [K] Aboslute temperature.
const double F = 96485.3415;	            // [s * A / mol] Faraday constant.
const double eps = 2.5e-2; 	            // [F/m] Electric permeability.
const double mu = D / (R * T);              // Mobility of ions.
const double z = 1;		            // Charge number.
const double K = z * mu * F;                // Constant for equation.
const double L =  F / eps;	            // Constant for equation.
const double VOLTAGE = 1;	            // [V] Applied voltage.
const scalar C0 = 1200;	                    // [mol/m^3] Anion and counterion concentration.

/* For Neumann boundary */
const double height = 180e-6;	            // [m] thickness of the domain.
const double E_FIELD = VOLTAGE / height;    // Boundary condtion for positive voltage electrode


/* Simulation parameters */
const int NSTEP = 50;                 // Number of time steps.
const double TAU = 0.1;               // Size of the time step.
const int P_INIT = 3;       	      // Initial polynomial degree of all mesh elements.
const int REF_INIT = 1;     	      // Number of initial refinements.
const bool MULTIMESH = false;	      // Multimesh?
const int TIME_DISCR = 2;             // 1 for implicit Euler, 2 for Crank-Nicolson.
const int VOLT_BOUNDARY = 1;          // 1 for Dirichlet, 2 for Neumann.

/* Nonadaptive solution parameters */
const double NEWTON_TOL = 1e-6;       // Stopping criterion for nonadaptive solution.

/* Adaptive solution parameters */
const bool SOLVE_ON_COARSE_MESH = false;  // true... Newton is done on coarse mesh in every adaptivity step.
                                      // false...Newton is done on coarse mesh only once, then projection 
                                      // of the fine mesh solution to coarse mesh is used.
const double NEWTON_TOL_COARSE = 0.01;// Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;  // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 100;      // Maximum allowed number of Newton iterations.

const int UNREF_FREQ = 5;             // every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;         // This is a quantitative parameter of the adapt(...) function and
                                      // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;               // Adaptive strategy:
                                      // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                      //   error is processed. If more elements have similar errors, refine
                                      //   all to keep the mesh symmetric.
                                      // STRATEGY = 1 ... refine all elements whose error is larger
                                      //   than THRESHOLD times maximum element error.
                                      // STRATEGY = 2 ... refine all elements whose error is larger
                                      //   than THRESHOLD.
                                      // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;  // Maximum allowed level of hanging nodes:
                                 // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                 // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                 // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                 // Note that regular meshes are not supported, this is due to
                                 // their notoriously bad performance.
const double CONV_EXP = 1.0;     // Default value is 1.0. This parameter influences the selection of
                                 // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int NDOF_STOP = 5000;		        // To prevent adaptivity from going on forever.
const double ERR_STOP = 0.1;          // Stopping criterion for adaptivity (rel. error tolerance between the
                                      // fine mesh and coarse mesh solution in percent).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Program parameters
const std::string USE_ADAPTIVE("adapt");

// Weak forms
#include "forms.cpp"


/*** Boundary types and conditions ***/

// Poisson takes Dirichlet and Neumann boundaries
BCType phi_bc_types(int marker) {
  return (marker == SIDE_MARKER || (marker == TOP_MARKER && VOLT_BOUNDARY == 2)) 
    ? BC_NATURAL : BC_ESSENTIAL;
}

// Nernst-Planck takes Neumann boundaries
BCType C_bc_types(int marker) {
  return BC_NATURAL;
}

// Diricleht Boundary conditions for Poisson equation.
scalar C_essential_bc_values(int ess_bdy_marker, double x, double y) {
  return 0;
}

// Diricleht Boundary conditions for Poisson equation.
scalar phi_essential_bc_values(int ess_bdy_marker, double x, double y) {
  return ess_bdy_marker == TOP_MARKER ? VOLTAGE : 0.0;
}

template<class Real, class Scalar>
Scalar linear_form_surf_top(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
  return -E_FIELD * int_v<Real, Scalar>(n, wt, v);
}

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
  Mesh Cmesh, phimesh, basemesh;
  H2DReader mloader;
  mloader.load("small.mesh", &basemesh);
  
  basemesh.refine_towards_boundary(TOP_MARKER, REF_INIT);
  basemesh.refine_towards_boundary(BOT_MARKER, REF_INIT - 1);
  Cmesh.copy(&basemesh);
  phimesh.copy(&basemesh);

  // Spaces for concentration and the voltage.
  H1Space C(&Cmesh, C_bc_types, C_essential_bc_values, P_INIT);
  H1Space phi(MULTIMESH ? &phimesh : &Cmesh, phi_bc_types, phi_essential_bc_values, P_INIT);
  int ndof = Space::get_num_dofs(Hermes::Tuple<Space*>(&C, &phi));

  Solution C_sln, C_ref_sln;
  Solution phi_sln, phi_ref_sln; 

  // Assign initial condition to mesh.
  Solution C_prev_time(&Cmesh, concentration_ic);
  Solution phi_prev_time(MULTIMESH ? &phimesh : &Cmesh, voltage_ic);

  // The weak form for 2 equations.
  WeakForm wf(2);
  // Add the bilinear and linear forms.
  wf.add_matrix_form(0, 0, callback(J_euler_DFcDYc), HERMES_UNSYM, HERMES_ANY, &phi_prev_time);
  wf.add_matrix_form(0, 1, callback(J_euler_DFcDYphi), HERMES_UNSYM, HERMES_ANY, &C_prev_time);
  wf.add_matrix_form(1, 0, callback(J_euler_DFphiDYc), HERMES_UNSYM);
  wf.add_matrix_form(1, 1, callback(J_euler_DFphiDYphi), HERMES_UNSYM);
  wf.add_vector_form(0, callback(Fc_euler), HERMES_ANY, Hermes::Tuple<MeshFunction*>(&C_prev_time, &phi_prev_time));
  wf.add_vector_form(1, callback(Fphi_euler), HERMES_ANY, Hermes::Tuple<MeshFunction*>(&C_prev_time, &phi_prev_time));

  // Neumann voltage boundary.
  wf.add_vector_form_surf(1, callback(linear_form_surf_top), TOP_MARKER);
 
  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  scalar* coeff_vec_coarse = new scalar[ndof];
  OGProjection::project_global(Hermes::Tuple<Space *>(&C, &phi), Hermes::Tuple<MeshFunction *>(&C_prev_time, &phi_prev_time), coeff_vec_coarse, matrix_solver);

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp_coarse(&wf, Hermes::Tuple<Space *>(&C, &phi), is_linear);

  // Set up the solver, matrix, and rhs for the coarse mesh according to the solver selection.
  SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
  Vector* rhs_coarse = create_vector(matrix_solver);
  Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh:");
  int it = 1;
  while (1)
  {
    // Obtain the number of degrees of freedom.
    int ndof = Space::get_num_dofs(Hermes::Tuple<Space *>(&C, &phi));

    // Assemble the Jacobian matrix and residual vector.
    dp_coarse.assemble(coeff_vec_coarse, matrix_coarse, rhs_coarse, false);

    // Multiply the residual vector with -1 since the matrix 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    for (int i = 0; i < ndof; i++) rhs_coarse->set(i, -rhs_coarse->get(i));
    
    // Calculate the l2-norm of residual vector.
    double res_l2_norm = get_l2_norm(rhs_coarse);

    // Info for user.
    info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(Hermes::Tuple<Space *>(&C, &phi)), res_l2_norm);

    // If l2 norm of the residual vector is in tolerance, or the maximum number 
    // of iteration has been hit, then quit.
    if (res_l2_norm < NEWTON_TOL_COARSE || it > NEWTON_MAX_ITER) break;

    // Solve the linear system and if successful, obtain the solution.
    if(!solver_coarse->solve())
      error ("Matrix solver failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) coeff_vec_coarse[i] += solver_coarse->get_solution()[i];
    
    if (it >= NEWTON_MAX_ITER)
      error ("Newton method did not converge.");

    it++;
  }

  // Translate the resulting coefficient vector into the Solution sln.
  Solution::vector_to_solutions(coeff_vec_coarse, Hermes::Tuple<Space *>(&C, &phi), Hermes::Tuple<Solution *>(&C_sln, &phi_sln));

  // Cleanup after the Newton loop on the coarse mesh.
  delete matrix_coarse;
  delete rhs_coarse;
  delete solver_coarse;
  delete [] coeff_vec_coarse;
  
  // Time stepping loop.
  int num_time_steps = (int)(NEWTON_MAX_ITER/TAU + 0.5);
  for(int ts = 1; ts <= num_time_steps; ts++)
  {
    // Periodic global derefinements.
    if (ts > 1 && ts % UNREF_FREQ == 0) 
    {
      info("Global mesh derefinement.");
      Cmesh.copy(&basemesh);
      if (MULTIMESH)
      {
        phimesh.copy(&basemesh);
      }
      C.set_uniform_order(P_INIT);
      phi.set_uniform_order(P_INIT);

      // Project on globally derefined mesh.
      info("Projecting previous fine mesh solution on derefined mesh.");
      OGProjection::project_global(Hermes::Tuple<Space *>(&C, &phi), Hermes::Tuple<Solution *>(&C_ref_sln, &phi_ref_sln), Hermes::Tuple<Solution *>(&C_sln, &phi_sln));
    }

    // Adaptivity loop:
    bool done = false; int as = 1;
    double err_est;
    do {
      info("Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh
      // and setup reference space.
      Hermes::Tuple<Space *>* ref_spaces = construct_refined_spaces(Hermes::Tuple<Space *>(&C, &phi));

      scalar* coeff_vec = new scalar[Space::get_num_dofs(*ref_spaces)];
      DiscreteProblem* dp = new DiscreteProblem(&wf, *ref_spaces, is_linear);
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Calculate initial coefficient vector for Newton on the fine mesh.
      if (as == 1) {
        info("Projecting coarse mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces, Hermes::Tuple<MeshFunction *>(&C_sln, &phi_sln), coeff_vec);
      }
      else {
        info("Projecting previous fine mesh solution to obtain coefficient vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces, Hermes::Tuple<MeshFunction *>(&C_ref_sln, &phi_ref_sln), coeff_vec);
      }

      // Newton's loop on the fine mesh.
      info("Solving on fine mesh:");
      int it = 1;
      while (1)
      {
        // Obtain the number of degrees of freedom.
        int ndof = Space::get_num_dofs(*ref_spaces);

        // Assemble the Jacobian matrix and residual vector.
        dp->assemble(coeff_vec, matrix, rhs, false);

        // Multiply the residual vector with -1 since the matrix 
        // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
        for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));
        
        // Calculate the l2-norm of residual vector.
        double res_l2_norm = get_l2_norm(rhs);

        // Info for user.
        info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(*ref_spaces), res_l2_norm);

        // If l2 norm of the residual vector is within tolerance, or the maximum number 
        // of iteration has been reached, then quit.
        if (res_l2_norm < NEWTON_TOL_FINE || it > NEWTON_MAX_ITER) break;

        // Solve the linear system.
        if(!solver->solve())
          error ("Matrix solver failed.\n");

        // Add \deltaY^{n+1} to Y^n.
        for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];
        
        if (it >= NEWTON_MAX_ITER)
          error ("Newton method did not converge.");

        it++;
      }

      // Store the result in ref_sln.
      Solution::vector_to_solutions(coeff_vec, *ref_spaces, Hermes::Tuple<Solution *>(&C_ref_sln, &phi_ref_sln));

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt* adaptivity = new Adapt(Hermes::Tuple<Space *>(&C, &phi), Hermes::Tuple<ProjNormType>(HERMES_H1_NORM, HERMES_H1_NORM));
      bool solutions_for_adapt = true;
      Hermes::Tuple<double> err_est_rel;
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::Tuple<Solution *>(&C_sln, &phi_sln), Hermes::Tuple<Solution *>(&C_ref_sln, &phi_ref_sln), solutions_for_adapt, 
                                 HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS, &err_est_rel) * 100;

      // Report results.
      info("ndof_coarse[0]: %d, ndof_fine[0]: %d",
           C.get_num_dofs(), (*ref_spaces)[0]->get_num_dofs());
      info("err_est_rel[0]: %g%%", err_est_rel[0]*100);
      info("ndof_coarse[1]: %d, ndof_fine[1]: %d",
           phi.get_num_dofs(), (*ref_spaces)[1]->get_num_dofs());
      info("err_est_rel[1]: %g%%", err_est_rel[1]*100);
      // Report results.
      info("ndof_coarse_total: %d, ndof_fine_total: %d, err_est_rel: %g%%", 
           Space::get_num_dofs(Hermes::Tuple<Space *>(&C, &phi)), Space::get_num_dofs(*ref_spaces), err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
      else 
      {
        info("Adapting the coarse mesh.");
        done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

        if (Space::get_num_dofs(Hermes::Tuple<Space *>(&C, &phi)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }
      
      info("Projecting fine mesh solution on new coarse mesh.");
        OGProjection::project_global(Hermes::Tuple<Space *>(&C, &phi), Hermes::Tuple<Solution *>(&C_ref_sln, &phi_ref_sln), Hermes::Tuple<Solution *>(&C_sln, &phi_sln));

      // Clean up.
      delete solver;
      delete matrix;
      delete rhs;
      delete adaptivity;
      if(done == false)
      for(int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i]->get_mesh();
      delete ref_spaces;
      delete dp;

    }
    while (done == false);

    // Copy last reference solution into sln_prev_time.
    C_prev_time.copy(&C_ref_sln);
    phi_prev_time.copy(&phi_ref_sln);
  }

  return ERROR_FAILURE;
}
