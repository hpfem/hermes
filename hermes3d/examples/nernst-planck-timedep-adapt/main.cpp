#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"

#include "timestep_controller.h"
//#include "timestep_controller.h"

//using namespace RefinementSelectors;


MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, 

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
double INIT_TAU = 0.1;
double *TAU = &INIT_TAU;                          // Size of the time step.
const int REF_INIT = 3;     	                  // Number of initial refinements.

const int P_INIT_X = 3,
          P_INIT_Y = 3,
          P_INIT_Z = 3;                           // Initial orders.
const bool MULTIMESH = false;	                  // Multimesh?
const int TIME_DISCR = 1;                         // 1 for implicit Euler, 2 for Crank-Nicolson.

const double NEWTON_TOL_COARSE = 0.01;            // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 0.05;              // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.

const int UNREF_FREQ = 1;                         // every UNREF_FREQth time step the mesh is unrefined.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const double ERR_STOP = 5;                        // Stopping criteria
const int NDOF_STOP = 8000;                       // To prevent adaptivity from going on forever.
/*const int STRATEGY = 0;                         // Adaptive strategy:
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
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
*/
// Weak forms
#include "forms.cpp"


/*** Boundary types and conditions ***/

// Boundary markers.
const int BDY_TOP = 1;
const int BDY_BOT = 2;



scalar ic_phi(double x, double y, double z, double &dx, double &dy, double &dz) {
  return 0.0;
}

scalar ic_C(double x, double y, double z, double &dx, double &dy, double &dz) {
  return C0;
}

// Boundary condition types. 
BCType bc_types_phi(int marker) {
  if (marker == BDY_TOP || marker == BDY_BOT) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

BCType bc_types_C(int marker) {
  return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values. 
scalar essential_bc_values_phi(int ess_bdy_marker, double x, double y, double z) {
  if (ess_bdy_marker == BDY_TOP)
    return VOLTAGE;
  else
    return 0.0;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values_C(int ess_bdy_marker, double x, double y, double z) {
    return 0.0;
}


int main (int argc, char* argv[]) {
  
  // Load the mesh. 
  Mesh basemesh;
  ExodusIIReader mesh_loader;
  if (!mesh_loader.load("coarse_mesh_full.e", &basemesh))
    error("Loading mesh file '%s' failed.\n", "coarse_mesh_full.e");

  Mesh C_mesh, phi_mesh;
  C_mesh.copy(basemesh);

  phi_mesh.copy(basemesh);

  H1Space C_space(&C_mesh, bc_types_C, essential_bc_values_C, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
  H1Space phi_space(MULTIMESH ? &phi_mesh : &C_mesh, bc_types_phi, essential_bc_values_phi, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));

  Solution C_prev_time(&C_mesh);
  C_prev_time.set_const(C0);
  Solution phi_prev_time(MULTIMESH ? &phi_mesh : &C_mesh);
  phi_prev_time.set_const(0.0);


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

  int ndof = Space::get_num_dofs(Hermes::vector<Space*>(&C_space, &phi_space));

  Solution C_sln(C_space.get_mesh());
  Solution phi_sln(phi_space.get_mesh());

  info("Projecting initial condition to obtain initial vector for the Newton's method.");
  scalar* coeff_vec_coarse = new scalar[ndof];
  OGProjection::project_global(Hermes::vector<Space *>(&C_space, &phi_space),
                               Hermes::vector<MeshFunction *>(&C_prev_time, &phi_prev_time),
                               coeff_vec_coarse, matrix_solver);



  bool is_linear = false;
  DiscreteProblem dp_coarse(&wf, Hermes::vector<Space *>(&C_space, &phi_space), is_linear);

  //Solution::vector_to_solutions(coeff_vec_coarse, dp_coarse.get_spaces(), Hermes::vector<Solution *>(&C_sln, &phi_sln), NULL);


  // Set up the solver, matrix, and rhs for the coarse mesh according to the solver selection.
  SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
  Vector* rhs_coarse = create_vector(matrix_solver);
  Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

  info("Solving on coarse mesh:");
  bool verbose = true;
  if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse,
      NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");

  info("Solved!");
  // Translate the resulting coefficient vector into the Solution sln.
  Solution::vector_to_solutions(coeff_vec_coarse, Hermes::vector<Space *>(&C_space, &phi_space),
                                Hermes::vector<Solution *>(&C_sln, &phi_sln));

  out_fn_vtk(&C_sln,"C_init_sln");
  out_fn_vtk(&phi_sln,"phi_init_sln");
  //out_fn_vtk(&sln, "sln", ts);

  Solution *C_ref_sln, *phi_ref_sln;

  PidTimestepController pid(T_FINAL, false, INIT_TAU);
  TAU = pid.timestep;
  info("Starting time iteration with the step %g", *TAU);
  do {
    pid.begin_step();

    if (pid.get_timestep_number() > 1 && pid.get_timestep_number() % UNREF_FREQ == 0)
    {
      info("Global mesh derefinement.");
      C_mesh.copy(basemesh);
      if (MULTIMESH)
      {
        phi_mesh.copy(basemesh);
      }
      C_space.set_uniform_order(Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
      phi_space.set_uniform_order(Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));


      }
      bool done = false; int as = 1;
      double err_est;
      do {
        info("Time step %d, adaptivity step %d:", pid.get_timestep_number(), as);

        // Construct globally refined reference mesh
        // and setup reference space.
        int order_increase = 1;
        Hermes::vector<Space *>* ref_spaces = construct_refined_spaces(Hermes::vector<Space *>(&C_space, &phi_space), 
                                                                       order_increase);
        scalar* coeff_vec = new scalar[Space::get_num_dofs(*ref_spaces)];
        DiscreteProblem* dp = new DiscreteProblem(&wf, *ref_spaces, is_linear);
        SparseMatrix* matrix = create_matrix(matrix_solver);
        Vector* rhs = create_vector(matrix_solver);
        Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);


        if (as == 1 && pid.get_timestep_number() == 1) {
          info("Projecting coarse mesh solution to obtain coefficient vector on new fine mesh.");
          OGProjection::project_global(*ref_spaces, Hermes::vector<MeshFunction *>(&C_sln, &phi_sln),
                                       coeff_vec, matrix_solver);
        }
        else {
          info("Projecting previous fine mesh solution to obtain coefficient vector on new fine mesh.");
          OGProjection::project_global(*ref_spaces, Hermes::vector<MeshFunction *>(C_ref_sln, phi_ref_sln),
                                       coeff_vec, matrix_solver);
        }
        if (as > 1) {
          // Now deallocate the previous mesh
          info("Deallocating the previous mesh");
          //delete C_ref_sln->get_mesh();
          //delete phi_ref_sln->get_mesh();
          //delete C_ref_sln;
          //delete phi_ref_sln;
        }
        /*TODO TEMP */
        if (pid.get_timestep_number() > 1) {

        delete C_ref_sln;
        delete phi_ref_sln;
        }

        info("Solving on fine mesh:");
        if (!solve_newton(coeff_vec, dp, solver, matrix, rhs,
              NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) error("Newton's iteration failed.");


        // Store the result in ref_sln.
        C_ref_sln = new Solution(ref_spaces->at(0)->get_mesh());
        phi_ref_sln = new Solution(ref_spaces->at(1)->get_mesh());

        Solution::vector_to_solutions(coeff_vec, *ref_spaces,
                                      Hermes::vector<Solution *>(C_ref_sln, phi_ref_sln));
        // Projecting reference solution onto the coarse mesh
        info("Projecting fine mesh solution on coarse mesh.");
        OGProjection::project_global(Hermes::vector<Space *>(&C_space, &phi_space),
                                     Hermes::vector<Solution *>(C_ref_sln, phi_ref_sln),
                                     Hermes::vector<Solution *>(&C_sln, &phi_sln),
                                     matrix_solver);


        info("Calculating error estimate.");
        Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&C_space, &phi_space),
            Hermes::vector<ProjNormType> (HERMES_H1_NORM, HERMES_H1_NORM));
        Hermes::vector<double> err_est_rel;
        bool solutions_for_adapt = true;

        double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&C_sln, &phi_sln),
            Hermes::vector<Solution *>(C_ref_sln, phi_ref_sln), solutions_for_adapt,
            HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS, &err_est_rel) * 100;

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
          adaptivity->adapt(THRESHOLD);

          info("Adapted...");

          if (Space::get_num_dofs(Hermes::vector<Space *>(&C_space, &phi_space)) >= NDOF_STOP)
            done = true;
          else
            // Increase the counter of performed adaptivity steps.
            as++;
        }



        //as++;
        delete solver;
        delete matrix;
        delete rhs;
        delete ref_spaces;
        delete dp;
        delete[] coeff_vec;
        done = true;
      } while (!done);
      out_fn_vtk(C_ref_sln,"C_sln", pid.get_timestep_number());
      out_fn_vtk(phi_ref_sln,"phi_sln", pid.get_timestep_number());

      pid.end_step(Hermes::vector<Solution*> (C_ref_sln, phi_ref_sln), Hermes::vector<Solution*> (&C_prev_time, &phi_prev_time));

      // Copy last reference solution into sln_prev_time.
      C_prev_time.copy(C_ref_sln);
      phi_prev_time.copy(phi_ref_sln);
  } while (pid.has_next());


  //View::wait();
  return 0;
}

