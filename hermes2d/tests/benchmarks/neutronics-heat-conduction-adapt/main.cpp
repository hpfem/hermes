#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include <cmath>
#include <iostream>

using namespace RefinementSelectors;

// This test makes sure that the benchmark "neutronics-heat-conduction-adapt" works correctly.

const bool SOLVE_ON_COARSE_MESH = false;   // true... Newton is done on coarse mesh in every adaptivity step.
                                           // false...Newton is done on coarse mesh only once, then projection
                                           //         of the fine mesh solution to coarse mesh is used.
const int INIT_GLOB_REF_NUM = 1;           // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;            // Number of initial refinements towards boundary.
const int P_INIT = 1;                      // Initial polynomial degree of all mesh elements

// Time-stepping:
const double TAU = 0.1;                    // Time step.
const double T_FINAL = 0.3;                // Time interval length.

// Adaptivity:
const int UNREF_FREQ = 10;                 // Every UNREF_FREQ time step the mesh is unrefined.
const double THRESHOLD = 0.3;              // This is a quantitative parameter of the adapt(...) function and
                                           // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                    // Adaptive strategy:
                                           // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                           //   error is processed. If more elements have similar errors, refine
                                           //   all to keep the mesh symmetric.
                                           // STRATEGY = 1 ... refine all elements whose error is larger
                                           //   than THRESHOLD times maximum element error.
                                           // STRATEGY = 2 ... refine all elements whose error is larger
                                           //   than THRESHOLD.
                                           // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ISO;     // Predefined list of element refinement candidates. Possible values are
                                           // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                           // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                           // See User Documentation for details.
const int MESH_REGULARITY = -1;            // Maximum allowed level of hanging nodes:
                                           // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                           // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                           // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                           // Note that regular meshes are not supported, this is due to
                                           // their notoriously bad performance.
const double CONV_EXP = 1.0;               // Default value is 1.0. This parameter influences the selection of
                                           // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_P = 6;                       // Maximum polynomial order allowed in hp-adaptivity
                                           // had to be limited due to complicated integrals.
const double ERR_STOP = 0.01;              // Stopping criterion for hp-adaptivity
                                           // (relative error between reference and coarse solution in percent).
const int NDOF_STOP = 10000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.
const int ORDER_INCREASE = 1;              // Increase in approximation order associated with the global refinement.


// Newton's method:
const double NEWTON_TOL_COARSE = 1.0e-2;   // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 5.0e-2;     // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;            // Maximum allowed number of Newton iterations.

// Linear system solvers for the coarse and refined problems, respectively.
// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
// (depending on which optional solver libraries you have installed and enabled in hermes2d/CMake.vars).
MatrixSolverType matrix_solver_coarse = SOLVER_UMFPACK;  
MatrixSolverType matrix_solver_fine = SOLVER_UMFPACK;

// Problem parameters.
const double CT = 1.0;
const double CF = 1.0;
const double rT = 1.0;
const double rF = 0.25;
const double LX = 100.0;          // Domain sizes in the x and y dimensions.
const double LY = 100.0;
const double invvel = 2.0e-4;     // Inverse of neutron velocity.
const double xsdiff = 1.268;      // Neutron diffusion coefficient.
const double Tref = 0.0;

const double nu = 2.41;           // Number of neutrons emitted per fission event.
const double xsfiss = 0.00191244; // Fission cross section.
const double kappa = 1.0e-6;      // Energy per fission.
const double rho = 1.0;           // Fuel density.
const double cp = 1.0;            // Fuel heat capacity.

// Miscellaneous:
double TIME = 0.0;                // Current time.

// Thermal conductivity dependence on temperature.
const  double k0 = 3.0e-3;
const  double k1 = 2.0e-4;
template<typename Real>
Real k(Real T) {
  return k0 + k1 * (T - Tref);
}

// Derivative of the thermal conductivity of fuel with respect to temperature.
template<typename Real>
Real dk_dT(Real T) {
  return k1;
}

// Removal cross section dependence on temperature.
const double xsa_ref = 0.0349778;
const double doppler_coeff = 1.0e-5;
template<typename Real>
Real xsrem(Real T) {
  return xsa_ref + doppler_coeff * (sqrt(T) - sqrt(Tref));
}

// Derivative of the removal cross section with respect to temperature.
template<typename Real>
Real dxsrem_dT(Real T) {
  return doppler_coeff / (2*sqrt(T));
}

// Time dependence of the temperature.
template<typename Real>
Real T_FTIME() {
//  return 1.0;
  return 1+tanh(rT*TIME);
}

template<typename Real>
Real DT_FTIME() {
//  return 0.0;
  return rT*(1-pow(tanh(rT*TIME),2));
}

// Time dependence of the neutron flux.
template<typename Real>
Real PHI_FTIME() {
//  return T_FTIME<Real>();
  return 1+exp(rF*TIME);
}

template<typename Real>
Real DPHI_FTIME() {
//  return DT_FTIME<Real>();
  return rF*exp(rF*TIME);
}

// Heat source.
template<typename Real>
Real qT(Real x, Real y) {
  Real dTdt = DT_FTIME<Real>();
  Real Tt = T_FTIME<Real>();
  Real PHIt = PHI_FTIME<Real>();
  
  Real PI_sqr = sqr(M_PI);
  Real sx = sin((M_PI*x)/LX);
  Real sy = sin((M_PI*y)/LY);
  
  return cp*CT*dTdt*rho*sx*sy - (CF*kappa*PHIt*x*xsfiss*y*sx*sy)/(LX*LY) - (-((CT*PI_sqr*Tt*sx*sy)/sqr(LX)) - (CT*PI_sqr*Tt*sx*sy)/sqr(LY)) * (k0 + k1*(-Tref + CT*Tt*sx*sy));
}

// Extraneous neutron source.
template<typename Real>
Real q(Real x, Real y) {
  Real PHIt = PHI_FTIME<Real>();
  Real dPHIdt = DPHI_FTIME<Real>();
  Real Tt = T_FTIME<Real>();
  
  Real PI_sqr = sqr(M_PI);
  Real sx = sin((M_PI*x)/LX);
  Real sy = sin((M_PI*y)/LY);
  
  return (CF*dPHIdt*invvel*x*y*sx*sy)/(LX*LY) - xsdiff*((2*CF*PHIt*M_PI*x*cos((M_PI*y)/LY)*sx)/(LX*sqr(LY)) + (2*CF*PHIt*M_PI*y*cos((M_PI*x)/LX)*sy)/(sqr(LX)*LY) - 
         (CF*PHIt*PI_sqr*x*y*sx*sy)/(LX*pow(LY,3)) - (CF*PHIt*PI_sqr*x*y*sx*sy)/(pow(LX,3)*LY)) - 
         (CF*PHIt*x*y*sx*sy*(nu*xsfiss-xsa_ref*(1 + doppler_coeff*(-sqrt(Tref) + sqrt(CT*Tt*sx*sy)))))/(LX*LY);
}

// Boundary markers.
const int BDY_DIRICHLET = 1;

// Weak forms.
#include "forms.cpp"

// Exact solutions.
#include "exact_solution.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh basemesh, mesh_T, mesh_phi;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_GLOB_REF_NUM; i++)
    basemesh.refine_all_elements();
  basemesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);
  
  // Create an individual mesh for each physical field.
  mesh_T.copy(&basemesh);
  mesh_phi.copy(&basemesh);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boudnary values.
  BCValues bc_values;
  bc_values.add_zero(BDY_DIRICHLET);

  // Create H1 spaces with default shapesets.
  H1Space space_T(&mesh_T, &bc_types, &bc_values, P_INIT);
  H1Space space_phi(&mesh_phi, &bc_types, &bc_values, P_INIT);

  Hermes::vector<Space*> spaces(&space_T, &space_phi);
  int ndof = Space::get_num_dofs(spaces); 
 
  // Solutions in the previous time step (converging within the time stepping loop).
  Solution T_prev_time, phi_prev_time;
  Hermes::vector<Solution*> prev_time_solutions(&T_prev_time, &phi_prev_time);
  
  // Solutions on the coarse and refined meshes in current time step (converging within the Newton's loop).
  Solution T_coarse, phi_coarse, T_fine, phi_fine;
  Hermes::vector<Solution*> coarse_mesh_solutions(&T_coarse, &phi_coarse);
  Hermes::vector<Solution*> fine_mesh_solutions(&T_fine, &phi_fine);
  
  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, jac_TT, jac_TT_ord);
  wf.add_matrix_form(0, 1, jac_Tphi, jac_Tphi_ord);
  wf.add_vector_form(0, res_T, res_T_ord, HERMES_ANY, &T_prev_time);
  wf.add_matrix_form(1, 0, jac_phiT, jac_phiT_ord);
  wf.add_matrix_form(1, 1, jac_phiphi, jac_phiphi_ord);
  wf.add_vector_form(1, res_phi, res_phi_ord, HERMES_ANY, &phi_prev_time);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est_T, graph_dof_exact_T, graph_cpu_est_T, graph_cpu_exact_T;
  SimpleGraph graph_dof_est_phi, graph_dof_exact_phi, graph_cpu_est_phi, graph_cpu_exact_phi;
  
  // Exact solutions for error evaluation.
  ExactSolution T_exact_solution(&mesh_T, T_exact);
  ExactSolution phi_exact_solution(&mesh_phi, phi_exact);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize the nonlinear system.
  DiscreteProblem dp(&wf, spaces, false);
  
  // Set initial conditions.
  T_prev_time.set_exact(&mesh_T, T_exact);
  phi_prev_time.set_exact(&mesh_phi, phi_exact);
  
  // Newton's loop on the initial coarse meshes.
  info("Solving on coarse meshes.");
  scalar* coeff_vec_coarse = new scalar[Space::get_num_dofs(spaces)];
  OGProjection::project_global(spaces, Hermes::vector<MeshFunction*>((MeshFunction*)&T_prev_time, (MeshFunction*)&phi_prev_time), 
                 coeff_vec_coarse, matrix_solver_coarse);
  
  // Indicate to all DiscreteProblem constructors that we solve a non-linear problem.
  bool is_linear = false;                 
  
  // Initialize the FE problem.
  DiscreteProblem dp_coarse(&wf, spaces, is_linear);

  // Setup the solvers for the coarse and fine mesh calculations, respectively.
  SparseMatrix* matrix_coarse = create_matrix(matrix_solver_coarse);
  Vector* rhs_coarse = create_vector(matrix_solver_coarse);
  Solver* solver_coarse = create_linear_solver(matrix_solver_coarse, matrix_coarse, rhs_coarse);
  //solver_coarse->set_factorization_scheme(HERMES_REUSE_MATRIX_REORDERING);
  SparseMatrix* matrix_fine = create_matrix(matrix_solver_fine);
  Vector* rhs_fine = create_vector(matrix_solver_fine);
  Solver* solver_fine = create_linear_solver(matrix_solver_fine, matrix_fine, rhs_fine);
  //solver_fine->set_factorization_scheme(HERMES_REUSE_MATRIX_REORDERING);
  
  // Perform Newton's iteration.
  info("Newton's solve on the initial coarse meshes.");
  bool verbose = false;
  if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) 
    error("Newton's iteration failed.");
  
  // Translate the resulting coefficient vector into the actual solutions. 
  Solution::vector_to_solutions(coeff_vec_coarse, spaces, coarse_mesh_solutions);

  delete [] coeff_vec_coarse;
  
  // Time stepping loop:
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= nstep; ts++)
  {
    // Update global time.
    TIME = ts*TAU;

    // Update time-dependent exact solutions.
    T_exact_solution.update(&mesh_T, T_exact);
    phi_exact_solution.update(&mesh_phi, phi_exact);

    // Periodic global derefinement.
    if (ts > 1) {
      if (ts % UNREF_FREQ == 0) {
        info("---- Time step %d - prior to adaptivity:", ts);
        info("Global mesh derefinement.");
        mesh_T.copy(&basemesh);
        mesh_phi.copy(&basemesh);
        space_T.set_uniform_order(P_INIT);
        space_phi.set_uniform_order(P_INIT);
        
        if (SOLVE_ON_COARSE_MESH) {
          // Newton's loop on the globally derefined meshes.
          scalar* coeff_vec_coarse = new scalar[Space::get_num_dofs(spaces)];
          info("Projecting previous fine mesh solution to obtain initial vector for Newton's iteration on globally derefined meshes.");
          OGProjection::project_global(spaces, Hermes::vector<MeshFunction*>((MeshFunction*)&T_fine, (MeshFunction*)&phi_fine), 
                         coeff_vec_coarse, matrix_solver_coarse);
          
          // Initialize the FE problem.
          DiscreteProblem dp_coarse(&wf, spaces, is_linear);

          // Perform Newton's iteration.
          info("Newton's solve on globally derefined meshes.");
          bool verbose = false; // Default is false.
          if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) 
            error("Newton's iteration failed.");
          
          // Translate the resulting coefficient vector into the actual solutions. 
          Solution::vector_to_solutions(coeff_vec_coarse, spaces, coarse_mesh_solutions);
          delete [] coeff_vec_coarse;
        } 
        else {
        // Projection onto the globally derefined meshes.
          info("Projecting fine mesh solutions from previous time step onto globally derefined meshes.");
          OGProjection::project_global(spaces, fine_mesh_solutions, coarse_mesh_solutions, matrix_solver_coarse); 
        }
      } 
    }

    // Adaptivity loop:
    bool done = false;
    int as = 0;
    do {
      as++;
      
      info("---- Time step %d, adaptivity step %d:", ts, as);

      // Construct globally refined reference mesh
      // and setup reference space.
      Hermes::vector<Space *>* ref_spaces = construct_refined_spaces(spaces, ORDER_INCREASE);

      // Newton's loop on the refined meshes.
      scalar* coeff_vec = new scalar[Space::get_num_dofs(*ref_spaces)];
      if (as == 1) {
        info("Projecting coarse mesh solution to obtain coefficients vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces, Hermes::vector<MeshFunction*>((MeshFunction*)&T_coarse, (MeshFunction*)&phi_coarse), 
                       coeff_vec, matrix_solver_fine);
      } else {
        info("Projecting previous fine mesh solution to obtain coefficients vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces, Hermes::vector<MeshFunction*>((MeshFunction*)&T_fine, (MeshFunction*)&phi_fine), 
                       coeff_vec, matrix_solver_fine);
        
        // Deallocate the previous fine mesh.
        delete T_fine.get_mesh();
        delete phi_fine.get_mesh();
      }
      
      // Initialize the FE problem.
      DiscreteProblem dp(&wf, *ref_spaces, is_linear);

      // Perform Newton's iteration.
      info("Newton's solve on fine meshes.");
      bool verbose = false;
      if (!solve_newton(coeff_vec, &dp, solver_fine, matrix_fine, rhs_fine, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) 
        error("Newton's iteration failed.");
            
      // Translate the resulting coefficient vector into the actual solutions. 
      Solution::vector_to_solutions(coeff_vec, *ref_spaces, fine_mesh_solutions);
      delete [] coeff_vec;
      
      if (SOLVE_ON_COARSE_MESH) {        
        // Newton's loop on the new coarse meshes.
        scalar* coeff_vec_coarse = new scalar[Space::get_num_dofs(spaces)];
        info("Projecting fine mesh solutions back onto coarse mesh to obtain initial vector for following Newton's iteration.");
        OGProjection::project_global(spaces, Hermes::vector<MeshFunction*>((MeshFunction*)&T_fine, (MeshFunction*)&phi_fine), 
                        coeff_vec_coarse, matrix_solver_coarse);
        
        // Initialize the FE problem.
        DiscreteProblem dp_coarse(&wf, spaces, is_linear);

        // Perform Newton's iteration.
        info("Newton's solve on coarse meshes.");
        bool verbose = false;
        if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) 
          error("Newton's iteration failed.");
        
        // Translate the resulting coefficient vector into the actual solutions. 
        Solution::vector_to_solutions(coeff_vec_coarse, spaces, coarse_mesh_solutions);
        delete [] coeff_vec_coarse;
      } 
      else {
        // Projection onto the new coarse meshes.
        info("Projecting fine mesh solutions back onto coarse meshes.");
        OGProjection::project_global(spaces, fine_mesh_solutions, coarse_mesh_solutions, matrix_solver_coarse); 
      }

      // Calculate element errors.
      info("Calculating error estimate and exact error."); 
      Adapt* adaptivity = new Adapt(spaces);

      // Calculate error estimate for each solution component and the total error estimate.
      Hermes::vector<double> err_est_rel;
      double err_est_rel_total = adaptivity->calc_err_est(coarse_mesh_solutions, fine_mesh_solutions, &err_est_rel) * 100;

      // Calculate exact error for each solution component and the total exact error.
      bool solutions_for_adapt = false;
      Hermes::vector<double> err_exact_rel;
      double err_exact_rel_total = adaptivity->calc_err_exact(coarse_mesh_solutions, 
                                   Hermes::vector<Solution *>(&T_exact_solution, &phi_exact_solution), 
                                   &err_exact_rel, solutions_for_adapt) * 100;

      info("T: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%, err_exact: %g %%", 
            space_T.get_num_dofs(), (*ref_spaces)[0]->get_num_dofs(), err_est_rel[0]*100, err_exact_rel[0]*100);
      info("phi: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%, err_exact: %g %%", 
            space_phi.get_num_dofs(), (*ref_spaces)[1]->get_num_dofs(), err_est_rel[1]*100, err_exact_rel[1]*100);
      
      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
      else {
        info("Adapting the coarse meshes.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector*> (&selector, &selector), THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (Space::get_num_dofs(spaces) >= NDOF_STOP) done = true; 
      }
      
      delete adaptivity;
      delete ref_spaces;
    }
    while (!done);

    // Make the fine mesh solution at current time level the previous time level solution in the following time step.
    T_prev_time.copy(&T_fine);
    phi_prev_time.copy(&phi_fine);
  }

  info("Coordinate ( 0.0,  0.0) T value = %lf", T_prev_time.get_pt_value(0.0, 0.0));
  info("Coordinate ( 25.0, 25.0) T value = %lf", T_prev_time.get_pt_value(25.0, 25.0));
  info("Coordinate ( 25.0, 75.0) T value = %lf", T_prev_time.get_pt_value(25.0, 75.0));
  info("Coordinate ( 75.0, 25.0) T value = %lf", T_prev_time.get_pt_value(75.0, 25.0));
  info("Coordinate ( 75.0, 75.0) T value = %lf", T_prev_time.get_pt_value(75.0, 75.0));

  info("Coordinate ( 0.0,  0.0) phi value = %lf", phi_prev_time.get_pt_value(0.0, 0.0));
  info("Coordinate ( 25.0, 25.0) phi value = %lf", phi_prev_time.get_pt_value(25.0, 25.0));
  info("Coordinate ( 25.0, 75.0) phi value = %lf", phi_prev_time.get_pt_value(25.0, 75.0));
  info("Coordinate ( 75.0, 25.0) phi value = %lf", phi_prev_time.get_pt_value(75.0, 25.0));
  info("Coordinate ( 75.0, 75.0) phi value = %lf", phi_prev_time.get_pt_value(75.0, 75.0));
  
  int success = 1;
  
  // Test convergence to correct results.  
  double eps = 1e-5;
  if (fabs(T_prev_time.get_pt_value(0.0, 0.0) - 0.000000) > eps) {
    info("Coordinate (  0,  0) T value = %lf\n", T_prev_time.get_pt_value(0.0, 0.0));
    success = 0;
  }
  
  if (fabs(T_prev_time.get_pt_value(25.0, 25.0) - 0.693146) > eps) {
    printf("Coordinate ( 25, 25) T value = %lf\n", T_prev_time.get_pt_value(25.0, 25.0));
    success = 0;
  }
  
  if (fabs(T_prev_time.get_pt_value(25.0, 75.0) - 0.693146) > eps) {
    printf("Coordinate ( 25, 75) T value = %lf\n", T_prev_time.get_pt_value(25.0, 75.0));
    success = 0;
  }
  
  if (fabs(T_prev_time.get_pt_value(75.0, 25.0) - 0.693146) > eps) {
    printf("Coordinate ( 75, 25) T value = %lf\n", T_prev_time.get_pt_value(75.0, 25.0));
    success = 0;
  }
  
  if (fabs(T_prev_time.get_pt_value(75.0, 75.0) - 0.693146) > eps) {
    printf("Coordinate ( 75, 75) T value = %lf\n", T_prev_time.get_pt_value(75.0, 75.0));
    success = 0;
  }
  
  if (fabs(phi_prev_time.get_pt_value(0.0, 0.0) - 0.000000) > eps) {
    printf("Coordinate (  0,  0) phi value = %lf\n", phi_prev_time.get_pt_value(0.0, 0.0));
    success = 0;
  }
  
  if (fabs(phi_prev_time.get_pt_value(25.0, 25.0) - 0.064914) > eps) {
    printf("Coordinate ( 25, 25) phi value = %lf\n", phi_prev_time.get_pt_value(25.0, 25.0));
    success = 0;
  }
  
  if (fabs(phi_prev_time.get_pt_value(25.0, 75.0) - 0.194752) > eps) {
    printf("Coordinate ( 25, 75) phi value = %lf\n", phi_prev_time.get_pt_value(25.0, 75.0));
    success = 0;
  }
  
  if (fabs(phi_prev_time.get_pt_value(75.0, 25.0) - 0.194752) > eps) {
    printf("Coordinate ( 75, 25) phi value = %lf\n", phi_prev_time.get_pt_value(75.0, 25.0));
    success = 0;
  }
  
  if (fabs(phi_prev_time.get_pt_value(75.0, 75.0) - 0.584277) > eps) {
    printf("Coordinate ( 75, 75) phi value = %lf\n", phi_prev_time.get_pt_value(75.0, 75.0));
    success = 0;
  }
  
  // Test adaptivity.
  int ndof_allowed_T = 85;
  int ndof_allowed_phi = 85;
  int ndof_T = Space::get_num_dofs(&space_T);
  int ndof_phi = Space::get_num_dofs(&space_phi);
  printf("ndof_actual_T = %d\n", ndof_T);
  printf("ndof_actual_phi = %d\n", ndof_phi);
  printf("ndof_allowed_T = %d\n", ndof_allowed_T);
  printf("ndof_allowed_phi = %d\n", ndof_allowed_phi);
  if ((ndof_T > ndof_allowed_T) || (ndof_phi > ndof_allowed_phi)) {
    // ndofs_T was 81 and ndof_phi was 81 at the time this test was created
    printf("Adaptivity failed.");
    success = 0;
  }
  
  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
