#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include <cmath>
#include <iostream>

// This test makes sure that the benchmark "neutronics-heat-conduction" works correctly.

const int INIT_GLOB_REF_NUM = 4;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;                   // Number of initial refinements towards boundary.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const double TAU = 0.1;                           // Time step.
const double T_FINAL = 1.0;                       // Time interval length.
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
const double INIT_COND_CONST = 3.0;               // Constant initial condition.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

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
double TIME = 0.0;                                    // Current time.
const int TIME_MAX_ITER = int(T_FINAL / TAU + 0.5);   // Number of time steps.

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

// Boundary markers.
const int BDY_DIRICHLET = 1; 

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

// Weak forms.
#include "forms.cpp"

// Exact solutions.
#include "exact_solution.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  info("TIME_MAX_ITER = %d", TIME_MAX_ITER);

  // Load the mesh file.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_DIRICHLET);

  // Enter Dirichlet boudnary values.
  BCValues bc_values;
  bc_values.add_zero(BDY_DIRICHLET);

  // Create H1 spaces with default shapesets.
  H1Space space_T(&mesh, &bc_types, &bc_values, P_INIT);
  H1Space space_phi(&mesh, &bc_types, &bc_values, P_INIT);
  Hermes::vector<Space*> spaces(&space_T, &space_phi);

  // Exact solutions for error evaluation.
  ExactSolution T_exact_solution(&mesh, T_exact),
                phi_exact_solution(&mesh, phi_exact);

  // Solutions in the previous time step.
  Solution T_prev_time, phi_prev_time;
  Hermes::vector<MeshFunction*> time_iterates(&T_prev_time, &phi_prev_time);
  
  // Solutions in the previous Newton's iteration.
  Solution T_prev_newton, phi_prev_newton;
  Hermes::vector<Solution*> newton_iterates(&T_prev_newton, &phi_prev_newton);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, jac_TT, jac_TT_ord);
  wf.add_matrix_form(0, 1, jac_Tphi, jac_Tphi_ord);
  wf.add_vector_form(0, res_T, res_T_ord, HERMES_ANY, &T_prev_time);
  wf.add_matrix_form(1, 0, jac_phiT, jac_phiT_ord);
  wf.add_matrix_form(1, 1, jac_phiphi, jac_phiphi_ord);
  wf.add_vector_form(1, res_phi, res_phi_ord, HERMES_ANY, &phi_prev_time);
  
  // Set initial conditions.
  T_prev_time.set_exact(&mesh, T_exact);
  phi_prev_time.set_exact(&mesh, phi_exact);
  
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  solver->set_factorization_scheme(HERMES_REUSE_MATRIX_REORDERING);

  // Time stepping.
  int t_step = 1;
  do {
    TIME += TAU;

    info("---- Time step %d, t = %g s:", t_step, TIME); t_step++;
    info("Projecting to obtain initial vector for the Newton's method.");

    scalar* coeff_vec = new scalar[Space::get_num_dofs(spaces)];
    OGProjection::project_global(spaces, time_iterates, coeff_vec, matrix_solver);
    Solution::vector_to_solutions(coeff_vec, Hermes::vector<Space*>(&space_T, &space_phi), 
                                  Hermes::vector<Solution*>(&T_prev_newton, &phi_prev_newton));
    
    // Initialize the FE problem.
    bool is_linear = false;
    DiscreteProblem dp(&wf, spaces, is_linear);

    // Perform Newton's iteration.
    info("Newton's iteration...");
    bool verbose = false;
    if(!solve_newton(coeff_vec, &dp, solver, matrix, rhs, NEWTON_TOL, NEWTON_MAX_ITER, verbose))
      error("Newton's iteration failed.");
        
    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solutions(coeff_vec, spaces, newton_iterates);
    delete [] coeff_vec;

    // Exact solution for comparison with computational results.
    T_exact_solution.update(&mesh, T_exact);
    phi_exact_solution.update(&mesh, phi_exact);
    
    // Calculate exact error.
    info("Calculating error (exact).");
    Hermes::vector<double> exact_errors;
    Adapt adaptivity_exact(spaces);
    bool solutions_for_adapt = false;
    adaptivity_exact.calc_err_exact(Hermes::vector<Solution *>(&T_prev_newton, &phi_prev_newton), Hermes::vector<Solution *>(&T_exact_solution, &phi_exact_solution), &exact_errors, solutions_for_adapt);
    
    double maxerr = std::max(exact_errors[0], exact_errors[1])*100;
    info("Exact solution error for T (H1 norm): %g %%", exact_errors[0]*100);
    info("Exact solution error for phi (H1 norm): %g %%", exact_errors[1]*100);
    info("Exact solution error (maximum): %g %%", maxerr);
    
    // Prepare previous time level solution for the next time step.
    T_prev_time.copy(&T_prev_newton);
    phi_prev_time.copy(&phi_prev_newton);
  }
  while (t_step <= TIME_MAX_ITER);

  // Cleanup.
  delete matrix;
  delete rhs;
  delete solver;
  
  info("Coordinate (  0,  0) T value = %lf", T_prev_time.get_pt_value(0.0, 0.0));
  info("Coordinate ( 25, 25) T value = %lf", T_prev_time.get_pt_value(25.0, 25.0));
  info("Coordinate ( 75, 25) T value = %lf", T_prev_time.get_pt_value(75.0, 25.0));
  info("Coordinate ( 25, 75) T value = %lf", T_prev_time.get_pt_value(25.0, 75.0));
  info("Coordinate ( 75, 75) T value = %lf", T_prev_time.get_pt_value(75.0, 75.0));

  info("Coordinate (  0,  0) phi value = %lf", phi_prev_time.get_pt_value(0.0, 0.0));
  info("Coordinate ( 25, 25) phi value = %lf", phi_prev_time.get_pt_value(25.0, 25.0));
  info("Coordinate ( 75, 25) phi value = %lf", phi_prev_time.get_pt_value(75.0, 25.0));
  info("Coordinate ( 25, 75) phi value = %lf", phi_prev_time.get_pt_value(25.0, 75.0));
  info("Coordinate ( 75, 75) phi value = %lf", phi_prev_time.get_pt_value(75.0, 75.0));
 
  int success = 1;
  double eps = 1e-5;
  if (fabs(T_prev_time.get_pt_value(0.0, 0.0) - 0.000000) > eps) {
    printf("Coordinate (  0,  0) T value = %lf\n", T_prev_time.get_pt_value(0.0, 0.0));
    success = 0;
  }

  if (fabs(T_prev_time.get_pt_value(25.0, 25.0) - 0.915885) > eps) {
    printf("Coordinate ( 25, 25) T value = %lf\n", T_prev_time.get_pt_value(25.0, 25.0));
    success = 0;
  }

  if (fabs(T_prev_time.get_pt_value(75.0, 25.0) - 0.915885) > eps) {
    printf("Coordinate ( 75, 25) T value = %lf\n", T_prev_time.get_pt_value(75.0, 25.0));
    success = 0;
  }

  if (fabs(T_prev_time.get_pt_value(25.0, 75.0) - 0.915885) > eps) {
    printf("Coordinate ( 25, 75) T value = %lf\n", T_prev_time.get_pt_value(25.0, 75.0));
    success = 0;
  }

  if (fabs(T_prev_time.get_pt_value(75.0, 75.0) - 0.915885) > eps) {
    printf("Coordinate ( 75, 75) T value = %lf\n", T_prev_time.get_pt_value(75.0, 75.0));
    success = 0;
  }

  if (fabs(phi_prev_time.get_pt_value(0.0, 0.0) - 0.000000) > eps) {
    printf("Coordinate (  0,  0) phi value = %lf\n", phi_prev_time.get_pt_value(0.0, 0.0));
    success = 0;
  }

  if (fabs(phi_prev_time.get_pt_value(25.0, 25.0) - 0.071349) > eps) {
    printf("Coordinate ( 25, 25) phi value = %lf\n", phi_prev_time.get_pt_value(25.0, 25.0));
    success = 0;
  }

  if (fabs(phi_prev_time.get_pt_value(75.0, 25.0) - 0.214063) > eps) {
    printf("Coordinate ( 75, 25) phi value = %lf\n", phi_prev_time.get_pt_value(75.0, 25.0));
    success = 0;
  }

  if (fabs(phi_prev_time.get_pt_value(25.0, 75.0) - 0.214063) > eps) {
    printf("Coordinate ( 25, 75) phi value = %lf\n", phi_prev_time.get_pt_value(25.0, 75.0));
    success = 0;
  }

  if (fabs(phi_prev_time.get_pt_value(75.0, 75.0) - 0.642226) > eps) {
    printf("Coordinate ( 75, 75) phi value = %lf\n", phi_prev_time.get_pt_value(75.0, 75.0));
    success = 0;
  }

  if (success == 1) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
