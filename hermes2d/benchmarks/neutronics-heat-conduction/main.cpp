#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include <cmath>
#include <iostream>

// Neutronics/heat conduction test case.
//
// Authors: Damien Lebrun-Grandie (Texas A&M University, USA),
//          Milan Hanus (University of West Bohemia, Pilsen, Czech Rep.).
// PDEs:
//
// 1/v d/dt phi = div(D grad phi) + (nu Sigma_f - Sigma_r(T)) phi + q
//
// rho c_p d/dt T = div(k(T) grad T) + kappa Sigma_f phi + qT
//
// Domain: rectangle (LX, LY).
//
// BC: homogeneous Dirichlet.
//
// The following parameters can be changed:

const int INIT_GLOB_REF_NUM = 2;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;                   // Number of initial refinements towards boundary.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const double TAU = 0.05;                          // Time step.
const double T_FINAL = 10.0;                      // Time interval length.
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
const double LX = 100.0;                          // Domain sizes in the x and y dimensions.
const double LY = 100.0;
const double invvel = 2.0e-4;                     // Inverse of neutron velocity.
const double xsdiff = 1.268;                      // Neutron diffusion coefficient.
const double Tref = 0.0;

const double nu = 2.41;                           // Number of neutrons emitted per fission event.
const double xsfiss = 0.00191244;                 // Fission cross section.
const double kappa = 1.0e-6;                      // Energy per fission.
const double rho = 1.0;                           // Fuel density.
const double cp = 1.0;                            // Fuel heat capacity.

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

  // Initialize solution views (their titles will be2 updated in each time step).
  ScalarView sview_T("", new WinGeom(0, 0, 500, 400));
  sview_T.fix_scale_width(50);
  ScalarView sview_phi("", new WinGeom(0, 500, 500, 400));
  sview_phi.fix_scale_width(50);
  ScalarView sview_T_exact("", new WinGeom(550, 0, 500, 400));
  sview_T_exact.fix_scale_width(50);
  ScalarView sview_phi_exact("", new WinGeom(550, 500, 500, 400));
  sview_phi_exact.fix_scale_width(50);
  char title[100]; // Character array to store the title for an actual view and time step.

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
    bool verbose = true;
    if(!solve_newton(coeff_vec, &dp, solver, matrix, rhs, NEWTON_TOL, NEWTON_MAX_ITER, verbose))
      error("Newton's iteration failed.");
        
    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solutions(coeff_vec, spaces, newton_iterates);
    delete [] coeff_vec;
    
    // Show the new time level solution.
    sprintf(title, "Approx. solution for T, t = %g s", TIME);
    sview_T.set_title(title); 
    sview_T.show(&T_prev_newton);
    
    sprintf(title, "Approx. solution for phi, t = %g s", TIME);
    sview_phi.set_title(title);
    sview_phi.show(&phi_prev_newton);

    // Exact solution for comparison with computational results.
    T_exact_solution.update(&mesh, T_exact);
    phi_exact_solution.update(&mesh, phi_exact);

    // Show exact solution.
    sview_T_exact.show(&T_exact_solution);
    sprintf(title, "Exact solution for T, t = %g s", TIME);
    sview_T_exact.set_title(title);
    
    sview_phi_exact.show(&phi_exact_solution);
    sprintf(title, "Exact solution for phi, t = %g s", TIME);
    sview_phi_exact.set_title(title);
    
    // Calculate exact error.
    info("Calculating error (exact).");
    Hermes::vector<double> exact_errors;
    Adapt adaptivity_exact(spaces, Hermes::vector<ProjNormType>(HERMES_H1_NORM, HERMES_H1_NORM));
    bool solutions_for_adapt = false;
    adaptivity_exact.calc_err_exact(Hermes::vector<Solution *>(&T_prev_newton, &phi_prev_newton), 
                                    Hermes::vector<Solution *>(&T_exact_solution, &phi_exact_solution), 
                                    &exact_errors, solutions_for_adapt);
    
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
  
  // Wait for all views to be closed.
  View::wait();

  return 0;
}
