#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function.h"
#include "math.h"
#include <iostream>


// Neutronics/heat conduction test case.
//
// Author: Damien Lebrun-Grandie (Texas A&M University).
//
// PDE:
//
// 1/v1 d/dt phi_1 = div(D_1 grad phi_1) + (nu Sigma_f1 - Sigma_r1) phi_1 + nu Sigma_f2 phi_2 
//                        + lambda_1 C1 + lambda_2 C2 + Q_1
//
// 1/v2 d/dt phi_2 = div(D_2 grad phi_2) - Sigma_r2 phi_2 + Sigma_12 phi_1 + Q_2
//
// d/dt C_i = beta_1i nu Sigma_f1 phi_1 + beta_2i nu Sigma_f2 phi_2 - lambda_i C_i    i = 1,2
//
// rho c_p d/dt T = div(k grad T) + kappa_1 Sigma_f1 phi_1 + kappa_2 Sigma_f2 phi_2 + Q_T
//
// Domain: rectangle (Lx, Ly).
//
// BC: homogeneous Dirichlet.
// 
// The following parameters can be changed:
//

const int INIT_GLOB_REF_NUM = 4;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;                   // Number of initial refinements towards boundary.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const double TAU = 0.1;                           // Time step.
const double T_FINAL = 10.0;                      // Time interval length.
const double NEWTON_TOL = 1e-6;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
const double INIT_COND_CONST = 3.0;               // Constant initial condition.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Problem parameters.

const double CT = 1.0;
const double CF = 1.0;
const double rT = 1.0;
const double rF = 0.25;
const double LX = 100.0;          // Domain sizes in the x and y dimensions
const double LY = 100.0;
const double invvel = 2.0e-4;     // Inverse of neutron velocity
const double xsdiff = 1.268;      // Diffusion coefficient
const double Tref = 0.0;          // Temperature at boundary

const double nu = 2.41;           // Number of neutrons emitted per fission event
const double xsfiss = 0.00191244; // Fission cross section
const double kappa = 1.0e-6;
const double rho = 1.0;           // Density
const double cp = 1.0;            // Heat capacity

const double PI = acos(-1.0);
const double normalization_const = 1.0;

const double energy_per_fission = kappa * xsfiss;

// Internal variables.
double TIME = 0.0;                              // Current time.
const int TIME_MAX_ITER = int(T_FINAL / TAU);   // Number of time steps.

// Thermal conductivity depends on temperature
const  double k0 = 3.0e-3;
const  double k1 = 2.0e-4;
template<typename Real>
Real k(Real T) {
  return k0 + k1 * (T - Tref);
}

// Derivative of the thermal conductivity
template<typename Real>
Real dk_dT(Real T) {
  return k1;
}

// Removal cross section depends on temperature
const double xsa_ref = 0.0349778;
const double doppler_coeff = 1.0e-5;
template<typename Real>
Real xsrem(Real T) {
  return xsa_ref + doppler_coeff * (sqrt(T) - sqrt(Tref));
}

// Derivative of the removal cross section with respect to temperature
template<typename Real>
Real dxsrem_dT(Real T) {
  return doppler_coeff / (2*sqrt(T+1.0e-10));
}

// Heat source
template<typename Real>
Real qT(Real x, Real y) {
  //  std::cout<<"entering qt"<<std::endl;
  return
rho*cp*CT*(1.0-pow(tanh(rT*TIME),2.0))*rT*sin(x/LX*PI)*sin(y/LY*PI)-k1*CT*CT*pow(1.0+tanh(rT*TIME),2.0)*pow(cos(x/LX*PI),2.0)/(LX*LX)*PI*PI*pow(sin(y/LY*PI),2.0)+(k0+k1*(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)-Tref))*CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)/(LX*LX)*PI*PI*sin(y/LY*PI)-k1*CT*CT*pow(1.0+tanh(rT*TIME),2.0)*pow(sin(x/LX*PI),2.0)*pow(cos(y/LY*PI),2.0)/(LY*LY)*PI*PI+(k0+k1*(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)-Tref))*CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)/(LY*LY)*PI*PI-normalization_const*energy_per_fission*xsfiss*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY;
}

// Extraneous neutron source
template<typename Real>
Real q(Real x, Real y) {
  //  std::cout<<"entering q"<<std::endl;
  return 
invvel*CF*rF*exp(rF*TIME)*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY-xsdiff*(-CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)/(LX*LX*LX)*PI*PI*sin(y/LY*PI)*x*y/LY+2.0*CF*(1.0+exp(rF*TIME))*cos(x/LX*PI)/(LX*LX)*PI*sin(y/LY*PI)*y/LY)-xsdiff*(-CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)/(LY*LY*LY)*PI*PI*x/LX*y+2.0*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*cos(y/LY*PI)/(LY*LY)*PI*x/LX)+(xsa_ref+doppler_coeff*(sqrt(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI))-sqrt(Tref)))*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY-nu*xsfiss*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY;
}

// Boundary condition types.
BCType bc_types_T(int marker)
{
  //return BC_ESSENTIAL;
  return BC_NATURAL;
}
 

// Dirichlet boundary condition values
scalar essential_bc_values_T(int ess_bdy_marker, double x, double y)
{
  return Tref;
}


BCType bc_types_phi(int marker)
{
  //return BC_ESSENTIAL;
  return BC_NATURAL;
}
 
scalar essential_bc_values_phi(int ess_bdy_marker, double x, double y)
{
  return 0.0;
}

// Weak forms.
#include "forms.cpp"

// Exact solutions.
# include "exact_solution.cpp"

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

  // Create H1 spaces with default shapesets.
  H1Space space_T(&mesh, bc_types_T, essential_bc_values_T, P_INIT);
  H1Space space_phi(&mesh, bc_types_phi, essential_bc_values_phi, P_INIT);
  Hermes::Tuple<Space*> spaces(&space_T, &space_phi);

  // Exact solutions for error evaluation.
  ExactSolution T_exact_solution(&mesh, T_exact),
                phi_exact_solution(&mesh, phi_exact);

  // Initialize solution views (their titles will be2 updated in each time step).
  ScalarView sview_T("", new WinGeom(0, 0, 500, 400));
  ScalarView sview_phi("", new WinGeom(0, 500, 500, 400));
  ScalarView sview_T_exact("", new WinGeom(550, 0, 500, 400));
  ScalarView sview_phi_exact("", new WinGeom(550, 500, 500, 400));
  char title[100]; // Character array to store the title for an actual view and time step.

  // Solutions in the previous time step.
  Solution T_prev_time, phi_prev_time;
  Hermes::Tuple<MeshFunction*> time_iterates(&T_prev_time, &phi_prev_time);
  
  // Solutions in the previous Newton's iteration.
  Solution T_prev_newton, phi_prev_newton;
  Hermes::Tuple<Solution*> newton_iterates(&T_prev_newton, &phi_prev_newton);

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

  // Time stepping.
  int t_step = 1;
  do {
    TIME += TAU;

    info("---- Time step %d, t = %g s:", t_step, TIME); t_step++;
    info("Projecting to obtain initial vector for the Newton's method.");

    scalar* coeff_vec = new scalar[Space::get_num_dofs(spaces)];
    OGProjection::project_global(spaces, time_iterates, coeff_vec, matrix_solver);
    Solution::vector_to_solutions(coeff_vec, Hermes::Tuple<Space*>(&space_T, &space_phi), 
                                  Hermes::Tuple<Solution*>(&T_prev_newton, &phi_prev_newton));

    // Newton's method.
    info("Newton's iteration...");
    bool verbose = true; // Default is false.


    // Initialize the FE problem.
    bool is_linear = false;
    DiscreteProblem dp(&wf, spaces, is_linear);

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Perform Newton's iteration.
    int it = 1;
    while (1)
    {
      // Obtain the number of degrees of freedom.
      int ndof = Space::get_num_dofs(spaces);

      // Assemble the Jacobian matrix and residual vector.
      dp.assemble(coeff_vec, matrix, rhs, false);

      // Multiply the residual vector with -1 since the matrix 
      // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
      for (int i = 0; i < ndof; i++) rhs->set(i, -rhs->get(i));
      
      // Calculate the l2-norm of residual vector.
      double res_l2_norm = get_l2_norm(rhs);

      // Info for user.
      info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(spaces), res_l2_norm);

      // If l2 norm of the residual vector is within tolerance, or the maximum number 
      // of iteration has been reached, then quit.
      if (res_l2_norm < NEWTON_TOL || it > NEWTON_MAX_ITER) break;

      // Solve the linear system.
      if(!solver->solve())
        error ("Matrix solver failed.\n");

      // Add \deltaY^{n+1} to Y^n.
      for (int i = 0; i < ndof; i++) coeff_vec[i] += solver->get_solution()[i];
      
      if (it >= NEWTON_MAX_ITER)
        error ("Newton method did not converge.");

      it++;
    }

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solutions(coeff_vec, spaces, newton_iterates);

    // Cleanup.
    delete [] coeff_vec;
    delete matrix;
    delete rhs;
    delete solver;
    
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
    Hermes::Tuple<double> exact_errors;
    Adapt adaptivity_exact(spaces, Hermes::Tuple<ProjNormType>(HERMES_H1_NORM, HERMES_H1_NORM));
    bool solutions_for_adapt = true;
    adaptivity_exact.calc_err_exact(Hermes::Tuple<Solution *>(&T_prev_newton, &phi_prev_newton), Hermes::Tuple<Solution *>(&T_exact_solution, &phi_exact_solution), solutions_for_adapt, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL, &exact_errors);
    
    double error = std::max(exact_errors[0], exact_errors[1]);
    info("Exact solution error for T (H1 norm): %g %%", exact_errors[0]);
    info("Exact solution error for phi (H1 norm): %g %%", exact_errors[1]);
    info("Exact solution error (maximum): %g %%", error);
    
    // Prepare previous time level solution for the next time step.
    T_prev_time.copy(&T_prev_newton);
    phi_prev_time.copy(&phi_prev_newton);
  }
  while (t_step <= TIME_MAX_ITER);

  // Wait for all views to be closed.
  View::wait();

  return 0;
}
