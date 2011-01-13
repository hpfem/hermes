#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"

#include "hermes2d.h"

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

#define SIDE_MARKER 1
#define TOP_MARKER 2
#define BOT_MARKER 3
#define NONCIRCULAR

/*** Fundamental coefficients ***/
const double D = 10e-11; 	            // [m^2/s] Diffusion coefficient
const double R = 8.31; 		            // [J/mol*K] Gas constant
const double T = 293; 		            // [K] Aboslute temperature
const double F = 96485.3415;	        // [s * A / mol] Faraday constant
const double eps = 2.5e-2; 	          // [F/m] Electric permeability
const double mu = D / (R * T);        // Mobility of ions
const double z = 1;		                // Charge number
const double K = z * mu * F;          // Constant for equation
const double L =  F / eps;	          // Constant for equation
const double VOLTAGE = 1;	            // [V] Applied voltage
const scalar C0 = 1200;	              // [mol/m^3] Anion and counterion concentration

/* For Neumann boundary */
const double height = 180e-6;	              // [m] thickness of the domain
const double E_FIELD = VOLTAGE / height;    // Boundary condtion for positive voltage electrode


/* Simulation parameters */
const int PROJ_TYPE = 1;              // For the projection of the initial condition 
                                      // on the initial mesh: 1 = H1 projection, 0 = L2 projection
const double TAU = 0.05;              // Size of the time step
const int T_FINAL = 5;                // Final time
const int P_INIT = 3;       	        // Initial polynomial degree of all mesh elements.
const int REF_INIT = 5;     	        // Number of initial refinements
const bool MULTIMESH = false;	        // Multimesh?
const int TIME_DISCR = 1;             // 1 for implicit Euler, 2 for Crank-Nicolson
const int VOLT_BOUNDARY = 1;          // 1 for Dirichlet, 2 for Neumann

/* Nonadaptive solution parameters */
const double NEWTON_TOL = 1e-6;       // Stopping criterion for nonadaptive solution
const int NEWTON_MAX_ITER = 20;       // Maximum allowed number of Newton iterations

// Weak forms
#include "../forms.cpp"

/*** Boundary types and conditions ***/

// Poisson takes Dirichlet and Neumann boundaries
BCType phi_bc_types(int marker) {
  return (marker == SIDE_MARKER)
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

scalar voltage_ic(double x, double y, double &dx, double &dy) {
  // y^2 function for the domain.
  //return (y+100e-6) * (y+100e-6) / (40000e-12);
  return 0.0;
}

scalar concentration_ic(double x, double y, double &dx, double &dy) {
  return C0;
}



int main (int argc, char* argv[]) {
  // load the mesh file
  Mesh Cmesh, phimesh, basemesh;

  H2DReader mloader;

#ifdef CIRCULAR
  mloader.load("../circular.mesh", &basemesh);
  basemesh.refine_all_elements(0);
  basemesh.refine_towards_boundary(TOP_MARKER, 2);
  basemesh.refine_towards_boundary(BOT_MARKER, 4);
  MeshView mview("Mesh", 0, 600, 400, 400);
  mview.show(&basemesh);
#else
  mloader.load("../small.mesh", &basemesh);
  basemesh.refine_all_elements(1);
  //basemesh.refine_all_elements(1); // when only p-adapt is used
  //basemesh.refine_all_elements(1); // when only p-adapt is used
  basemesh.refine_towards_boundary(TOP_MARKER, REF_INIT);
  //basemesh.refine_towards_boundary(BOT_MARKER, (REF_INIT - 1) + 8); // when only p-adapt is used
  basemesh.refine_towards_boundary(BOT_MARKER, REF_INIT - 1);
#endif

  Cmesh.copy(&basemesh);
  phimesh.copy(&basemesh);

  // Spaces for concentration and the voltage
  H1Space Cspace(&Cmesh, C_bc_types, C_essential_bc_values, P_INIT);
  H1Space phispace(MULTIMESH ? &phimesh : &Cmesh, phi_bc_types, phi_essential_bc_values, P_INIT);

  // The weak form for 2 equations
  WeakForm wf(2);

  Solution C_prev_time,    // prveious time step solution, for the time integration
    C_prev_newton,   // solution convergin during the Newton's iteration
    phi_prev_time,
    phi_prev_newton;

  // Add the bilinear and linear forms
  // generally, the equation system is described:
  if (TIME_DISCR == 1) {  // Implicit Euler.
    wf.add_vector_form(0, callback(Fc_euler), H2D_ANY,
		  Hermes::vector<MeshFunction*>(&C_prev_time, &C_prev_newton, &phi_prev_newton));
    wf.add_vector_form(1, callback(Fphi_euler), H2D_ANY, Hermes::vector<MeshFunction*>(&C_prev_newton, &phi_prev_newton));
    wf.add_matrix_form(0, 0, callback(J_euler_DFcDYc), H2D_UNSYM, H2D_ANY, &phi_prev_newton);
    wf.add_matrix_form(0, 1, callback(J_euler_DFcDYphi), H2D_UNSYM, H2D_ANY, &C_prev_newton);
    wf.add_matrix_form(1, 0, callback(J_euler_DFphiDYc), H2D_UNSYM);
    wf.add_matrix_form(1, 1, callback(J_euler_DFphiDYphi), H2D_UNSYM);
  } else {
    wf.add_vector_form(0, callback(Fc_cranic), H2D_ANY, 
		  Hermes::vector<MeshFunction*>(&C_prev_time, &C_prev_newton, &phi_prev_newton, &phi_prev_time));
    wf.add_vector_form(1, callback(Fphi_cranic), H2D_ANY, Hermes::vector<MeshFunction*>(&C_prev_newton, &phi_prev_newton));
    wf.add_matrix_form(0, 0, callback(J_cranic_DFcDYc), H2D_UNSYM, H2D_ANY, Hermes::vector<MeshFunction*>(&phi_prev_newton, &phi_prev_time));
    wf.add_matrix_form(0, 1, callback(J_cranic_DFcDYphi), H2D_UNSYM, H2D_ANY, Hermes::vector<MeshFunction*>(&C_prev_newton, &C_prev_time));
    wf.add_matrix_form(1, 0, callback(J_cranic_DFphiDYc), H2D_UNSYM);
    wf.add_matrix_form(1, 1, callback(J_cranic_DFphiDYphi), H2D_UNSYM);
  }

  // Nonlinear solver
  NonlinSystem nls(&wf, Hermes::vector<Space*>(&Cspace, &phispace));

  phi_prev_time.set_exact(MULTIMESH ? &phimesh : &Cmesh, voltage_ic);
  C_prev_time.set_exact(&Cmesh, concentration_ic);

  C_prev_newton.copy(&C_prev_time);
  phi_prev_newton.copy(&phi_prev_time);

  // Project the function init_cond() on the FE space
  // to obtain initial coefficient vector for the Newton's method.
  info("Projecting initial conditions to obtain initial vector for the Newton'w method.");
  nls.project_global(Hermes::vector<MeshFunction*>(&C_prev_time, &phi_prev_time), 
      Hermes::vector<Solution*>(&C_prev_newton, &phi_prev_newton));
  
  
  //VectorView vview("electric field [V/m]", 0, 0, 600, 600);
  ScalarView Cview("Concentration [mol/m3]", 0, 0, 800, 800);
  ScalarView phiview("Voltage [V]", 650, 0, 600, 600);
  phiview.show(&phi_prev_time);
  Cview.show(&C_prev_time);
  char title[100];

  int nstep = (int) (T_FINAL / TAU + 0.5);
  for (int n = 1; n <= nstep; n++) {
    verbose("\n---- Time step %d ----", n);
    bool verbose = true; // Default is false.
    if (!nls.solve_newton(Hermes::vector<Solution*>(&C_prev_newton, &phi_prev_newton),
        NEWTON_TOL, NEWTON_MAX_ITER, verbose)) 
          error("Newton's method did not converge.");
    sprintf(title, "time step = %i", n);
    phiview.set_title(title);
    phiview.show(&phi_prev_newton);
    Cview.set_title(title);
    Cview.show(&C_prev_newton);
    phi_prev_time.copy(&phi_prev_newton);
    C_prev_time.copy(&C_prev_newton);
  }
  View::wait();

  return 0;
}
/// \}
