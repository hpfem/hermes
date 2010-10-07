#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include <cmath>
#include <iostream>

using namespace RefinementSelectors;

// Neutronics/heat conduction test case (adaptive).
//
// Author: Damien Lebrun-Grandie (Texas A&M University).
//
// PDE:
//
// 1/v d/dt phi = div(D grad phi) + nu Sigma_f phi_1 + q
//
// rho c_p d/dt T = div(k grad T) + kappa Sigma_f phi + qT
//
// Domain: rectangle (Lx, Ly).
//
// BC: homogeneous Dirichlet.
//
// The following parameters can be changed:

const bool SOLVE_ON_COARSE_MESH = false;   // true... Newton is done on coarse mesh in every adaptivity step.
                                           // false...Newton is done on coarse mesh only once, then projection
                                           //         of the fine mesh solution to coarse mesh is used.
const int INIT_GLOB_REF_NUM = 2;           // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 0;            // Number of initial refinements towards boundary.
const int P_INIT = 2;                      // Initial polynomial degree of all mesh elements
const int PROJ_TYPE = 1;                   // For the projection of the initial condition.
                                           // on the initial mesh: 1 = H1 projection, 0 = L2 projection.
// Time-stepping:
const double TAU = 0.1;                    // Time step.
const double T_FINAL = TAU;                // Time interval length.

// Adaptivity:
const int UNREF_FREQ = 1;                  // Every UNREF_FREQ time step the mesh is unrefined.
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
const CandList CAND_LIST = H2D_HP_ANISO;   // Predefined list of element refinement candidates. Possible values are
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
const double ERR_STOP = 0.001;             // Stopping criterion for hp-adaptivity
                                           // (relative error between reference and coarse solution in percent).
const int NDOF_STOP = 100000;              // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.
const int ORDER_INCREASE = 1;              // The two following parameters are used in the constructor of the class RefSystem
const int REFINEMENT = 1;                  //   Default values are 1


// Newton's method:
const double NEWTON_TOL_COARSE = 1.0e-6;   // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 5.0e-6;     // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 100;           // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Problem parameters.
const double CT = 1.0;
const double CF = 1.0;
const double rT = 1.0;
const double rF = 0.25;
const double LX = 100.0;          // Domain sizes in the x and y dimensions.
const double LY = 100.0;
const double invvel = 2.0e-4;     // Inverse of neutron velocity.
const double xsdiff = 1.268;      // Diffusion coefficient.
const double Tref = 0.0;

const double nu = 2.41;           // Number of neutrons emitted per fission event.
const double xsfiss = 0.00191244; // Fission cross section.
const double kappa = 1.0e-6;
const double rho = 1.0;           // Density.
const double cp = 1.0;            // Heat capacity.

const double PI = acos(-1.0);
const double normalization_const = 1.0;

const double energy_per_fission = kappa * xsfiss;
const double PI_ = PI;

// Miscellaneous:
double TIME = 0.0;                // Current time.

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

// Derivative of the removal cross section with respect to temperature.
template<typename Real>
Real dxsrem_dT(Real T) {
  return doppler_coeff / (2*sqrt(T));
}

// Time dependence of the temperature.
template<typename Real>
Real T_FTIME(Real x, Real y) {
  return 1.0;
  return 1+tanh(rT*TIME);
}

template<typename Real>
Real DT_FTIME(Real x, Real y) {
  return 0.0;
  return rT*(1-pow(tanh(rT*TIME),2));
}

// Time dependence of the neutron flux.
template<typename Real>
Real PHI_FTIME(Real x, Real y) {
  return T_FTIME(x, y);
  return 1+exp(rF*TIME);
}

template<typename Real>
Real DPHI_FTIME(Real x, Real y) {
  return DT_FTIME(x, y);
  return rF*(1+exp(rF*TIME));
}

// Heat source.
template<typename Real>
Real qT(Real x, Real y) {
  return CT*DT_FTIME(x,y)*cp*rho*sin((PI_*x)/LX)*sin((PI_*y)/LY)+CT*1/(LX*LX)*(PI_*PI_)*T_FTIME(x,y)*sin((PI_*x)/LX)*sin((PI_*y)/LY)*(k0-k1*(Tref-CT*T_FTIME(x,y)*sin((PI_*x)/LX)*sin((PI_*y)/LY)))+CT*1/(LY*LY)*(PI_*PI_)*T_FTIME(x,y)*sin((PI_*x)/LX)*sin((PI_*y)/LY)*(k0-k1*(Tref-CT*T_FTIME(x,y)*sin((PI_*x)/LX)*sin((PI_*y)/LY)))-(CT*CT)*1/(LX*LX)*(PI_*PI_)*(T_FTIME(x,y)*T_FTIME(x,y))*k1*pow(cos((PI_*x)/LX),2.0)*pow(sin((PI_*y)/LY),2.0)-(CT*CT)*1/(LY*LY)*(PI_*PI_)*(T_FTIME(x,y)*T_FTIME(x,y))*k1*pow(cos((PI_*y)/LY),2.0)*pow(sin((PI_*x)/LX),2.0)-(CF*PHI_FTIME(x,y)*kappa*x*xsfiss*y*sin((PI_*x)/LX)*sin((PI_*y)/LY))/(LX*LY);
}

// Extraneous neutron source.
template<typename Real>
Real q(Real x, Real y) {
  return -xsdiff*((CF*1/(LY*LY)*PHI_FTIME(x,y)*PI_*x*cos((PI_*y)/LY)*sin((PI_*x)/LX)*2.0)/LX-(CF*1/(LY*LY*LY)*PHI_FTIME(x,y)*(PI_*PI_)*x*y*sin((PI_*x)/LX)*sin((PI_*y)/LY))/LX)-xsdiff*((CF*1/(LX*LX)*PHI_FTIME(x,y)*PI_*y*cos((PI_*x)/LX)*sin((PI_*y)/LY)*2.0)/LY-(CF*1/(LX*LX*LX)*PHI_FTIME(x,y)*(PI_*PI_)*x*y*sin((PI_*x)/LX)*sin((PI_*y)/LY))/LY)+(CF*DPHI_FTIME(x,y)*invvel*x*y*sin((PI_*x)/LX)*sin((PI_*y)/LY))/(LX*LY)+(CF*PHI_FTIME(x,y)*x*y*sin((PI_*x)/LX)*sin((PI_*y)/LY)*(xsa_ref-doppler_coeff*(sqrt(Tref)-sqrt(CT*T_FTIME(x,y)*sin((PI_*x)/LX)*sin((PI_*y)/LY)))))/(LX*LY)-(CF*PHI_FTIME(x,y)*nu*x*xsfiss*y*sin((PI_*x)/LX)*sin((PI_*y)/LY))/(LX*LY);
}

// Boundary condition types.
BCType bc_types_T(int marker)
{
  return BC_ESSENTIAL;
}

BCType bc_types_phi(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values_T(int ess_bdy_marker, double x, double y)
{
  return 0.0;
}
 
scalar essential_bc_values_phi(int ess_bdy_marker, double x, double y)
{
  return 0.0;
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

  // Load the mesh.
  Mesh basemesh, mesh_T, mesh_phi;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_GLOB_REF_NUM; i++)
    basemesh.refine_all_elements();
  basemesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);
  
  // Create a special mesh for each physical field.
  mesh_T.copy(&basemesh);
  mesh_phi.copy(&basemesh);

  // Create H1 spaces with default shapesets.
  H1Space space_T(&mesh_T, bc_types_T, essential_bc_values_T, P_INIT);
  H1Space space_phi(&mesh_phi, bc_types_phi, essential_bc_values_phi, P_INIT);
  Tuple<Space*> spaces(&space_T, &space_phi);
  int ndof = Space::get_num_dofs(spaces); 
 
  // Solutions in the previous time step (converging within the time stepping loop).
  Solution T_prev_time, phi_prev_time;
  Tuple<Solution*> prev_time_solutions(&T_prev_time, &phi_prev_time);
  
  // Solutions on the coarse and refined meshes in current time step (converging within the Newton's loop).
  Solution T_coarse, phi_coarse, T_fine, phi_fine;
  Tuple<Solution*> coarse_mesh_solutions(&T_coarse, &phi_coarse);
  Tuple<Solution*> fine_mesh_solutions(&T_fine, &phi_fine);
  
  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, jac_TT, jac_TT_ord);
  wf.add_matrix_form(0, 1, jac_Tphi, jac_Tphi_ord);
  wf.add_vector_form(0, res_T, res_T_ord, HERMES_ANY, &T_prev_time);
  wf.add_matrix_form(1, 0, jac_phiT, jac_phiT_ord);
  wf.add_matrix_form(1, 1, jac_phiphi, jac_phiphi_ord);
  wf.add_vector_form(1, res_phi, res_phi_ord, HERMES_ANY, &phi_prev_time);

  // Initialize solution views (their titles will be updated in each time step).
  ScalarView view_T("", 460, 0, 450, 350);
  view_T.fix_scale_width(80);
  ScalarView view_T_exact("", 0, 0, 450, 350);
  view_T_exact.fix_scale_width(80);
  view_T_exact.show_mesh(false);
  ScalarView view_phi("", 460, 400, 450, 350);
  view_phi.fix_scale_width(80);
  ScalarView view_phi_exact("", 0, 400, 450, 350);
  view_phi_exact.fix_scale_width(80);
  view_phi_exact.show_mesh(false);

  // Initialize mesh views (their titles will be updated in each time step).
  OrderView ordview_T_coarse("", 920, 0, 450, 350);
  ordview_T_coarse.fix_scale_width(80);
  OrderView ordview_phi_coarse("", 920, 400, 450, 350);
  ordview_phi_coarse.fix_scale_width(80);
  
  char title[100]; // Character array to store the title for an actual view and time step.

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est_T, graph_dof_exact_T, graph_cpu_est_T, graph_cpu_exact_T;
  SimpleGraph graph_dof_est_phi, graph_dof_exact_phi, graph_cpu_est_phi, graph_cpu_exact_phi;
  
  // Exact solutions for error evaluation.
  ExactSolution T_exact_solution(&mesh_T, T_exact);
  ExactSolution phi_exact_solution(&mesh_phi, phi_exact);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize the nonlinear system.
  FeProblem dp(&wf, spaces);
  Tuple<ProjNormType> proj_norms(HERMES_H1_NORM, HERMES_H1_NORM);
  
  // Set initial conditions.
  T_prev_time.set_exact(&mesh_T, T_exact);
  phi_prev_time.set_exact(&mesh_phi, phi_exact);
  
  // Newton's loop on the initial coarse meshes.
  info("Solving on coarse meshes.");
  scalar* coeff_vec = new scalar[Space::get_num_dofs(spaces)];
  OGProjection::project_global(spaces, Tuple<MeshFunction*>((MeshFunction*)&T_prev_time, (MeshFunction*)&phi_prev_time), 
                 coeff_vec, matrix_solver, proj_norms);
  bool verbose = true; // Default is false.
  bool did_converge = solve_newton(spaces, &wf, coeff_vec, matrix_solver, 
                                   NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose); 
  if (!did_converge)
    error("Newton's method did not converge.");

  // Translate the resulting coefficient vector into the actual solutions. 
  Solution::vector_to_solutions(coeff_vec, Tuple<Space*>(&space_T, &space_phi), 
                                Tuple<Solution*>(&T_coarse, &phi_coarse));  
  
  // Time stepping loop:
  int nstep = (int)(T_FINAL/TAU + 0.5);
  for(int ts = 1; ts <= nstep; ts++)
  {
    // Update global time.
    TIME = ts*TAU;

    // Update time-dependent exact solutions.
    T_exact_solution.update(&mesh_T, T_exact);
    phi_exact_solution.update(&mesh_phi, phi_exact);
   
    // Show exact solution.
    view_T_exact.show(&T_exact_solution);
    sprintf(title, "T (exact), t = %g s", TIME);
    view_T_exact.set_title(title);    
    view_phi_exact.show(&phi_exact_solution);
    sprintf(title, "phi (exact), t = %g s", TIME);
    view_phi_exact.set_title(title);

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
          info("Solving on globally derefined meshes, starting from the latest fine mesh solutions.");
          OGProjection::project_global(spaces, Tuple<MeshFunction*>((MeshFunction*)&T_fine, (MeshFunction*)&phi_fine), 
                         coeff_vec, matrix_solver, proj_norms);
          did_converge = solve_newton(spaces, &wf, coeff_vec, matrix_solver, 
                                      NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose); 
          if (!did_converge)
            error("Newton's method did not converge.");
          
          // Translate the resulting coefficient vector into the actual solutions. 
          Solution::vector_to_solutions(coeff_vec, Tuple<Space*>(&space_T, &space_phi), 
                                        Tuple<Solution*>(&T_coarse, &phi_coarse));            
        } else {
          // Projection onto the globally derefined meshes.
          info("Projecting the latest fine mesh solution onto globally derefined meshes.");
          OGProjection::project_global(spaces, fine_mesh_solutions, coarse_mesh_solutions, matrix_solver, proj_norms); 
        }
      } 
    }
      

    // Adaptivity loop:
    
    bool done = false;
    int as = 0;
    do {
      as++;
      
      info("---- Time step %d, adaptivity step %d:", ts, as);
      
      // Visualize intermediate solutions and mesh during adaptivity.  
      view_T.show(&T_coarse);
      sprintf(title, "T (fine mesh), t = %g s, adapt step %d", TIME, as);
      view_T.set_title(title);
      
      view_phi.show(&phi_coarse);
      sprintf(title, "phi (fine mesh), t = %g s, adapt step %d", TIME, as);
      view_phi.set_title(title);
      
      ordview_T_coarse.show(&space_T);
      sprintf(title, "T mesh (coarse), t = %g, adapt step %d", TIME, as);
      ordview_T_coarse.set_title(title);
      
      ordview_phi_coarse.show(&space_phi);
      sprintf(title, "phi mesh (coarse), t = %g, adapt step %d", TIME, as);
      ordview_phi_coarse.set_title(title);

      // Construct globally refined reference meshes and setup reference spaces.
      int num_fields = 2;         // Number of physical fields being solved for (T, phi).
      int order_increase = 1;     // Increase in polynomial orders for the reference spaces.
      Tuple<Mesh *> ref_meshes;   // Reference meshes.
      Tuple<Space *> ref_spaces;  // Reference spaces.
      for (int i = 0; i < num_fields; i++) {
        ref_meshes.push_back(new Mesh());
        Mesh *ref_mesh = ref_meshes.back();
        ref_mesh->copy(spaces[i]->get_mesh());
        ref_mesh->refine_all_elements();
        ref_spaces.push_back(spaces[i]->dup(ref_mesh));
        ref_spaces[i]->copy_orders(spaces[i], order_increase);
      }     
      
      // Newton's loop on the refined meshes.
      if (as == 1) {
        info("Solving on fine meshes, starting from previous coarse mesh solutions.");
        OGProjection::project_global(ref_spaces, Tuple<MeshFunction*>((MeshFunction*)&T_coarse, (MeshFunction*)&phi_coarse), 
                       coeff_vec, matrix_solver, proj_norms);
      } else {
        info("Solving on fine meshes, starting from previous fine mesh solutions.");
        OGProjection::project_global(ref_spaces, Tuple<MeshFunction*>((MeshFunction*)&T_fine, (MeshFunction*)&phi_fine), 
                       coeff_vec, matrix_solver, proj_norms);
      }
      if( !solve_newton(ref_spaces, &wf, coeff_vec, matrix_solver, 
                        NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose) )
        error("Newton's method did not converge."); 
      
      // Translate the resulting coefficient vector into the actual solutions. 
      Solution::vector_to_solutions(coeff_vec, Tuple<Space*>(ref_spaces[0], ref_spaces[1]), 
                                    Tuple<Solution*>(&T_fine, &phi_fine));

      // Calculate error estimates and exact errors.
      info("Calculating errors.");
      double T_err_est = calc_rel_error(&T_coarse, &T_fine, HERMES_H1_NORM) * 100;
      double phi_err_est = calc_rel_error(&phi_coarse, &phi_fine, HERMES_H1_NORM) * 100;
      double T_err_exact = calc_rel_error(&T_coarse, &T_exact_solution, HERMES_H1_NORM) * 100;
      double phi_err_exact = calc_rel_error(&phi_coarse, &phi_exact_solution, HERMES_H1_NORM) * 100;
      info("T: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%, err_exact: %g %%", 
            space_T.Space::get_num_dofs(), ref_spaces[0]->Space::get_num_dofs(), T_err_est, T_err_exact);
      info("phi: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%, err_exact: %g %%", 
            space_phi.Space::get_num_dofs(), ref_spaces[1]->Space::get_num_dofs(), phi_err_est, phi_err_exact);
 
      // Calculate element errors and total error estimate for adaptivity.      
      Adapt hp(spaces, proj_norms);
      hp.set_solutions(coarse_mesh_solutions, fine_mesh_solutions);
      double err_est_rel_total = hp.calc_err_est(HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS) * 100;

/*
      if (ts==1) {
	// Add entries to DOF convergence graph.
	graph_dof_exact_T.add_values(space_T.Space::get_num_dofs(), T_err_exact);
	graph_dof_exact_T.save("conv_dof_exact_T.dat");
	graph_dof_est_T.add_values(space_T.Space::get_num_dofs(), T_err_est);
	graph_dof_est_T.save("conv_dof_est_T.dat");
	graph_dof_exact_phi.add_values(space_phi.Space::get_num_dofs(), phi_err_exact);
	graph_dof_exact_phi.save("conv_dof_exact_phi.dat");
	graph_dof_est_phi.add_values(space_phi.Space::get_num_dofs(), phi_err_est);
	graph_dof_est_phi.save("conv_dof_est_phi.dat");

	// Add entries to CPU convergence graph.
	graph_cpu_exact_T.add_values(cpu_time.accumulated(), T_err_exact);
	graph_cpu_exact_T.save("conv_cpu_exact_T.dat");
	graph_cpu_est_T.add_values(cpu_time.accumulated(), T_err_est);
	graph_cpu_est_T.save("conv_cpu_est_T.dat");
	graph_cpu_exact_phi.add_values(cpu_time.accumulated(), phi_err_exact);
	graph_cpu_exact_phi.save("conv_cpu_exact_phi.dat");
	graph_cpu_est_phi.add_values(cpu_time.accumulated(), phi_err_est);
	graph_cpu_est_phi.save("conv_cpu_est_phi.dat");
      }
*/
      
      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) done = true;
      else {
        info("Adapting the coarse meshes.");
        done = hp.adapt(Tuple<RefinementSelectors::Selector*> (&selector, &selector), THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (Space::get_num_dofs(spaces) >= NDOF_STOP) done = true; 
        
        if (!done) {
          if (SOLVE_ON_COARSE_MESH) {        
            // Newton's loop on the new coarse meshes.
            info("Solving on coarse meshes, starting from the latest fine mesh solutions.");
            OGProjection::project_global(spaces, Tuple<MeshFunction*>((MeshFunction*)&T_fine, (MeshFunction*)&phi_fine), 
                           coeff_vec, matrix_solver, proj_norms);
            did_converge = solve_newton(spaces, &wf, coeff_vec, matrix_solver, 
                                        NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose); 
            if (!did_converge)
              error("Newton's method did not converge.");
            
            // Translate the resulting coefficient vector into the actual solutions. 
            Solution::vector_to_solutions(coeff_vec, Tuple<Space*>(ref_spaces[0], ref_spaces[1]), 
                                          Tuple<Solution*>(&T_coarse, &phi_coarse));
          } else {
            // Projection onto the new coarse meshes.
            info("Projecting the latest fine mesh solution onto new coarse meshes.");
            OGProjection::project_global(spaces, fine_mesh_solutions, coarse_mesh_solutions, matrix_solver, proj_norms); 
          }
        }
      }
      
      // Free reference meshes and spaces.
      for (int i = 0; i < num_fields; i++) {
        ref_spaces[i]->free(); // This does not free the associated mesh, we must do it separately.
        ref_meshes[i]->free();
      }
    }
    while (!done);

    // Make the fine mesh solution at current time level the previous time level solution in the following time step.
    T_prev_time.copy(&T_fine);
    phi_prev_time.copy(&phi_fine);
    
    // Compute exact error.
    double T_error = calc_rel_error(&T_prev_time, &T_exact_solution, HERMES_H1_NORM) * 100;
    double phi_error = calc_rel_error(&phi_prev_time, &phi_exact_solution, HERMES_H1_NORM) * 100;
    info("Exact solution error for T (H1 norm): %g %%", T_error);
    info("Exact solution error for phi (H1 norm): %g %%", phi_error);
  }
  
  delete coeff_vec;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
