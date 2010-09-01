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

// Time-stepping:
const double TAU = 0.1;                    // Time step.
const double T_FINAL = 10.0;               // Time interval length.

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
const double ERR_STOP = 0.01;              // Stopping criterion for hp-adaptivity
                                           // (relative error between reference and coarse solution in percent).
const int NDOF_STOP = 60000;               // Adaptivity process stops when the number of degrees of freedom grows
                                           // over this limit. This is to prevent h-adaptivity to go on forever.


// Newton's method:
const double NEWTON_TOL_COARSE = 1.0e-2;          // Stopping criterion for Newton on coarse mesh.
const double NEWTON_TOL_FINE = 5.0e-2;            // Stopping criterion for Newton on fine mesh.
const int NEWTON_MAX_ITER = 20;                   // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC, SOLVER_MUMPS
                                                  // and more are coming.

// Problem parameters.
const double CT = 1.0;
const double CF = 1.0;
const double rT = 1.0;
const double rF = 0.25;
const double LX = 100.0;          // Domain sizes in the x and y dimensions.
const double LY = 100.0;
const double invvel = 2.0e-4;     // Inverse of neutron velocity.
const double xsdiff = 1.268;      // Diffusion coefficient.
const double Tref = 0.0;          // Temperature at boundary.

const double nu = 2.41;           // Number of neutrons emitted per fission event.
const double xsfiss = 0.00191244; // Fission cross section.
const double kappa = 1.0e-6;
const double rho = 1.0;           // Density.
const double cp = 1.0;            // Heat capacity.

const double PI = acos(-1.0);
const double normalization_const = 1.0;

const double energy_per_fission = kappa * xsfiss;

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
  return xsa_ref + doppler_coeff * (sqrt(T+1.0e-10) - sqrt(Tref));
  //return xsa_ref + doppler_coeff * (sqrt(T) - sqrt(Tref));
}

// Derivative of the removal cross section with respect to temperature
template<typename Real>
Real dxsrem_dT(Real T) {
  return doppler_coeff / (2*sqrt(T+1.0e-10));
  //return doppler_coeff / (2*sqrt(T));
}

// Heat source.
template<typename Real>
Real qT(Real x, Real y) {
  return
rho*cp*CT*(1.0-pow(tanh(rT*TIME),2.0))*rT*sin(x/LX*PI)*sin(y/LY*PI)-k1*CT*CT*pow(1.0+tanh(rT*TIME),2.0)*pow(cos(x/LX*PI),2.0)/(LX*LX)*PI*PI*pow(sin(y/LY*PI),2.0)+(k0+k1*(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)-Tref))*CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)/(LX*LX)*PI*PI*sin(y/LY*PI)-k1*CT*CT*pow(1.0+tanh(rT*TIME),2.0)*pow(sin(x/LX*PI),2.0)*pow(cos(y/LY*PI),2.0)/(LY*LY)*PI*PI+(k0+k1*(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)-Tref))*CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI)/(LY*LY)*PI*PI-normalization_const*energy_per_fission*xsfiss*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY;
}

// Extraneous neutron source.
template<typename Real>
Real q(Real x, Real y) {
  return 
invvel*CF*rF*exp(rF*TIME)*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY-xsdiff*(-CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)/(LX*LX*LX)*PI*PI*sin(y/LY*PI)*x*y/LY+2.0*CF*(1.0+exp(rF*TIME))*cos(x/LX*PI)/(LX*LX)*PI*sin(y/LY*PI)*y/LY)-xsdiff*(-CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)/(LY*LY*LY)*PI*PI*x/LX*y+2.0*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*cos(y/LY*PI)/(LY*LY)*PI*x/LX)+(xsa_ref+doppler_coeff*(sqrt(CT*(1.0+tanh(rT*TIME))*sin(x/LX*PI)*sin(y/LY*PI))-sqrt(Tref)))*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY-nu*xsfiss*CF*(1.0+exp(rF*TIME))*sin(x/LX*PI)*sin(y/LY*PI)*x/LX*y/LY;
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
  return Tref;
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
  // Load the mesh.
  Mesh basemesh, mesh_T, mesh_phi;
  H2DReader mloader;
  mloader.load("domain.mesh", &basemesh);

  // Perform initial mesh refinements.
  for (int i=0; i < INIT_GLOB_REF_NUM; i++)  basemesh.refine_all_elements();
  basemesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);
  
  // Create a special mesh for each physical field.
  mesh_T.copy(&basemesh);
  mesh_phi.copy(&basemesh);

  // Create H1 spaces with default shapesets.
  H1Space space_T(&mesh_T, bc_types_T, essential_bc_values_T, P_INIT);
  H1Space space_phi(&mesh_phi, bc_types_phi, essential_bc_values_phi, P_INIT);
  Tuple<Space*> spaces(&space_T, &space_phi);
  int ndof = get_num_dofs(spaces); 
 
  // Solutions in the previous time step (converging within the time stepping loop).
  Solution T_prev_time, phi_prev_time;
  Tuple<MeshFunction*> prev_time_meshfns(&T_prev_time, &phi_prev_time);
  
  // Solutions on the coarse and refined meshes in current time step (converging within the Newton's loop).
  Solution T_coarse, phi_coarse, T_fine, phi_fine;
  Tuple<MeshFunction*> coarse_mesh_meshfns(&T_coarse, &phi_coarse);
  Tuple<MeshFunction*> fine_mesh_meshfns(&T_fine, &phi_fine);
  Tuple<Solution*> coarse_mesh_solutions(&T_coarse, &phi_coarse);
  Tuple<Solution*> fine_mesh_solutions(&T_fine, &phi_fine);
  
  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, jac_TT, jac_TT_ord);
  wf.add_matrix_form(0, 1, jac_Tphi, jac_Tphi_ord);
  wf.add_vector_form(0, res_T, res_T_ord, H2D_ANY, &T_prev_time);
  wf.add_matrix_form(1, 0, jac_phiT, jac_phiT_ord);
  wf.add_matrix_form(1, 1, jac_phiphi, jac_phiphi_ord);
  wf.add_vector_form(1, res_phi, res_phi_ord, H2D_ANY, &phi_prev_time);

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
  
  // Exact solutions for error evaluation.
  ExactSolution T_exact_solution(&mesh_T, T_exact);
  ExactSolution phi_exact_solution(&mesh_phi, phi_exact);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize the nonlinear system.
  DiscreteProblem dp(&wf, spaces);
  Tuple<int> proj_norms(H2D_H1_NORM, H2D_H1_NORM);
  
  // Set initial conditions.
  T_prev_time.set_exact(&mesh_T, T_exact);
  phi_prev_time.set_exact(&mesh_phi, phi_exact);
  
  // Newton's loop on the initial coarse meshes.
  info("Solving on coarse meshes.");
  Vector* coeff_vec = new AVector();
  project_global(spaces, proj_norms, prev_time_meshfns, Tuple<Solution*>(), coeff_vec);
  bool verbose = true; // Default is false.
  bool did_converge = solve_newton(spaces, &wf, coeff_vec, matrix_solver, 
                                   NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose); 
  if (!did_converge)
    error("Newton's method did not converge.");

  // Translate the resulting coefficient vector into the actual solutions. 
  T_coarse.set_fe_solution(&space_T, coeff_vec);
  phi_coarse.set_fe_solution(&space_phi, coeff_vec);
  
  
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
          project_global(spaces, proj_norms, fine_mesh_meshfns, Tuple<Solution*>(), coeff_vec);
          did_converge = solve_newton(spaces, &wf, coeff_vec, matrix_solver, 
                                      NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose); 
          if (!did_converge)
            error("Newton's method did not converge.");
          
          // Translate the resulting coefficient vector into the actual solutions. 
          T_coarse.set_fe_solution(&space_T, coeff_vec);
          phi_coarse.set_fe_solution(&space_phi, coeff_vec);
          
        } else {
          // Projection onto the globally derefined meshes.
          info("Projecting the latest fine mesh solution onto globally derefined meshes.");
          project_global(spaces, proj_norms, fine_mesh_meshfns, coarse_mesh_solutions); 
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
        project_global(ref_spaces, proj_norms, coarse_mesh_meshfns, Tuple<Solution*>(), coeff_vec);
      } else {
        info("Solving on fine meshes, starting from previous fine mesh solutions.");
        project_global(ref_spaces, proj_norms, fine_mesh_meshfns, Tuple<Solution*>(), coeff_vec);
      }
      if( !solve_newton(ref_spaces, &wf, coeff_vec, matrix_solver, 
                        NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose) )
        error("Newton's method did not converge."); 
      
      // Translate the resulting coefficient vector into the actual solutions. 
      T_fine.set_fe_solution(ref_spaces[0], coeff_vec);
      phi_fine.set_fe_solution(ref_spaces[1], coeff_vec);

      // Calculate error estimates and exact errors.
      info("Calculating errors.");
      double T_err_est = calc_rel_error(&T_coarse, &T_fine, H2D_H1_NORM) * 100;
      double phi_err_est = calc_rel_error(&phi_coarse, &phi_fine, H2D_H1_NORM) * 100;
      double T_err_exact = calc_rel_error(&T_coarse, &T_exact_solution, H2D_H1_NORM) * 100;
      double phi_err_exact = calc_rel_error(&phi_coarse, &phi_exact_solution, H2D_H1_NORM) * 100;
      info("T: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%, err_exact: %g %%", 
            space_T.get_num_dofs(), ref_spaces[0]->get_num_dofs(), T_err_est, T_err_exact);
      info("phi: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%, err_exact: %g %%", 
            space_phi.get_num_dofs(), ref_spaces[1]->get_num_dofs(), phi_err_est, phi_err_exact);
 
      // Calculate element errors and total error estimate for adaptivity.      
      Adapt hp(spaces, proj_norms);
      hp.set_solutions(coarse_mesh_solutions, fine_mesh_solutions);
      hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_ABS) * 100;

      double err_est = 0.0, norm_est = 0.0;
      for (int i = 0; i < num_fields; i++) {
        double cur_field_err_est = calc_abs_error( coarse_mesh_solutions[i], fine_mesh_solutions[i], proj_norms[i] );
        double cur_field_norm_est = calc_norm( fine_mesh_solutions[i], proj_norms[i] );
        err_est += sqr(cur_field_err_est);
        norm_est += sqr(cur_field_norm_est);
      }

      err_est = sqrt(err_est/norm_est) * 100.;  // Get the percentual relative error estimate.
      
      // If err_est too large, adapt the mesh.
      if (err_est < ERR_STOP) done = true;
      else {
        info("Adapting the coarse meshes.");
        done = hp.adapt(Tuple<RefinementSelectors::Selector*> (&selector, &selector), THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (get_num_dofs(spaces) >= NDOF_STOP) done = true; 
        
        if (!done) {
          if (SOLVE_ON_COARSE_MESH) {        
            // Newton's loop on the new coarse meshes.
            info("Solving on coarse meshes, starting from the latest fine mesh solutions.");
            project_global(spaces, proj_norms, fine_mesh_meshfns, Tuple<Solution*>(), coeff_vec);
            did_converge = solve_newton(spaces, &wf, coeff_vec, matrix_solver, 
                                        NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose); 
            if (!did_converge)
              error("Newton's method did not converge.");
            
            // Translate the resulting coefficient vector into the actual solutions. 
            T_coarse.set_fe_solution(&space_T, coeff_vec);
            phi_coarse.set_fe_solution(&space_phi, coeff_vec);
          } else {
            // Projection onto the new coarse meshes.
            info("Projecting the latest fine mesh solution onto new coarse meshes.");
            project_global(spaces, proj_norms, fine_mesh_meshfns, coarse_mesh_solutions, NULL); 
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
    double T_error = calc_rel_error(&T_prev_time, &T_exact_solution, H2D_H1_NORM) * 100;
    double phi_error = calc_rel_error(&phi_prev_time, &phi_exact_solution, H2D_H1_NORM) * 100;
    info("Exact solution error for T (H1 norm): %g %%", T_error);
    info("Exact solution error for phi (H1 norm): %g %%", phi_error);
  }
  
  delete coeff_vec;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
