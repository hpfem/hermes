#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
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
const double T_FINAL = 10*TAU;             // Time interval length.

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
//  return 1.0;
  return 1+tanh(rT*TIME);
}

template<typename Real>
Real DT_FTIME(Real x, Real y) {
//  return 0.0;
  return rT*(1-pow(tanh(rT*TIME),2));
}

// Time dependence of the neutron flux.
template<typename Real>
Real PHI_FTIME(Real x, Real y) {
//  return T_FTIME(x, y);
  return 1+exp(rF*TIME);
}

template<typename Real>
Real DPHI_FTIME(Real x, Real y) {
//  return DT_FTIME(x, y);
  return rF*exp(rF*TIME);
}

// Heat source.
template<typename Real>
Real qT(Real x, Real y) {
  return CT*DT_FTIME(x,y)*cp*rho*sin((M_PI*x)/LX)*sin((M_PI*y)/LY)+CT*1/(LX*LX)*(M_PI*M_PI)*T_FTIME(x,y)*sin((M_PI*x)/LX)*sin((M_PI*y)/LY)*(k0-k1*(Tref-CT*T_FTIME(x,y)*sin((M_PI*x)/LX)*sin((M_PI*y)/LY)))+CT*1/(LY*LY)*(M_PI*M_PI)*T_FTIME(x,y)*sin((M_PI*x)/LX)*sin((M_PI*y)/LY)*(k0-k1*(Tref-CT*T_FTIME(x,y)*sin((M_PI*x)/LX)*sin((M_PI*y)/LY)))-(CT*CT)*1/(LX*LX)*(M_PI*M_PI)*(T_FTIME(x,y)*T_FTIME(x,y))*k1*pow(cos((M_PI*x)/LX),2.0)*pow(sin((M_PI*y)/LY),2.0)-(CT*CT)*1/(LY*LY)*(M_PI*M_PI)*(T_FTIME(x,y)*T_FTIME(x,y))*k1*pow(cos((M_PI*y)/LY),2.0)*pow(sin((M_PI*x)/LX),2.0)-(CF*PHI_FTIME(x,y)*kappa*x*xsfiss*y*sin((M_PI*x)/LX)*sin((M_PI*y)/LY))/(LX*LY);
}

// Extraneous neutron source.
template<typename Real>
Real q(Real x, Real y) {
  return -xsdiff*((CF*1/(LY*LY)*PHI_FTIME(x,y)*M_PI*x*cos((M_PI*y)/LY)*sin((M_PI*x)/LX)*2.0)/LX-(CF*1/(LY*LY*LY)*PHI_FTIME(x,y)*(M_PI*M_PI)*x*y*sin((M_PI*x)/LX)*sin((M_PI*y)/LY))/LX)-xsdiff*((CF*1/(LX*LX)*PHI_FTIME(x,y)*M_PI*y*cos((M_PI*x)/LX)*sin((M_PI*y)/LY)*2.0)/LY-(CF*1/(LX*LX*LX)*PHI_FTIME(x,y)*(M_PI*M_PI)*x*y*sin((M_PI*x)/LX)*sin((M_PI*y)/LY))/LY)+(CF*DPHI_FTIME(x,y)*invvel*x*y*sin((M_PI*x)/LX)*sin((M_PI*y)/LY))/(LX*LY)+(CF*PHI_FTIME(x,y)*x*y*sin((M_PI*x)/LX)*sin((M_PI*y)/LY)*(xsa_ref-doppler_coeff*(sqrt(Tref)-sqrt(CT*T_FTIME(x,y)*sin((M_PI*x)/LX)*sin((M_PI*y)/LY)))))/(LX*LY)-(CF*PHI_FTIME(x,y)*nu*x*xsfiss*y*sin((M_PI*x)/LX)*sin((M_PI*y)/LY))/(LX*LY);
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
  Hermes::Tuple<Space*> spaces(&space_T, &space_phi);
  int ndof = Space::get_num_dofs(spaces); 
 
  // Solutions in the previous time step (converging within the time stepping loop).
  Solution T_prev_time, phi_prev_time;
  Hermes::Tuple<Solution*> prev_time_solutions(&T_prev_time, &phi_prev_time);
  
  // Solutions on the coarse and refined meshes in current time step (converging within the Newton's loop).
  Solution T_coarse, phi_coarse, T_fine, phi_fine;
  Hermes::Tuple<Solution*> coarse_mesh_solutions(&T_coarse, &phi_coarse);
  Hermes::Tuple<Solution*> fine_mesh_solutions(&T_fine, &phi_fine);
  
  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, jac_TT, jac_TT_ord);
  wf.add_matrix_form(0, 1, jac_Tphi, jac_Tphi_ord);
  wf.add_vector_form(0, res_T, res_T_ord, HERMES_ANY, &T_prev_time);
  wf.add_matrix_form(1, 0, jac_phiT, jac_phiT_ord);
  wf.add_matrix_form(1, 1, jac_phiphi, jac_phiphi_ord);
  wf.add_vector_form(1, res_phi, res_phi_ord, HERMES_ANY, &phi_prev_time);

  // Initialize solution views (their titles will be updated in each time step).
  ScalarView view_T("", new WinGeom(460, 0, 450, 350));
  view_T.fix_scale_width(80);
  ScalarView view_T_exact("", new WinGeom(0, 0, 450, 350));
  view_T_exact.fix_scale_width(80);
  view_T_exact.show_mesh(false);
  ScalarView view_phi("", new WinGeom(460, 400, 450, 350));
  view_phi.fix_scale_width(80);
  ScalarView view_phi_exact("", new WinGeom(0, 400, 450, 350));
  view_phi_exact.fix_scale_width(80);
  view_phi_exact.show_mesh(false);

  // Initialize mesh views (their titles will be updated in each time step).
  OrderView ordview_T_coarse("", new WinGeom(920, 0, 450, 350));
  ordview_T_coarse.fix_scale_width(80);
  OrderView ordview_phi_coarse("", new WinGeom(920, 400, 450, 350));
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
  DiscreteProblem dp(&wf, spaces);
  Hermes::Tuple<ProjNormType> proj_norms(HERMES_H1_NORM, HERMES_H1_NORM);
  
  // Set initial conditions.
  T_prev_time.set_exact(&mesh_T, T_exact);
  phi_prev_time.set_exact(&mesh_phi, phi_exact);
  
  // Newton's loop on the initial coarse meshes.
  info("Solving on coarse meshes.");
  scalar* coeff_vec_coarse = new scalar[Space::get_num_dofs(spaces)];
  OGProjection::project_global(spaces, Hermes::Tuple<MeshFunction*>((MeshFunction*)&T_prev_time, (MeshFunction*)&phi_prev_time), 
                 coeff_vec_coarse, matrix_solver, proj_norms);
  bool verbose = true; // Default is false.
  
  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp_coarse(&wf, spaces, is_linear);

  // Set up the solver_coarse, matrix_coarse, and rhs_coarse according to the solver_coarse selection.
  SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
  Vector* rhs_coarse = create_vector(matrix_solver);
  Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

  // Perform Newton's iteration.
  if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) 
    error("Newton's iteration failed.");
  
  // Translate the resulting coefficient vector into the actual solutions. 
  Solution::vector_to_solutions(coeff_vec_coarse, spaces, coarse_mesh_solutions);

  delete [] coeff_vec_coarse;
  delete rhs_coarse;
  delete matrix_coarse;
  delete solver_coarse;
  
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
          scalar* coeff_vec_coarse = new scalar[Space::get_num_dofs(spaces)];
          info("Projecting previous fine mesh solution to obtain initial vector for Newton's iteration on globally derefined meshes.");
          OGProjection::project_global(spaces, Hermes::Tuple<MeshFunction*>((MeshFunction*)&T_fine, (MeshFunction*)&phi_fine), 
                         coeff_vec_coarse, matrix_solver, proj_norms);
          
          // Initialize the FE problem.
          bool is_linear = false;
          DiscreteProblem dp_coarse(&wf, spaces, is_linear);

          // Set up the solver_coarse, matrix_coarse, and rhs_coarse according to the solver_coarse selection.
          SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
          Vector* rhs_coarse = create_vector(matrix_solver);
          Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

          // Perform Newton's iteration.
          if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) 
            error("Newton's iteration failed.");
          
          // Translate the resulting coefficient vector into the actual solutions. 
          Solution::vector_to_solutions(coeff_vec_coarse, spaces, coarse_mesh_solutions);
          delete [] coeff_vec_coarse;
          delete rhs_coarse;
          delete matrix_coarse;
          delete solver_coarse;
        } 
        else {
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

      // Construct globally refined reference mesh
      // and setup reference space.
      Hermes::Tuple<Space *>* ref_spaces = construct_refined_spaces(spaces);

      // Newton's loop on the refined meshes.
      scalar* coeff_vec = new scalar[Space::get_num_dofs(*ref_spaces)];
      if (as == 1) {
        info("Projecting coarse mesh solution to obtain coefficients vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces, Hermes::Tuple<MeshFunction*>((MeshFunction*)&T_coarse, (MeshFunction*)&phi_coarse), 
                       coeff_vec, matrix_solver, proj_norms);
      } else {
        info("Projecting previous fine mesh solution to obtain coefficients vector on new fine mesh.");
        OGProjection::project_global(*ref_spaces, Hermes::Tuple<MeshFunction*>((MeshFunction*)&T_fine, (MeshFunction*)&phi_fine), 
                       coeff_vec, matrix_solver, proj_norms);
        
        // Deallocate the previous fine mesh.
        delete T_fine.get_mesh();
        delete phi_fine.get_mesh();
      }
      
      // Initialize the FE problem.
      bool is_linear = false;
      DiscreteProblem dp(&wf, *ref_spaces, is_linear);

      // Set up the solver, matrix, and rhs according to the solver selection.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Perform Newton's iteration.
      if (!solve_newton(coeff_vec, &dp, solver, matrix, rhs, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose)) 
        error("Newton's iteration failed.");
            
      // Translate the resulting coefficient vector into the actual solutions. 
      Solution::vector_to_solutions(coeff_vec, *ref_spaces, fine_mesh_solutions);
      delete [] coeff_vec;
      delete rhs;
      delete matrix;
      delete solver;

      // Calculate element errors.
      info("Calculating error estimate and exact error."); 
      Adapt* adaptivity = new Adapt(spaces, proj_norms);

      // Calculate error estimate for each solution component and the total error estimate.
      bool solutions_for_adapt = true;
      Hermes::Tuple<double> err_est_rel;
      double err_est_rel_total = adaptivity->calc_err_est(coarse_mesh_solutions, fine_mesh_solutions, solutions_for_adapt, 
                                 HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS, &err_est_rel) * 100;

      // Calculate exact error for each solution component and the total exact error.
      solutions_for_adapt = false;
      Hermes::Tuple<double> err_exact_rel;
      double err_exact_rel_total = adaptivity->calc_err_exact(coarse_mesh_solutions, Hermes::Tuple<Solution *>(&T_exact_solution, &phi_exact_solution), solutions_for_adapt, 
                                 HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS, &err_exact_rel) * 100;

      info("T: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%, err_exact: %g %%", 
            space_T.get_num_dofs(), (*ref_spaces)[0]->get_num_dofs(), err_est_rel[0]*100, err_exact_rel[0]*100);
      info("phi: ndof_coarse: %d, ndof_fine: %d, err_est: %g %%, err_exact: %g %%", 
            space_phi.get_num_dofs(), (*ref_spaces)[1]->get_num_dofs(), err_est_rel[1]*100, err_exact_rel[1]*100);
 
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
        done = adaptivity->adapt(Hermes::Tuple<RefinementSelectors::Selector*> (&selector, &selector), THRESHOLD, STRATEGY, MESH_REGULARITY);
        if (Space::get_num_dofs(spaces) >= NDOF_STOP) done = true; 
        
        if (!done) {
          if (SOLVE_ON_COARSE_MESH) {        
            // Newton's loop on the new coarse meshes.
            scalar* coeff_vec_coarse = new scalar[Space::get_num_dofs(spaces)];
            info("Solving on coarse meshes, starting from the latest fine mesh solutions.");
            OGProjection::project_global(spaces, Hermes::Tuple<MeshFunction*>((MeshFunction*)&T_fine, (MeshFunction*)&phi_fine), 
                           coeff_vec_coarse, matrix_solver, proj_norms);
            
            // Initialize the FE problem.
            bool is_linear = false;
            DiscreteProblem dp_coarse(&wf, spaces, is_linear);

            // Set up the solver_coarse, matrix_coarse, and rhs_coarse according to the solver_coarse selection.
            SparseMatrix* matrix_coarse = create_matrix(matrix_solver);
            Vector* rhs_coarse = create_vector(matrix_solver);
            Solver* solver_coarse = create_linear_solver(matrix_solver, matrix_coarse, rhs_coarse);

            // Perform Newton's iteration.
            if (!solve_newton(coeff_vec_coarse, &dp_coarse, solver_coarse, matrix_coarse, rhs_coarse, NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) 
              error("Newton's iteration failed.");
            
            // Translate the resulting coefficient vector into the actual solutions. 
            Solution::vector_to_solutions(coeff_vec_coarse, spaces, coarse_mesh_solutions);
            delete [] coeff_vec_coarse;
            delete rhs_coarse;
            delete matrix_coarse;
            delete solver_coarse;
          } 
          else {
            // Projection onto the new coarse meshes.
            info("Projecting the latest fine mesh solution onto new coarse meshes.");
            OGProjection::project_global(spaces, fine_mesh_solutions, coarse_mesh_solutions, matrix_solver, proj_norms); 
          }
        }
      }
      delete adaptivity;
      delete ref_spaces;
    }
    while (!done);

    // Make the fine mesh solution at current time level the previous time level solution in the following time step.
    T_prev_time.copy(&T_fine);
    phi_prev_time.copy(&phi_fine);
  }
  
  // Wait for all views to be closed.
  View::wait();
  return 0;
}
