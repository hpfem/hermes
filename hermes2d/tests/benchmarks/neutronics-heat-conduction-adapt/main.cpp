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
const int STRATEGY = 0;                    // Adaptive strategy:
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
const double ERR_STOP = 0.01;              // Stopping criterion for hp-adaptivity
                                           // (relative error between reference and coarse solution in percent).
const int NDOF_STOP = 10000;               // Adaptivity process stops when the number of degrees of freedom grows
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
# include "forms.cpp"

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
  int it = 1;
  while (1)
  {
    // Obtain the number of degrees of freedom.
    int ndof = Space::get_num_dofs(spaces);

    // Assemble the Jacobian matrix_coarse and residual vector.
    dp_coarse.assemble(coeff_vec_coarse, matrix_coarse, rhs_coarse, false);

    // Multiply the residual vector with -1 since the matrix_coarse 
    // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
    for (int i = 0; i < ndof; i++) rhs_coarse->set(i, -rhs_coarse->get(i));
    
    // Calculate the l2-norm of residual vector.
    double res_l2_norm = get_l2_norm(rhs_coarse);

    // Info for user.
    info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(spaces), res_l2_norm);

    // If l2 norm of the residual vector is within tolerance, or the maximum number 
    // of iteration has been reached, then quit.
    if (res_l2_norm < NEWTON_TOL_COARSE || it > NEWTON_MAX_ITER) break;

    // Solve the linear system.
    if(!solver_coarse->solve())
      error ("matrix_coarse solver_coarse failed.\n");

    // Add \deltaY^{n+1} to Y^n.
    for (int i = 0; i < ndof; i++) coeff_vec_coarse[i] += solver_coarse->get_solution()[i];
    
    if (it >= NEWTON_MAX_ITER)
      error ("Newton method did not converge.");

    it++;
  }
  
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
          info("Solving on globally derefined meshes, starting from the latest fine mesh solutions.");
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
          int it = 1;
          while (1)
          {
            // Obtain the number of degrees of freedom.
            int ndof = Space::get_num_dofs(spaces);

            // Assemble the Jacobian matrix_coarse and residual vector.
            dp_coarse.assemble(coeff_vec_coarse, matrix_coarse, rhs_coarse, false);

            // Multiply the residual vector with -1 since the matrix_coarse 
            // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
            for (int i = 0; i < ndof; i++) rhs_coarse->set(i, -rhs_coarse->get(i));
            
            // Calculate the l2-norm of residual vector.
            double res_l2_norm = get_l2_norm(rhs_coarse);

            // Info for user.
            info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(spaces), res_l2_norm);

            // If l2 norm of the residual vector is within tolerance, or the maximum number 
            // of iteration has been reached, then quit.
            if (res_l2_norm < NEWTON_TOL_COARSE || it > NEWTON_MAX_ITER) break;

            // Solve the linear system.
            if(!solver_coarse->solve())
              error ("matrix_coarse solver_coarse failed.\n");

            // Add \deltaY^{n+1} to Y^n.
            for (int i = 0; i < ndof; i++) coeff_vec_coarse[i] += solver_coarse->get_solution()[i];
            
            if (it >= NEWTON_MAX_ITER)
              error ("Newton method did not converge.");

            it++;
          }
          
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

      // Construct globally refined reference mesh
      // and setup reference space.
      Hermes::Tuple<Space *>* ref_spaces = construct_refined_spaces(spaces);

      // Newton's loop on the refined meshes.
      scalar* coeff_vec = new scalar[Space::get_num_dofs(*ref_spaces)];
      if (as == 1) {
        info("Solving on fine meshes, starting from previous coarse mesh solutions.");
        OGProjection::project_global(*ref_spaces, Hermes::Tuple<MeshFunction*>((MeshFunction*)&T_coarse, (MeshFunction*)&phi_coarse), 
                       coeff_vec, matrix_solver, proj_norms);
      } else {
        info("Solving on fine meshes, starting from previous fine mesh solutions.");
        OGProjection::project_global(*ref_spaces, Hermes::Tuple<MeshFunction*>((MeshFunction*)&T_fine, (MeshFunction*)&phi_fine), 
                       coeff_vec, matrix_solver, proj_norms);
      }
      
      // Initialize the FE problem.
      bool is_linear = false;
      DiscreteProblem dp(&wf, *ref_spaces, is_linear);

      // Set up the solver, matrix, and rhs according to the solver selection.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

      // Perform Newton's iteration.
      int it = 1;
      while (1)
      {
        // Obtain the number of degrees of freedom.
        int ndof = Space::get_num_dofs(*ref_spaces);

        // Assemble the Jacobian matrix and residual vector.
        dp.assemble(coeff_vec, matrix, rhs, false);

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
            int it = 1;
            while (1)
            {
              // Obtain the number of degrees of freedom.
              int ndof = Space::get_num_dofs(spaces);

              // Assemble the Jacobian matrix_coarse and residual vector.
              dp.assemble(coeff_vec_coarse, matrix_coarse, rhs_coarse, false);

              // Multiply the residual vector with -1 since the matrix_coarse 
              // equation reads J(Y^n) \deltaY^{n+1} = -F(Y^n).
              for (int i = 0; i < ndof; i++) rhs_coarse->set(i, -rhs_coarse->get(i));
              
              // Calculate the l2-norm of residual vector.
              double res_l2_norm = get_l2_norm(rhs_coarse);

              // Info for user.
              info("---- Newton iter %d, ndof %d, res. l2 norm %g", it, Space::get_num_dofs(spaces), res_l2_norm);

              // If l2 norm of the residual vector is within tolerance, or the maximum number 
              // of iteration has been reached, then quit.
              if (res_l2_norm < NEWTON_TOL_COARSE || it > NEWTON_MAX_ITER) break;

              // Solve the linear system.
              if(!solver_coarse->solve())
                error ("matrix_coarse solver_coarse failed.\n");

              // Add \deltaY^{n+1} to Y^n.
              for (int i = 0; i < ndof; i++) coeff_vec_coarse[i] += solver_coarse->get_solution()[i];
              
              if (it >= NEWTON_MAX_ITER)
                error ("Newton method did not converge.");

              it++;
            }
            
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
      for(int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i]->get_mesh();
      delete ref_spaces;
    }
    while (!done);

    // Make the fine mesh solution at current time level the previous time level solution in the following time step.
    T_prev_time.copy(&T_fine);
    phi_prev_time.copy(&phi_fine);
  }

  info("Coordinate ( 25.0,  25.0) T value = %lf", T_prev_time.get_pt_value(25.0, 25.0));
  info("Coordinate ( 25.0,  75.0) T value = %lf", T_prev_time.get_pt_value(25.0, 75.0));
  info("Coordinate ( 75.0, 25.0) T value = %lf", T_prev_time.get_pt_value(75.0, 25.0));
  info("Coordinate ( 75.0, 75.0) T value = %lf", T_prev_time.get_pt_value(75.0, 75.0));
  info("Coordinate ( 50.0, 50.0) T value = %lf", T_prev_time.get_pt_value(50.0, 50.0));

  info("Coordinate ( 25.0, 25.0) phi value = %lf", phi_prev_time.get_pt_value(25.0, 25.0));
  info("Coordinate ( 25.0, 75.0) phi value = %lf", phi_prev_time.get_pt_value(25.0, 75.0));
  info("Coordinate ( 75.0, 25.0) phi value = %lf", phi_prev_time.get_pt_value(75.0, 25.0));
  info("Coordinate ( 75.0, 75.0) phi value = %lf", phi_prev_time.get_pt_value(75.0, 75.0));
  info("Coordinate ( 50.0, 50.0) phi value = %lf", phi_prev_time.get_pt_value(50.0, 50.0));

  int ndof_allowed_T = 130;
  int ndof_allowed_phi = 300;
  int ndof_T = Space::get_num_dofs(&space_T);
  int ndof_phi = Space::get_num_dofs(&space_phi);
  printf("ndof_actual_T = %d\n", ndof_T);
  printf("ndof_actual_phi = %d\n", ndof_phi);
  printf("ndof_allowed_T = %d\n", ndof_allowed_T);
  printf("ndof_allowed_phi = %d\n", ndof_allowed_phi);
  if ((ndof_T <= ndof_allowed_T) && (ndof_phi <= ndof_allowed_phi)) {  
    // ndofs_T was 121 and ndof_phi was 267 at the time this test was created
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
