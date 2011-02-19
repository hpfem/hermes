#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Teuchos;
using namespace RefinementSelectors;

// This example solves the compressible Euler equations using the JFNK
// method implemented in the NOX package of the Trilinos library.
//
// Equations: Compressible Euler equations, perfect gas state equation.
//
// Domain: GAMM channel, see mesh file GAMM-channel.mesh
//
// BC: Normal velocity component is zero on solid walls.
//     Subsonic state prescribed on inlet and outlet.
//
// IC: Constant subsonic state identical to inlet. 
//
// The following parameters can be changed:
// Visualization.
const bool HERMES_VISUALIZATION = true;           // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_OUTPUT = true;                     // Set to "true" to enable VTK output.
const unsigned int EVERY_NTH_STEP = 1;            // Set visual output for every nth step.

// The calculation will stop when this indicator of a steady state is reached.
const double TIME_DER_LIMIT = 5E-2;

// Use of preconditioning.
const bool PRECONDITIONING = true;

const Ord2 P_INIT = Ord2(0,0);                    // Initial polynomial degree.                      
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.                       
const int INIT_REF_NUM_BOUNDARY = 1;              // Number of initial anisotropic mesh refinements towards the horizontal parts of the boundary.
double TAU = 1E-2;                                // Time step.

// Adaptivity.
const int UNREF_FREQ = 1;                         // Every UNREF_FREQth time step the mesh is unrefined.
int REFINEMENT_COUNT = 0;                         // Number of mesh refinements between two unrefinements.
                                                  // The mesh is not unrefined unless there has been a refinement since
                                                  // last unrefinement.
const double THRESHOLD = 0.1;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //  error is processed. If more elements have similar errors, refine
                                                  //  all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //  than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //  than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_H_ANISO;           // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1;                        // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.5;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                     // Adaptivity process stops when the number of degrees of freedom grows over
                                                  // this limit. This is mainly to prevent h-adaptivity to go on forever.
                                                  // Matrix solver for orthogonal projections.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Equation parameters.
const double P_EXT = 2.5;                               // Exterior pressure (dimensionless).
const double RHO_EXT = 1.0;                             // Inlet density (dimensionless).   
const double V1_EXT = 1.25;                             // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;                              // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;                               // Kappa.
// Numerical flux.
// For numerical fluxes, please see hermes2d/src/numerical_flux.h
NumericalFlux num_flux(KAPPA);

// Utility functions for the Euler equations.
#include "../euler-util.cpp"

// Calculated exterior energy.
double ENERGY_EXT = calc_energy(RHO_EXT, RHO_EXT*V1_EXT, RHO_EXT*V2_EXT, P_EXT);

// Boundary markers.
const int BDY_SOLID_WALL = 1;
const int BDY_INLET_OUTLET = 2;

// Constant initial state (matching the supersonic inlet state).
double ic_density(double x, double y, scalar& dx, scalar& dy)
{
  return RHO_EXT;
}
double ic_density_vel_x(double x, double y, scalar& dx, scalar& dy)
{
  return RHO_EXT * V1_EXT;
}
double ic_density_vel_y(double x, double y, scalar& dx, scalar& dy)
{
  return RHO_EXT * V2_EXT;
}
double ic_energy(double x, double y, scalar& dx, scalar& dy)
{
  return calc_energy(RHO_EXT, RHO_EXT*V1_EXT, RHO_EXT*V2_EXT, P_EXT);
}

// Weak forms.
#include "forms.cpp"

// Filters.
#include "filters.cpp"

// Filter for entropy which uses the constants defined above.
static void calc_entropy_estimate_func(int n, Hermes::vector<scalar*> scalars, scalar* result)
{
  for (int i = 0; i < n; i++)
    result[i] = std::log((calc_pressure(scalars.at(0)[i], scalars.at(1)[i], scalars.at(2)[i], scalars.at(3)[i]) / P_EXT)
    / pow((scalars.at(0)[i] / RHO_EXT), KAPPA));
};

// Time is zero at the beginning.
double t = 0;

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh basemesh;
  Mesh mesh;
  H2DReader mloader;
  mloader.load("GAMM-channel.mesh", &basemesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) 
    basemesh.refine_all_elements();
  basemesh.refine_towards_boundary(1, INIT_REF_NUM_BOUNDARY);
  mesh.copy(&basemesh);

  // Boundary condition types;
  BCTypes bc_types;

  // Initialize boundary condition types and spaces with default shapesets.
  bc_types.add_bc_neumann(Hermes::vector<int>(BDY_SOLID_WALL, BDY_INLET_OUTLET));
  L2Space space_rho(&mesh, &bc_types, P_INIT);
  L2Space space_rho_v_x(&mesh, &bc_types, P_INIT);
  L2Space space_rho_v_y(&mesh, &bc_types, P_INIT);
  L2Space space_e(&mesh, &bc_types, P_INIT);

  // Initialize solutions, set initial conditions.
  Solution sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e;
  Solution rsln_rho, rsln_rho_v_x, rsln_rho_v_y, rsln_e;
  sln_rho.set_exact(&mesh, ic_density);
  sln_rho_v_x.set_exact(&mesh, ic_density_vel_x);
  sln_rho_v_y.set_exact(&mesh, ic_density_vel_y);
  sln_e.set_exact(&mesh, ic_energy);
  prev_rho.set_exact(&mesh, ic_density);
  prev_rho_v_x.set_exact(&mesh, ic_density_vel_x);
  prev_rho_v_y.set_exact(&mesh, ic_density_vel_y);
  prev_e.set_exact(&mesh, ic_energy);

  // Solutions for the time derivative estimate.
  Solution sln_temp_rho, sln_temp_rho_v_x, sln_temp_rho_v_y, sln_temp_e;

  // Initialize weak formulation.
  bool is_matrix_free = true;
  WeakForm wf(4, is_matrix_free);
  
  // Volumetric linear forms.
  wf.add_vector_form(0, callback(linear_form_0_time));
  wf.add_vector_form(1, callback(linear_form_1_time));
  wf.add_vector_form(2, callback(linear_form_2_time));
  wf.add_vector_form(3, callback(linear_form_3_time));

  // Volumetric linear forms.
  // Linear forms coming from the linearization by taking the Eulerian fluxes' Jacobian matrices 
  // from the previous time step.
  // Unnecessary for FVM.
  if(P_INIT.order_h > 0 || P_INIT.order_v > 0) {
    // First flux.
    wf.add_vector_form(0, callback(linear_form_0_1), HERMES_ANY);
    
    wf.add_vector_form(1, callback(linear_form_1_0_first_flux), HERMES_ANY);
    wf.add_vector_form(1, callback(linear_form_1_1_first_flux), HERMES_ANY);
    wf.add_vector_form(1, callback(linear_form_1_2_first_flux), HERMES_ANY);
    wf.add_vector_form(1, callback(linear_form_1_3_first_flux), HERMES_ANY);
    wf.add_vector_form(2, callback(linear_form_2_0_first_flux), HERMES_ANY);
    wf.add_vector_form(2, callback(linear_form_2_1_first_flux), HERMES_ANY);
    wf.add_vector_form(2, callback(linear_form_2_2_first_flux), HERMES_ANY);
    wf.add_vector_form(2, callback(linear_form_2_3_first_flux), HERMES_ANY);
    wf.add_vector_form(3, callback(linear_form_3_0_first_flux), HERMES_ANY);
    wf.add_vector_form(3, callback(linear_form_3_1_first_flux), HERMES_ANY);
    wf.add_vector_form(3, callback(linear_form_3_2_first_flux), HERMES_ANY);
    wf.add_vector_form(3, callback(linear_form_3_3_first_flux), HERMES_ANY);

    // Second flux.
    wf.add_vector_form(0, callback(linear_form_0_2), HERMES_ANY);
    wf.add_vector_form(1, callback(linear_form_1_0_second_flux), HERMES_ANY);
    wf.add_vector_form(1, callback(linear_form_1_1_second_flux), HERMES_ANY);
    wf.add_vector_form(1, callback(linear_form_1_2_second_flux), HERMES_ANY);
    wf.add_vector_form(1, callback(linear_form_1_3_second_flux), HERMES_ANY);
    wf.add_vector_form(2, callback(linear_form_2_0_second_flux), HERMES_ANY);
    wf.add_vector_form(2, callback(linear_form_2_1_second_flux), HERMES_ANY);
    wf.add_vector_form(2, callback(linear_form_2_2_second_flux), HERMES_ANY);
    wf.add_vector_form(2, callback(linear_form_2_3_second_flux), HERMES_ANY);
    wf.add_vector_form(3, callback(linear_form_3_0_second_flux), HERMES_ANY);
    wf.add_vector_form(3, callback(linear_form_3_1_second_flux), HERMES_ANY);
    wf.add_vector_form(3, callback(linear_form_3_2_second_flux), HERMES_ANY);
    wf.add_vector_form(3, callback(linear_form_3_3_second_flux), HERMES_ANY);
  }

  // Volumetric linear forms coming from the time discretization.
  wf.add_vector_form(0, linear_form_time, linear_form_order, HERMES_ANY, &prev_rho);
  wf.add_vector_form(1, linear_form_time, linear_form_order, HERMES_ANY, &prev_rho_v_x);
  wf.add_vector_form(2, linear_form_time, linear_form_order, HERMES_ANY, &prev_rho_v_y);
  wf.add_vector_form(3, linear_form_time, linear_form_order, HERMES_ANY, &prev_e);

  // Surface linear forms - inner edges coming from the DG formulation.
  wf.add_vector_form_surf(0, linear_form_interface_0, linear_form_order, H2D_DG_INNER_EDGE);
  wf.add_vector_form_surf(1, linear_form_interface_1, linear_form_order, H2D_DG_INNER_EDGE);
  wf.add_vector_form_surf(2, linear_form_interface_2, linear_form_order, H2D_DG_INNER_EDGE);
  wf.add_vector_form_surf(3, linear_form_interface_3, linear_form_order, H2D_DG_INNER_EDGE);

  // Surface linear forms - inlet / outlet edges.
  wf.add_vector_form_surf(0, bdy_flux_inlet_outlet_comp_0, linear_form_order, BDY_INLET_OUTLET);
  wf.add_vector_form_surf(1, bdy_flux_inlet_outlet_comp_1, linear_form_order, BDY_INLET_OUTLET);
  wf.add_vector_form_surf(2, bdy_flux_inlet_outlet_comp_2, linear_form_order, BDY_INLET_OUTLET);
  wf.add_vector_form_surf(3, bdy_flux_inlet_outlet_comp_3, linear_form_order, BDY_INLET_OUTLET);

  // Surface linear forms - Solid wall edges.
  wf.add_vector_form_surf(0, bdy_flux_solid_wall_comp_0, linear_form_order, BDY_SOLID_WALL);
  wf.add_vector_form_surf(1, bdy_flux_solid_wall_comp_1, linear_form_order, BDY_SOLID_WALL);
  wf.add_vector_form_surf(2, bdy_flux_solid_wall_comp_2, linear_form_order, BDY_SOLID_WALL);
  wf.add_vector_form_surf(3, bdy_flux_solid_wall_comp_3, linear_form_order, BDY_SOLID_WALL);

  if(PRECONDITIONING) {
    // Preconditioning forms.
    wf.add_matrix_form(0, 0, callback(bilinear_form_precon));
    wf.add_matrix_form(1, 1, callback(bilinear_form_precon));
    wf.add_matrix_form(2, 2, callback(bilinear_form_precon));
    wf.add_matrix_form(3, 3, callback(bilinear_form_precon));
  }
  
  // Filters for visualization of pressure and the two components of velocity.
  SimpleFilter pressure(calc_pressure_func, Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter u(calc_u_func, Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter w(calc_w_func, Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter Mach_number(calc_Mach_func, Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter entropy_estimate(calc_entropy_estimate_func, Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));
  VectorView vview("Velocity", new WinGeom(700, 400, 600, 300));

  // Initialize refinement selector.
  L2ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Iteration number.
  int iteration = 0;
  
  // Output of the approximate time derivative.
  std::ofstream time_der_out("time_der");

  // Select preconditioner.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));

  // Approximation of the time derivative.
  // Set to a big value at the beginning.
  double time_derivative = std::numeric_limits<double>::max();

  while(time_derivative > TIME_DER_LIMIT) {
    info("---- Time step %d, time %3.5f.", iteration++, t);
    t += TAU;

    // Periodic global derefinements.
    if (iteration > 1 && iteration % UNREF_FREQ == 0 && REFINEMENT_COUNT > 0) {
      REFINEMENT_COUNT = 0;
      info("Global mesh derefinement.");
      mesh.unrefine_all_elements();
      space_rho.set_uniform_order_internal(P_INIT);
      space_rho_v_x.set_uniform_order_internal(P_INIT);
      space_rho_v_y.set_uniform_order_internal(P_INIT);
      space_e.set_uniform_order_internal(P_INIT);
    }

    // Adaptivity loop:
    int as = 1; 
    bool done = false;
    do {
      info("---- Adaptivity step %d:", as);

      // Construct globally refined reference mesh and setup reference space.
      // Global polynomial order increase = 0;
      int order_increase = 0;
      Hermes::vector<Space *>* ref_spaces = construct_refined_spaces(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), order_increase);

      // Project the previous time level solution onto the new fine mesh
      // in order to obtain initial vector for NOX. 
      info("Projecting initial solution on the FE mesh.");
      scalar* coeff_vec = new scalar[Space::get_num_dofs(*ref_spaces)];
      OGProjection::project_global(*ref_spaces, Hermes::vector<MeshFunction *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec);
      
      // Initialize the FE problem.
      bool is_linear = false;
      DiscreteProblem dp(&wf, *ref_spaces, is_linear);
  
      // If the FE problem is in fact a FV problem.
      if(P_INIT.order_h == 0 && P_INIT.order_v == 0)
        dp.set_fvm();
      
      // Initialize NOX solver.
      NoxSolver solver(&dp);
      solver.set_ls_tolerance(1E-2);
      solver.disable_abs_resid();
      solver.set_conv_rel_resid(1.00);

      if(PRECONDITIONING)
        solver.set_precond(pc);

      info("Assembling by DiscreteProblem, solving by NOX.");
      solver.set_init_sln(coeff_vec);
      if (solver.solve())
        Solution::vector_to_solutions(solver.get_solution(), *ref_spaces, 
          Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e));
      else
        error("NOX failed.");
      
      info("Number of nonlin iterations: %d (norm of residual: %g)", 
        solver.get_num_iters(), solver.get_residual());
      info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
        solver.get_num_lin_iters(), solver.get_achieved_tol());

      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection::project_global(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), 
                     Hermes::vector<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), matrix_solver, 
                     Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM)); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate.");
      Adapt* adaptivity = new Adapt(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), Hermes::vector<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      double err_est_rel_total = adaptivity->calc_err_est(Hermes::vector<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e),
							  Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e)) * 100;

      // Report results.
      info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
        Space::get_num_dofs(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e)), Space::get_num_dofs(*ref_spaces), err_est_rel_total);

      // If err_est too large, adapt the mesh.
      if (err_est_rel_total < ERR_STOP) 
        done = true;
      else {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(Hermes::vector<RefinementSelectors::Selector *>(&selector, &selector, &selector, &selector), 
                                 THRESHOLD, STRATEGY, MESH_REGULARITY);

        REFINEMENT_COUNT++;
        if (Space::get_num_dofs(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
          &space_rho_v_y, &space_e)) >= NDOF_STOP) 
          done = true;
        else
          // Increase the counter of performed adaptivity steps.
          as++;
      }

      // Clean up.
      delete adaptivity;
      if(!done)
        for(unsigned int i = 0; i < ref_spaces->size(); i++)
          delete (*ref_spaces)[i]->get_mesh();
      for(unsigned int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i];
    }
    while (done == false);

    // Approximate the time derivative of the solution.
    Adapt *adapt_for_time_der_calc = new Adapt(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
    bool solutions_for_adapt = false;
    time_derivative = iteration == 1 ? time_derivative : 
      adapt_for_time_der_calc->calc_err_est(Hermes::vector<Solution *>(&rsln_rho, &rsln_rho_v_x, &rsln_rho_v_y, &rsln_e), 
                                            Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
                                            (Hermes::vector<double>*) NULL, solutions_for_adapt) / TAU;
    delete adapt_for_time_der_calc;
    
    // Info about the approximate time derivative.
    if(iteration > 1) {
      info("Approximate the norm time derivative : %g.", time_derivative);
      time_der_out << iteration << '\t' << time_derivative << std::endl;
    }
    
    // Copy the solutions into previous time level ones.
    prev_rho.copy(&rsln_rho);
    prev_rho_v_x.copy(&rsln_rho_v_x);
    prev_rho_v_y.copy(&rsln_rho_v_y);
    prev_e.copy(&rsln_e);

    // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0) {
      // Hermes visualization.
      if(HERMES_VISUALIZATION) {
        pressure.reinit();
        u.reinit();
        w.reinit();
        Mach_number.reinit();
        entropy_estimate.reinit();
        pressure_view.show(&pressure);
        entropy_production_view.show(&entropy_estimate);
        Mach_number_view.show(&Mach_number);
        vview.show(&u, &w);
      }
      // Output solution in VTK format.
      if(VTK_OUTPUT) {
        pressure.reinit();
        Mach_number.reinit();
        Linearizer lin;
        char filename[40];
        sprintf(filename, "pressure-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&pressure, filename, "Pressure", false);
        sprintf(filename, "pressure-3D-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&pressure, filename, "Pressure", true);
        sprintf(filename, "Mach number-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&Mach_number, filename, "MachNumber", false);
        sprintf(filename, "Mach number-3D-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&Mach_number, filename, "MachNumber", true);
      }
    }
  }
  
  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();
  vview.close();

  time_der_out.close();
  return 0;
}
