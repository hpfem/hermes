#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Teuchos;
using namespace RefinementSelectors;

// This example solves the compressible Euler equations coupled with an advection-diffution equation
// using a basic piecewise-constant finite volume method for the flow and continuous FEM for the concentration
// being advected by the flow.
//
// Equations: Compressible Euler equations, perfect gas state equation, advection-diffusion equation.
//
// Domains: GAMM channel, see mesh file GAMM-channel-4-bnds.mesh
//          a rectangular channel, see mesh file channel-4-bnds.mesh
//
// BC: Normal velocity component is zero on solid walls.
//     Subsonic state prescribed on inlet and outlet.
//     Various conditions for the concentration.
//     See the parameter INITIAL_CONCENTRATION_STATE.
//
// IC: Constant subsonic state identical to inlet. 
//     Various conditions for the concentration.
//     See the parameter INITIAL_CONCENTRATION_STATE.
//
// The following parameters can be changed.
// Some of them are not constants to be able to change them via command line arguments:

// This is the main parameter.
// Meaning: 0 - concentration is kept constant at the bottom of the domain.
//              at the beginning, concentration is equal throughout the whole
//              domain and is equal to the value at the bottom.
//          1 - concentration is kept constant at the bottom of the domain.
//              at the beginning, concentration is zero throughout the domain.
//          2 - concentration is kept constant at the inlet part of the domain.
//              at the beginning, concentration is zero throughout the domain.
// If not said otherwise, zero Neumann condition is imposed on all parts of the boundary.
unsigned int INITIAL_CONCENTRATION_STATE = 1;

// Use of preconditioning.
const bool PRECONDITIONING = true;

// Visualization.
const bool HERMES_VISUALIZATION = true;           // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_OUTPUT = false;                    // Set to "true" to enable VTK output.
const unsigned int EVERY_NTH_STEP = 1;            // Set visual output for every nth step.

const Ord2 P_INIT_FLOW = Ord2(0,0);               // Polynomial degree for the Euler equations (for the flow).
const Ord2 P_INIT_CONCENTRATION = Ord2(1,1);      // Polynomial degree for the concentration.
double TAU = 5E-2;                                // Time step.
const MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

unsigned int INIT_REF_NUM_FLOW = 2;               // Number of initial uniform mesh refinements of the mesh for the flow.
unsigned int INIT_REF_NUM_CONCENTRATION = 2;      // Number of initial uniform mesh refinements of the mesh for the concentration.
unsigned int INIT_REF_NUM_CONCENTRATION_BDY = 5;  // Number of initial mesh refinements of the mesh for the concentration towards the 
                                                  // part of the boundary where the concentration is prescribed.


// Equation parameters.
const double P_EXT = 2.5;         // Exterior pressure (dimensionless).
const double RHO_EXT = 1.0;       // Inlet density (dimensionless).   
const double V1_EXT = 1.25;       // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;        // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;         // Kappa.
// Numerical flux.
// For numerical fluxes, please see hermes2d/src/numerical_flux.h
NumericalFlux num_flux(KAPPA);

// Utility functions for the Euler equations.
#include "../euler-util.cpp"

// Calculated exterior energy.
double ENERGY_EXT = calc_energy(RHO_EXT, RHO_EXT*V1_EXT, RHO_EXT*V2_EXT, P_EXT);

// Diffusion parameter (diffusivity).
const double EPSILON = 1E-3;

// Boundary (initial) value of the concentration.
const double CONCENTRATION_EXT = 1.0;

// Boundary markers.
const int BDY_INLET = 1;
const int BDY_OUTLET = 2;
const int BDY_SOLID_WALL_BOTTOM = 3;
const int BDY_SOLID_WALL_TOP = 4;

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

double ic_concentration(double x, double y, scalar& dx, scalar& dy)
{
  if(INITIAL_CONCENTRATION_STATE == 0)
    return CONCENTRATION_EXT;
  else
    return 0.0;
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
  // Provide a possibility to change INITIAL_CONCENTRATION_STATE through an argument.
  if(argc > 1)
    INITIAL_CONCENTRATION_STATE = atoi(argv[1]);

  if(argc > 2)
    INIT_REF_NUM_FLOW = atoi(argv[2]);

  if(argc > 3)
    INIT_REF_NUM_CONCENTRATION = atoi(argv[3]);
  
  if(argc > 4)
    INIT_REF_NUM_CONCENTRATION_BDY = atoi(argv[4]);

  // Load the mesh.
  Mesh basemesh;
  H2DReader mloader;
  if(INITIAL_CONCENTRATION_STATE == 1)
    mloader.load("GAMM-channel-4-bnds.mesh", &basemesh);
  else
    mloader.load("channel-4-bnds.mesh", &basemesh);


  // Initialize the meshes.
  Mesh mesh_flow, mesh_concentration;
  mesh_flow.copy(&basemesh);
  mesh_concentration.copy(&basemesh);

  for(unsigned int i = 0; i < INIT_REF_NUM_CONCENTRATION; i++)
    mesh_concentration.refine_all_elements();

  if(INITIAL_CONCENTRATION_STATE != 2)
    mesh_concentration.refine_towards_boundary(3, INIT_REF_NUM_CONCENTRATION_BDY);
  
  //mesh_concentration.refine_towards_boundary(1, INIT_REF_NUM_CONCENTRATION_BDY, false);

  for(unsigned int i = 0; i < INIT_REF_NUM_FLOW; i++)
    mesh_flow.refine_all_elements();

  // Initialize boundary condition types and spaces with default shapesets.
  BCTypes bc_types_euler;
  bc_types_euler.add_bc_neumann(Hermes::vector<int>(BDY_SOLID_WALL_TOP, BDY_SOLID_WALL_BOTTOM, BDY_INLET, BDY_OUTLET));

  BCTypes bc_types_concentration;
  BCValues bc_values_concentration;

  switch(INITIAL_CONCENTRATION_STATE) {
  case 0:
    bc_types_concentration.add_bc_neumann(Hermes::vector<int>(BDY_INLET, BDY_OUTLET, BDY_SOLID_WALL_TOP));
    bc_types_concentration.add_bc_dirichlet(Hermes::vector<int>(BDY_SOLID_WALL_BOTTOM));
    bc_values_concentration.add_const(Hermes::vector<int>(BDY_SOLID_WALL_BOTTOM), CONCENTRATION_EXT);
    break;
  case 1:
    bc_types_concentration.add_bc_neumann(Hermes::vector<int>(BDY_INLET, BDY_OUTLET, BDY_SOLID_WALL_TOP));
    bc_types_concentration.add_bc_dirichlet(Hermes::vector<int>(BDY_SOLID_WALL_BOTTOM));
    bc_values_concentration.add_const(Hermes::vector<int>(BDY_SOLID_WALL_BOTTOM), CONCENTRATION_EXT);
    break;
  case 2:
    bc_types_concentration.add_bc_neumann(Hermes::vector<int>(BDY_SOLID_WALL_BOTTOM, BDY_OUTLET, BDY_SOLID_WALL_TOP));
    bc_types_concentration.add_bc_dirichlet(Hermes::vector<int>(BDY_INLET));
    bc_values_concentration.add_const(Hermes::vector<int>(BDY_INLET), CONCENTRATION_EXT);
    break;
  }

  L2Space space_rho(&mesh_flow, &bc_types_euler, P_INIT_FLOW);
  L2Space space_rho_v_x(&mesh_flow, &bc_types_euler, P_INIT_FLOW);
  L2Space space_rho_v_y(&mesh_flow, &bc_types_euler, P_INIT_FLOW);
  L2Space space_e(&mesh_flow, &bc_types_euler, P_INIT_FLOW);
  // Space for concentration.
  H1Space space_c(&mesh_concentration, &bc_types_concentration, &bc_values_concentration, P_INIT_CONCENTRATION);

  // Initialize solutions, set initial conditions.
  Solution sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e, sln_c, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e, prev_c;
  sln_rho.set_exact(&mesh_flow, ic_density);
  sln_rho_v_x.set_exact(&mesh_flow, ic_density_vel_x);
  sln_rho_v_y.set_exact(&mesh_flow, ic_density_vel_y);
  sln_e.set_exact(&mesh_flow, ic_energy);
  sln_c.set_exact(&mesh_concentration, ic_concentration);
  prev_rho.set_exact(&mesh_flow, ic_density);
  prev_rho_v_x.set_exact(&mesh_flow, ic_density_vel_x);
  prev_rho_v_y.set_exact(&mesh_flow, ic_density_vel_y);
  prev_e.set_exact(&mesh_flow, ic_energy);
  prev_c.set_exact(&mesh_concentration, ic_concentration);

  // Initialize weak formulation.
  bool is_matrix_free = true;
  WeakForm wf(5, is_matrix_free);

  // Volumetric linear forms.
  wf.add_vector_form(0, callback(linear_form_0_time));
  wf.add_vector_form(1, callback(linear_form_1_time));
  wf.add_vector_form(2, callback(linear_form_2_time));
  wf.add_vector_form(3, callback(linear_form_3_time));
  wf.add_vector_form(4, callback(linear_form_4_time));

  // Volumetric linear forms.
  // Linear forms coming from the linearization by taking the Eulerian fluxes' Jacobian matrices 
  // from the previous time step.
  // Unnecessary for FVM.
  if(P_INIT_FLOW.order_h > 0 || P_INIT_FLOW.order_v > 0) {
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
  wf.add_vector_form(4, callback(linear_form_time_concentration), HERMES_ANY, &prev_c);

  // Surface linear forms - inner edges coming from the DG formulation.
  wf.add_vector_form_surf(0, linear_form_interface_0, linear_form_order, H2D_DG_INNER_EDGE);
  wf.add_vector_form_surf(1, linear_form_interface_1, linear_form_order, H2D_DG_INNER_EDGE);
  wf.add_vector_form_surf(2, linear_form_interface_2, linear_form_order, H2D_DG_INNER_EDGE);
  wf.add_vector_form_surf(3, linear_form_interface_3, linear_form_order, H2D_DG_INNER_EDGE);

  // Surface linear forms - inlet / outlet edges.
  wf.add_vector_form_surf(0, bdy_flux_inlet_outlet_comp_0, linear_form_order, BDY_INLET);
  wf.add_vector_form_surf(1, bdy_flux_inlet_outlet_comp_1, linear_form_order, BDY_INLET);
  wf.add_vector_form_surf(2, bdy_flux_inlet_outlet_comp_2, linear_form_order, BDY_INLET);
  wf.add_vector_form_surf(3, bdy_flux_inlet_outlet_comp_3, linear_form_order, BDY_INLET);

  wf.add_vector_form_surf(0, bdy_flux_inlet_outlet_comp_0, linear_form_order, BDY_OUTLET);
  wf.add_vector_form_surf(1, bdy_flux_inlet_outlet_comp_1, linear_form_order, BDY_OUTLET);
  wf.add_vector_form_surf(2, bdy_flux_inlet_outlet_comp_2, linear_form_order, BDY_OUTLET);
  wf.add_vector_form_surf(3, bdy_flux_inlet_outlet_comp_3, linear_form_order, BDY_OUTLET);
  
  // Surface linear forms - Solid wall edges.
  wf.add_vector_form_surf(0, bdy_flux_solid_wall_comp_0, linear_form_order, BDY_SOLID_WALL_BOTTOM);
  wf.add_vector_form_surf(1, bdy_flux_solid_wall_comp_1, linear_form_order, BDY_SOLID_WALL_BOTTOM);
  wf.add_vector_form_surf(2, bdy_flux_solid_wall_comp_2, linear_form_order, BDY_SOLID_WALL_BOTTOM);
  wf.add_vector_form_surf(3, bdy_flux_solid_wall_comp_3, linear_form_order, BDY_SOLID_WALL_BOTTOM);

  wf.add_vector_form_surf(0, bdy_flux_solid_wall_comp_0, linear_form_order, BDY_SOLID_WALL_TOP);
  wf.add_vector_form_surf(1, bdy_flux_solid_wall_comp_1, linear_form_order, BDY_SOLID_WALL_TOP);
  wf.add_vector_form_surf(2, bdy_flux_solid_wall_comp_2, linear_form_order, BDY_SOLID_WALL_TOP);
  wf.add_vector_form_surf(3, bdy_flux_solid_wall_comp_3, linear_form_order, BDY_SOLID_WALL_TOP);

  // Forms for concentration.
  wf.add_vector_form(4, callback(volume_linear_form_concentration_grad_grad), HERMES_ANY);
  
  wf.add_vector_form(4, callback(volume_linear_form_concentration_convective), HERMES_ANY);

  wf.add_vector_form_surf(4, callback(surface_linear_form_concentration_inlet_outlet), BDY_INLET);

  wf.add_vector_form_surf(4, callback(surface_linear_form_concentration_inlet_outlet), BDY_OUTLET);

  wf.add_vector_form_surf(4, callback(surface_linear_form_concentration_solid_wall), BDY_SOLID_WALL_TOP);

  wf.add_vector_form_surf(4, callback(inner_linear_form_concentration), H2D_DG_INNER_EDGE);

  if(PRECONDITIONING) {
    // Preconditioning forms.
    wf.add_matrix_form(0, 0, callback(bilinear_form_precon));
    wf.add_matrix_form(1, 1, callback(bilinear_form_precon));
    wf.add_matrix_form(2, 2, callback(bilinear_form_precon));
    wf.add_matrix_form(3, 3, callback(bilinear_form_precon));
    wf.add_matrix_form(4, 4, callback(bilinear_form_precon));
  }

  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c), is_linear);
  // Project the initial solution on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution on the FE mesh.");
  scalar* coeff_vec = new scalar[Space::get_num_dofs(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c))];
  OGProjection::project_global(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c), 
    Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c), coeff_vec);
  
  // Filters for visualization of pressure and the two components of velocity.
  /*
  SimpleFilter pressure(calc_pressure_func, Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter u(calc_u_func, Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter w(calc_w_func, Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter Mach_number(calc_Mach_func, Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter entropy_estimate(calc_entropy_estimate_func, Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));
  VectorView vview("Velocity", new WinGeom(700, 400, 600, 300));
  */

  ScalarView s1("w0", new WinGeom(0, 0, 600, 300));
  ScalarView s2("w1", new WinGeom(700, 0, 600, 300));
  ScalarView s3("w2", new WinGeom(0, 400, 600, 300));
  ScalarView s4("w3", new WinGeom(700, 400, 600, 300));
  ScalarView s5("Concentration", new WinGeom(350, 200, 600, 300));
  
  // Iteration number.
  int iteration = 0;
  
  // Output of the approximate time derivative.
  std::ofstream time_der_out("time_der");

  // Initialize NOX solver.
  NoxSolver solver(&dp);
  solver.set_ls_tolerance(1E-2);
  solver.disable_abs_resid();
  solver.set_conv_rel_resid(1.00);

  // Select preconditioner.
  if(PRECONDITIONING) {
    RCP<Precond> pc = rcp(new MlPrecond("sa"));
    solver.set_precond(pc);
  }

  for(t = 0.0; t < 3.0; t += TAU)
  {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    OGProjection::project_global(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c), 
    Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c), coeff_vec);

    info("Assembling by DiscreteProblem, solving by NOX.");
    solver.set_init_sln(coeff_vec);
    if (solver.solve())
      Solution::vector_to_solutions(solver.get_solution(), Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c), 
      Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c));
    else
      error("NOX failed.");

    /*
    // Visualization.
    pressure.reinit();
    u.reinit();
    w.reinit();
    Mach_number.reinit();
    entropy_estimate.reinit();
    pressure_view.show(&pressure);
    entropy_production_view.show(&entropy_estimate);
    Mach_number_view.show(&Mach_number);
    vview.show(&u, &w);
    */

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      solver.get_num_lin_iters(), solver.get_achieved_tol());

    // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0) {
      // Hermes visualization.
      if(HERMES_VISUALIZATION) {
        s1.show(&prev_rho);
        s2.show(&prev_rho_v_x);
        s3.show(&prev_rho_v_y);
        s4.show(&prev_e);
        s5.show(&prev_c);
      }
      // Output solution in VTK format.
      if(VTK_OUTPUT) {
        Linearizer lin;
        char filename[40];
        sprintf(filename, "w0-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_rho, filename, "w0", false);
        sprintf(filename, "w1-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_rho_v_x, filename, "w1", false);
        sprintf(filename, "w2-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_rho_v_y, filename, "w2", false);
        sprintf(filename, "w3-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_e, filename, "w3", false);
        sprintf(filename, "concentration-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_c, filename, "concentration", false);
      }
    }
  }
  
  s1.close();
  s2.close();
  s3.close();
  s4.close();
  s5.close();
  time_der_out.close();
  return 0;
}