#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Teuchos;
using namespace RefinementSelectors;

//  This example solves the compressible Euler equations using the JFNK
//  method implemented in the NOX package of the Trilinos library.
//
//  Equations: Compressible Euler equations, perfect gas state equation.
//
//  Domain: GAMM channel, see mesh file GAMM-channel.mesh
//
//  BC: Normal velocity component is zero on solid walls.
//      Subsonic state prescribed on inlet and outlet.
//
//  IC: Constant subsonic state identical to inlet. 
//
//  The following parameters can be changed:

// Calculation of approximation of time derivative (and its output).
// Setting this option to false saves the computation time.
const bool CALC_TIME_DER = true;

const Ord2 P_INIT = Ord2(0,0);                    // Initial polynomial degree.                      
const int INIT_REF_NUM = 4;                       // Number of initial uniform mesh refinements.                       
double TAU = 1E-2;                                // Time step.

// Equation parameters.
double P_EXT = 2.5;                               // Exterior pressure (dimensionless).
double RHO_EXT = 1.0;                             // Inlet density (dimensionless).   
double V1_EXT = 1.25;                             // Inlet x-velocity (dimensionless).
double V2_EXT = 0.0;                              // Inlet y-velocity (dimensionless).
double KAPPA = 1.4;                               // Kappa.

double t = 0;

// Boundary markers.
const int BDY_SOLID_WALL = 1;
const int BDY_INLET_OUTLET = 2;

// Numerical flux.
// For numerical fluxes, please see hermes2d/src/numerical_flux.h
NumericalFlux num_flux(KAPPA);

// Inlet/outlet boundary conditions.
double bc_density(double y)
{
  return RHO_EXT;
}

// Density * velocity in the x coordinate boundary condition.
double bc_density_vel_x(double y)
{
  return RHO_EXT * V1_EXT;
}

// Density * velocity in the y coordinate boundary condition.
double bc_density_vel_y(double y)
{
  return V2_EXT;
}

// Calculation of the pressure on the boundary.
double bc_pressure(double y)
{
  return P_EXT;
}

// Energy boundary condition.
double bc_energy(double y)
{
  double rho = bc_density(y);
  double rho_v_x = bc_density_vel_x(y);
  double rho_v_y = bc_density_vel_y(y);
  double pressure = bc_pressure(y);
  return pressure/(num_flux.kappa - 1.) + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / 2*rho;
}

// Calculates energy from other quantities.
// FIXME: this should be in the src/ directory, not here.
double calc_energy(double rho, double rho_v_x, double rho_v_y, double pressure)
{
  return pressure/(num_flux.kappa - 1.) + (rho_v_x*rho_v_x+rho_v_y*rho_v_y) / 2*rho;
}

// Calculates pressure from other quantities.
// FIXME: this should be in the src/ directory, not here.
double calc_pressure(double rho, double rho_v_x, double rho_v_y, double energy)
{
  return (num_flux.kappa - 1.) * (energy - (rho_v_x*rho_v_x + rho_v_y*rho_v_y) / (2*rho));
}

// Calculates speed of sound.
// FIXME: this should be in the src/ directory, not here.
double calc_sound_speed(double rho, double rho_v_x, double rho_v_y, double energy)
{
  return std::sqrt(num_flux.kappa * calc_pressure(rho, rho_v_x, rho_v_y, energy) / rho);
}

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

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("GAMM-channel.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_neumann(Hermes::vector<int>(BDY_SOLID_WALL, BDY_INLET_OUTLET));

  // Create L2 spaces with default shapesets.
  L2Space space_rho(&mesh, &bc_types, P_INIT);
  L2Space space_rho_v_x(&mesh, &bc_types, P_INIT);
  L2Space space_rho_v_y(&mesh, &bc_types, P_INIT);
  L2Space space_e(&mesh, &bc_types, P_INIT);

  // Initialize solutions, set initial conditions.
  Solution prev_time_rho, prev_time_rho_v_x, prev_time_rho_v_y, prev_time_e;
  prev_time_rho.set_exact(&mesh, ic_density);
  prev_time_rho_v_x.set_exact(&mesh, ic_density_vel_x);
  prev_time_rho_v_y.set_exact(&mesh, ic_density_vel_y);
  prev_time_e.set_exact(&mesh, ic_energy);

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

  // Matrix form for preconditioning.
  wf.add_matrix_form(0, 0, callback(bilinear_form_precon));
  wf.add_matrix_form(1, 1, callback(bilinear_form_precon));
  wf.add_matrix_form(2, 2, callback(bilinear_form_precon));
  wf.add_matrix_form(3, 3, callback(bilinear_form_precon));

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
  wf.add_vector_form(0, linear_form, linear_form_order, HERMES_ANY, &prev_time_rho);
  wf.add_vector_form(1, linear_form, linear_form_order, HERMES_ANY, &prev_time_rho_v_x);
  wf.add_vector_form(2, linear_form, linear_form_order, HERMES_ANY, &prev_time_rho_v_y);
  wf.add_vector_form(3, linear_form, linear_form_order, HERMES_ANY, &prev_time_e);

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


  // Initialize the FE problem.
  bool is_linear = false;
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), is_linear);
  
  // If the FE problem is in fact a FV problem.
  if(P_INIT.order_h == 0 && P_INIT.order_v == 0)
    dp.set_fvm();

  // Project the initial solution on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution on the FE mesh.");
  scalar* coeff_vec = new scalar[Space::get_num_dofs(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e))];
  OGProjection::project_global(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), 
    Hermes::vector<MeshFunction*>(&prev_time_rho, &prev_time_rho_v_x, &prev_time_rho_v_y, &prev_time_e), coeff_vec);

  // Filters for visualization of pressure and the two components of velocity.
  SimpleFilter pressure(calc_pressure_func, Hermes::vector<MeshFunction*>(&prev_time_rho, &prev_time_rho_v_x, &prev_time_rho_v_y, &prev_time_e));
  SimpleFilter u(calc_u_func, Hermes::vector<MeshFunction*>(&prev_time_rho, &prev_time_rho_v_x, &prev_time_rho_v_y, &prev_time_e));
  SimpleFilter w(calc_w_func, Hermes::vector<MeshFunction*>(&prev_time_rho, &prev_time_rho_v_x, &prev_time_rho_v_y, &prev_time_e));
  SimpleFilter Mach_number(calc_Mach_func, Hermes::vector<MeshFunction*>(&prev_time_rho, &prev_time_rho_v_x, &prev_time_rho_v_y, &prev_time_e));
  SimpleFilter entropy_estimate(calc_entropy_estimate_func, Hermes::vector<MeshFunction*>(&prev_time_rho, &prev_time_rho_v_x, &prev_time_rho_v_y, &prev_time_e));

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));
  VectorView vview("Velocity", new WinGeom(700, 400, 600, 300));

  // Iteration number.
  int iteration = 0;
  
  // Output of the approximate time derivative.
  std::ofstream time_der_out("time_der");

  // Initialize NOX solver.
  NoxSolver solver(&dp);

  // Select preconditioner.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  solver.set_precond(pc);

  for(t = 0.0; t < 10; t += TAU)
  {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    OGProjection::project_global(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), 
    Hermes::vector<MeshFunction*>(&prev_time_rho, &prev_time_rho_v_x, &prev_time_rho_v_y, &prev_time_e), coeff_vec);

    info("Assembling by DiscreteProblem, solving by NOX.");
    solver.set_init_sln(coeff_vec);
    if (solver.solve())
      Solution::vector_to_solutions(solver.get_solution(), Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), 
      Hermes::vector<Solution *>(&prev_time_rho, &prev_time_rho_v_x, &prev_time_rho_v_y, &prev_time_e));
    else
      error("NOX failed.");

    // Approximate the time derivative of the solution.
    if(CALC_TIME_DER) {
      Adapt *adapt_for_time_der_calc = new Adapt(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e));
      bool solutions_for_adapt = false;
      double difference = iteration == 1 ? 0 : 
        adapt_for_time_der_calc->calc_err_est(Hermes::vector<Solution *>(&sln_temp_rho, &sln_temp_rho_v_x, &sln_temp_rho_v_y, &sln_temp_e), 
                                              Hermes::vector<Solution *>(&prev_time_rho, &prev_time_rho_v_x, &prev_time_rho_v_y, &prev_time_e), 
                                              (Hermes::vector<double>*) NULL, solutions_for_adapt, 
                                              HERMES_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_ABS) / TAU;
      delete adapt_for_time_der_calc;
      sln_temp_rho.copy(&prev_time_rho);
      sln_temp_rho_v_x.copy(&prev_time_rho_v_x);
      sln_temp_rho_v_y.copy(&prev_time_rho_v_y);
      sln_temp_e.copy(&prev_time_e);

      // Info about the approximate time derivative.
      if(iteration > 1)
      {
        info("Approximate the norm time derivative : %g.", difference);
        time_der_out << iteration << '\t' << difference << std::endl;
      }
    }

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

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      solver.get_num_lin_iters(), solver.get_achieved_tol());
  }
  
  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();
  vview.close();

  time_der_out.close();
  return 0;
}
