#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

//  This example solves the compressible Euler equations using a basic
//  piecewise-constant finite volume method.
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

// Experimental caching of vector valued (vector) forms.
#define HERMES_USE_VECTOR_VALUED_FORMS

// Calculation of approximation of time derivative (and its output).
// Setting this option to false saves the computation time.
const bool CALC_TIME_DER = true;

const int INIT_REF_NUM = 4;                       // Number of initial uniform mesh refinements.                       
double CFL = 0.8;                                 // CFL value.
double TAU = 1E-4;                                // Time step.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Equation parameters.
double P_EXT = 2.5;         // Exterior pressure (dimensionless).
double RHO_EXT = 1.0;       // Inlet density (dimensionless).   
double V1_EXT = 1.25;       // Inlet x-velocity (dimensionless).
double V2_EXT = 0.0;        // Inlet y-velocity (dimensionless).
double KAPPA = 1.4;         // Kappa.

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
  /* So far the DG assembling is very slow for higher order polynomials,
      so only constant functions are used here.
  if(argc < 3)
    error("Too few arguments in example-euler-gamm-explicit");

  Ord2 P_INIT= Ord2(atoi(argv[1]), atoi(argv[2]));
  */

  Ord2 P_INIT= Ord2(0, 0);

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("GAMM-channel.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, 1);
  mesh.refine_element(1053);
  mesh.refine_element(1054);
  mesh.refine_element(1087);
  mesh.refine_element(1088);
  mesh.refine_element(1117);
  mesh.refine_element(1118);
  mesh.refine_element(1151);
  mesh.refine_element(1152);

  // Enter boundary markers.  
  BCTypes bc_types;
  bc_types.add_bc_neumann(Hermes::vector<int>(BDY_SOLID_WALL, BDY_INLET_OUTLET));

  // Create L2 spaces with default shapesets.
  L2Space space_rho(&mesh, &bc_types, P_INIT);
  L2Space space_rho_v_x(&mesh, &bc_types, P_INIT);
  L2Space space_rho_v_y(&mesh, &bc_types, P_INIT);
  L2Space space_e(&mesh, &bc_types, P_INIT);

  // Initialize solutions, set initial conditions.
  Solution sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e;
  sln_rho.set_exact(&mesh, ic_density);
  sln_rho_v_x.set_exact(&mesh, ic_density_vel_x);
  sln_rho_v_y.set_exact(&mesh, ic_density_vel_y);
  sln_e.set_exact(&mesh, ic_energy);
  prev_rho.set_exact(&mesh, ic_density);
  prev_rho_v_x.set_exact(&mesh, ic_density_vel_x);
  prev_rho_v_y.set_exact(&mesh, ic_density_vel_y);
  prev_e.set_exact(&mesh, ic_energy);

  // Initialize weak formulation.
  WeakForm wf(4);

  // Bilinear forms coming from time discretization by explicit Euler's method.
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0_time));
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1_time));
  wf.add_matrix_form(2, 2, callback(bilinear_form_2_2_time));
  wf.add_matrix_form(3, 3, callback(bilinear_form_3_3_time));

  // Volumetric linear forms.
  // Linear forms coming from the linearization by taking the Eulerian fluxes' Jacobian matrices 
  // from the previous time step.
  // First flux.
  // Unnecessary for FVM.
  if(P_INIT.order_h > 0 || P_INIT.order_v > 0) {
    wf.add_vector_form(0, callback(linear_form_0_1), HERMES_ANY, Hermes::vector<MeshFunction*>(&prev_rho_v_x));
    
    wf.add_vector_form(1, callback(linear_form_1_0_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_1_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_2_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_3_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(2, callback(linear_form_2_0_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_1_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_2_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_3_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_0_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_1_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_2_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_3_first_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    // Second flux.
    
    wf.add_vector_form(0, callback(linear_form_0_2), HERMES_ANY, Hermes::vector<MeshFunction*>(&prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_0_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_1_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_2_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_3_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(2, callback(linear_form_2_0_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_1_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_2_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_3_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_0_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_1_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_2_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_3_second_flux), HERMES_ANY, 
                       Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  }

  // Volumetric linear forms coming from the time discretization.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form(0, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(1, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(2, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(3, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form(0, linear_form, linear_form_order, HERMES_ANY, &prev_rho);
  wf.add_vector_form(1, linear_form, linear_form_order, HERMES_ANY, &prev_rho_v_x);
  wf.add_vector_form(2, linear_form, linear_form_order, HERMES_ANY, &prev_rho_v_y);
  wf.add_vector_form(3, linear_form, linear_form_order, HERMES_ANY, &prev_e);
#endif

  // Surface linear forms - inner edges coming from the DG formulation.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form_surf(0, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form_surf(0, linear_form_interface_0, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, linear_form_interface_1, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, linear_form_interface_2, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, linear_form_interface_3, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#endif

  // Surface linear forms - inlet / outlet edges.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form_surf(0, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form_surf(0, bdy_flux_inlet_outlet_comp_0, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_inlet_outlet_comp_1, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_inlet_outlet_comp_2, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_inlet_outlet_comp_3, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#endif
  
  // Surface linear forms - Solid wall edges.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form_surf(0, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form_surf(0, bdy_flux_solid_wall_comp_0, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_solid_wall_comp_1, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_solid_wall_comp_2, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_solid_wall_comp_3, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#endif

  // Initialize the FE problem.
  bool is_linear = true;
  
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), is_linear);
  
  // If the FE problem is in fact a FV problem.
  if(P_INIT.order_h == 0 && P_INIT.order_v == 0)
    dp.set_fvm();
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  dp.use_vector_valued_forms();
#endif

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

  // Iteration number.
  int iteration = 0;
  
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // For testing purposes.
  double l2_norms[5][4];
  for(unsigned int i = 0; i < 5; i++)
    for(unsigned int j = 0; j < 4; j++)
      l2_norms[i][j] = 0.0;
  double point_values[5][3];
  for(unsigned int i = 0; i < 5; i++)
    for(unsigned int j = 0; j < 3; j++)
      point_values[i][j] = 0.0;
  // Calculate the special point where we will evaluate the solution.
  double x = 0.75;
  double y = sqrt((double)(1.-x*x)) + 0.001;

  for(unsigned int time_step = 0; time_step < 5; time_step++)
  {
    iteration++;

    bool rhs_only = (iteration == 1 ? false : true);
    // Assemble stiffness matrix and rhs or just rhs.
    if (rhs_only == false) info("Assembling the stiffness matrix and right-hand side vector.");
    else info("Assembling the right-hand side vector (only).");
    dp.assemble(matrix, rhs, rhs_only);

        
    // Solve the matrix problem.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution::vector_to_solutions(solver->get_solution(), Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), Hermes::vector<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
    else
    error ("Matrix solver failed.\n");

    // Approximate the time derivative of the solution.
    if(CALC_TIME_DER) {
      Adapt *adapt_for_time_der_calc = new Adapt(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
      bool solutions_for_adapt = false;
      double difference = iteration == 1 ? 0 : 
        adapt_for_time_der_calc->calc_err_est(Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
					      Hermes::vector<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), 
                                              (Hermes::vector<double>*) NULL, solutions_for_adapt, 
                                              HERMES_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_ABS) / TAU;
      delete adapt_for_time_der_calc;
    }

    // Determine the time step according to the CFL condition.
    // Only mean values on an element of each solution component are taken into account.
    double *solution_vector = solver->get_solution();
    double min_condition = 0;
    Element *e;
    for (int _id = 0, _max = mesh.get_max_element_id(); _id < _max; _id++) \
          if (((e) = mesh.get_element_fast(_id))->used) \
            if ((e)->active)
    {
      AsmList al;
      space_rho.get_element_assembly_list(e, &al);
      double rho = solution_vector[al.dof[0]];
      space_rho_v_x.get_element_assembly_list(e, &al);
      double v1 = solution_vector[al.dof[0]] / rho;
      space_rho_v_y.get_element_assembly_list(e, &al);
      double v2 = solution_vector[al.dof[0]] / rho;
      space_e.get_element_assembly_list(e, &al);
      double energy = solution_vector[al.dof[0]];
      
      double condition = e->get_area() / (std::sqrt(v1*v1 + v2*v2) + calc_sound_speed(rho, rho*v1, rho*v2, energy));
      
      if(condition < min_condition || min_condition == 0.)
        min_condition = condition;
    }
    if(TAU > min_condition)
      TAU = min_condition;
    if(TAU < min_condition * 0.9)
      TAU = min_condition;

    // Storing the testing values.
    for(unsigned int j = 0; j < 4; j++)
      for(unsigned int k = j*space_rho.get_num_dofs(); k < (j+1)*space_rho.get_num_dofs(); k++)
        l2_norms[time_step][j] += solver->get_solution()[k];    

    point_values[time_step][0] = sln_rho_v_x.get_pt_value(0.5, 0.001);
    point_values[time_step][1] = sln_rho_v_x.get_pt_value(x, y);
    point_values[time_step][2] = sln_rho_v_x.get_pt_value(1.5, 0.001);

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&sln_rho);
    prev_rho_v_x.copy(&sln_rho_v_x);
    prev_rho_v_y.copy(&sln_rho_v_y);
    prev_e.copy(&sln_e);

// If used, we need to clean the vector valued form caches.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
    DiscreteProblem::empty_form_caches();
#endif
  }
  bool okay = true;
  switch(P_INIT.order_h* 10 + P_INIT.order_v) {
  case 0:
    if(std::abs(l2_norms[0][0] - 888.0) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[0][1] - 1110) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[0][2]) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[0][3] - 6243.75) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[1][0] - 887.99997637865545) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[1][1] - 1109.9997956458228) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[1][2] - 3.1927018090871903e-008) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[1][3] - 6243.7496921971369) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[2][0] - 887.99993429457072) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[2][1] - 1109.9994322038613) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[2][2] + 5.3556469633245445e-008) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[2][3] - 6243.7491437826511) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[3][0] - 887.99987376550200) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[3][1] - 1109.9989102977672) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[3][2] + 3.6958140470412712e-007) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[3][3] - 6243.7483549661320) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[4][0] - 887.99979481320088) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[4][1] - 1109.9982305630808) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[4][2] + 1.0303296924184822e-006) > 1E-8)
      okay = false;
    if(std::abs(l2_norms[4][3] - 6243.7473260085462) > 1E-8)
      okay = false;

    // points
    if(std::abs(point_values[0][0] - 1.25) > 1E-8)
      okay = false;
    if(std::abs(point_values[0][1] - 1.25) > 1E-8)
      okay = false;
    if(std::abs(point_values[0][2] - 1.2459744738974898) > 1E-8)
      okay = false;
    if(std::abs(point_values[1][0] - 1.2499951035194972) > 1E-8)
      okay = false;
    if(std::abs(point_values[1][1] - 1.25) > 1E-8)
      okay = false;
    if(std::abs(point_values[1][2] - 1.2428402692519325) > 1E-8)
      okay = false;
    if(std::abs(point_values[2][0] - 1.2499864002215795) > 1E-8)
      okay = false;
    if(std::abs(point_values[2][1] - 1.25) > 1E-8)
      okay = false;
    if(std::abs(point_values[2][2] - 1.2397180160697001) > 1E-8)
      okay = false;
    if(std::abs(point_values[3][0] - 1.2499739085927257) > 1E-8)
      okay = false;
    if(std::abs(point_values[3][1] - 1.25) > 1E-8)
      okay = false;
    if(std::abs(point_values[3][2] - 1.2366079101139087) > 1E-8)
      okay = false;
    if(std::abs(point_values[4][0] - 1.2499576472516911) > 1E-8)
      okay = false;
    if(std::abs(point_values[4][1] - 1.25) > 1E-8)
      okay = false;
    if(std::abs(point_values[4][2] - 1.2335101392959738) > 1E-8)
      okay = false;
    break;
  }

  if (okay) {      // ndofs was 908 at the time this test was created
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  } 
}
