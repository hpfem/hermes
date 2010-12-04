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

const int P_INIT = 0;                             // Initial polynomial degree.                      
const int INIT_REF_NUM = 4;                       // Number of initial uniform mesh refinements.                       
double CFL = 0.8;                                 // CFL value.
double TAU = 1E-4;                                // Time step.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC, SOLVER_MUMPS, 
                                                  // SOLVER_PARDISO, SOLVER_SUPERLU, SOLVER_AMESOS, SOLVER_AZTECOO

// Equation parameters.
double P_EXT = 2.5;         // Exterior pressure (dimensionless).
double RHO_EXT = 1.0;       // Inlet density (dimensionless).   
double V1_EXT = 1.25;       // Inlet x-velocity (dimensionless).
double V2_EXT = 0.0;        // Inlet y-velocity (dimensionless).
double KAPPA = 1.4;         // Kappa.

double t = 0;

// Mesh boundary markers.
#define BDY_SOLID_WALL 1
#define BDY_INLET_OUTLET 2

// Numerical flux.
// For numerical fluxes, please see hermes2d/src/numerical_flux.h
NumericalFlux num_flux(KAPPA);

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_NATURAL;
}

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

int main(int argc, char* argv[])
{
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

  // Initialize spaces with default shapesets.
  L2Space space_rho(&mesh, bc_types, NULL, P_INIT);
  L2Space space_rho_v_x(&mesh, bc_types, NULL, P_INIT);
  L2Space space_rho_v_y(&mesh, bc_types, NULL, P_INIT);
  L2Space space_e(&mesh, bc_types, NULL, P_INIT);

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
  if(P_INIT > 0) {
    wf.add_vector_form(0, callback(linear_form_0_1), HERMES_ANY, Hermes::Tuple<MeshFunction*>(&prev_rho_v_x));
    wf.add_vector_form(1, callback(linear_form_1_0_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_1_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_2_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_3_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(2, callback(linear_form_2_0_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_1_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_2_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_3_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_0_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_1_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_2_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_3_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    // Second flux.
    
    wf.add_vector_form(0, callback(linear_form_0_2), HERMES_ANY, Hermes::Tuple<MeshFunction*>(&prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_0_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_1_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_2_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_3_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(2, callback(linear_form_2_0_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_1_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_2_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_3_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_0_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_1_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_2_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_3_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  }

  // Volumetric linear forms coming from the time discretization.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form(0, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(1, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(2, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form(3, linear_form_vector, linear_form_order, HERMES_ANY, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form(0, linear_form, linear_form_order, HERMES_ANY, &prev_rho);
  wf.add_vector_form(1, linear_form, linear_form_order, HERMES_ANY, &prev_rho_v_x);
  wf.add_vector_form(2, linear_form, linear_form_order, HERMES_ANY, &prev_rho_v_y);
  wf.add_vector_form(3, linear_form, linear_form_order, HERMES_ANY, &prev_e);
#endif

  // Surface linear forms - inner edges coming from the DG formulation.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form_surf(0, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form_surf(0, linear_form_interface_0, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, linear_form_interface_1, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, linear_form_interface_2, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, linear_form_interface_3, linear_form_order, H2D_DG_INNER_EDGE, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#endif

  // Surface linear forms - inlet / outlet edges.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form_surf(0, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_inlet_outlet_comp_vector, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form_surf(0, bdy_flux_inlet_outlet_comp_0, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_inlet_outlet_comp_1, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_inlet_outlet_comp_2, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_inlet_outlet_comp_3, linear_form_order, BDY_INLET_OUTLET, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#endif
  
  // Surface linear forms - Solid wall edges.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  wf.add_vector_form_surf(0, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_solid_wall_comp_vector, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#else
  wf.add_vector_form_surf(0, bdy_flux_solid_wall_comp_0, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(1, bdy_flux_solid_wall_comp_1, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(2, bdy_flux_solid_wall_comp_2, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  wf.add_vector_form_surf(3, bdy_flux_solid_wall_comp_3, linear_form_order, BDY_SOLID_WALL, 
                          Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
#endif

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, Hermes::Tuple<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), is_linear);
  // If the FE problem is in fact a FV problem.
  if(P_INIT == 0)
    dp.set_fvm();
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
  dp.use_vector_valued_forms();
#endif


  // Filters for visualization of pressure and the two components of velocity.
  SimpleFilter pressure(calc_pressure_func, Hermes::Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter u(calc_u_func, Hermes::Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter w(calc_w_func, Hermes::Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));

  //VectorView vview("Velocity", new WinGeom(0, 0, 600, 300));
  //ScalarView sview("Pressure", new WinGeom(700, 0, 600, 300));

  ScalarView s1("w1", new WinGeom(0, 0, 620, 300));
  s1.fix_scale_width(80);
  ScalarView s2("w2", new WinGeom(625, 0, 600, 300));
  s2.fix_scale_width(50);
  ScalarView s3("w3", new WinGeom(0, 350, 620, 300));
  s3.fix_scale_width(80);
  ScalarView s4("w4", new WinGeom(625, 350, 600, 300));
  s4.fix_scale_width(50);

  // Iteration number.
  int iteration = 0;
  
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Output of the approximate time derivative.
  std::ofstream time_der_out("time_der");

  for(t = 0.0; t < 10; t += TAU)
  {
    info("---- Time step %d, time %3.5f.", iteration, t);

    iteration++;

    bool rhs_only = (iteration == 1 ? false : true);
    // Assemble stiffness matrix and rhs or just rhs.
    if (rhs_only == false) info("Assembling the stiffness matrix and right-hand side vector.");
    else info("Assembling the right-hand side vector (only).");
    dp.assemble(matrix, rhs, rhs_only);

        
    // Solve the matrix problem.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution::vector_to_solutions(solver->get_solution(), Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), Hermes::Tuple<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
    else
    error ("Matrix solver failed.\n");

    // Approximate the time derivative of the solution.
    if(CALC_TIME_DER) {
      Adapt *adapt_for_time_der_calc = new Adapt(Hermes::Tuple<Space *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::Tuple<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      bool solutions_for_adapt = false;
      double difference = adapt_for_time_der_calc->calc_err_est(Hermes::Tuple<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
        Hermes::Tuple<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), solutions_for_adapt, HERMES_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_ABS) / TAU;
      delete adapt_for_time_der_calc;

      // Info about the approximate time derivative.
      if(iteration > 1)
      {
        info("Approximate the norm time derivative : %g.", difference);
        time_der_out << iteration << '\t' << difference << std::endl;
      }
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

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&sln_rho);
    prev_rho_v_x.copy(&sln_rho_v_x);
    prev_rho_v_y.copy(&sln_rho_v_y);
    prev_e.copy(&sln_e);

    // Visualization.
    /*
    pressure.reinit();
    u.reinit();
    w.reinit();
    sview.show(&pressure);
    vview.show(&u, &w);
    */

    s1.show(&sln_rho);
    s2.show(&sln_rho_v_x);
    s3.show(&sln_rho_v_y);
    s4.show(&sln_e);
    
    // If used, we need to clean the vector valued form caches.
#ifdef HERMES_USE_VECTOR_VALUED_FORMS
    DiscreteProblem::empty_form_caches();
#endif
  }
  
  time_der_out.close();
  return 0;
}
