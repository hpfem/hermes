#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This example solves the compressible Euler equations coupled with an advection-diffution equation
// using a basic piecewise-constant finite volume method for the flow and continuous FEM for the concentration
// being advected by the flow.
//
// Equations: Compressible Euler equations, perfect gas state equation, advection-diffusion equation.
//
// Domains: Various
//
// BC: Various.
//
// IC: Various.
//
// The following parameters can be changed:
// Visualization.
const bool HERMES_VISUALIZATION = true;               // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;                  // Set to "true" to enable VTK output.
const unsigned int EVERY_NTH_STEP = 1;                // Set visual output for every nth step.

// Shock capturing.
bool SHOCK_CAPTURING = false;

// Quantitative parameter of the discontinuity detector.
double DISCONTINUITY_DETECTOR_PARAM = 1.0;

// Stability for the concentration part.
double ADVECTION_STABILITY_CONSTANT = 1.0;
const double DIFFUSION_STABILITY_CONSTANT = 1.0;

const int P_INIT_FLOW = 0;                             // Polynomial degree for the Euler equations (for the flow).
const int P_INIT_CONCENTRATION = 1;                    // Polynomial degree for the concentration.
double CFL_NUMBER = 1.0;                               // CFL value.
double time_step = 1E-5, util_time_step;               // Initial and utility time step.
const MatrixSolverType matrix_solver = SOLVER_UMFPACK; // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                       // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

unsigned int INIT_REF_NUM_FLOW = 3;                    // Number of initial uniform mesh refinements of the mesh for the flow.
unsigned int INIT_REF_NUM_CONCENTRATION = 3;           // Number of initial uniform mesh refinements of the mesh for the concentration.
unsigned int INIT_REF_NUM_CONCENTRATION_BDY = 1;       // Number of initial mesh refinements of the mesh for the concentration towards the 
                                                       // part of the boundary where the concentration is prescribed.
// Equation parameters.
const double P_EXT = 2.5;                              // Exterior pressure (dimensionless).
const double RHO_EXT = 1.0;                            // Inlet density (dimensionless).   
const double V1_EXT = 1.25;                             // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;                            // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;                              // Kappa.
const double CONCENTRATION_EXT = 0.1;                  // Concentration on the boundary.
const double CONCENTRATION_EXT_STARTUP_TIME = 0.0;     // Start time of the concentration on the boundary.

const double EPSILON = 0.01;                           // Diffusivity.

// Boundary markers.
const std::string BDY_INLET = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_BOTTOM = "3";
const std::string BDY_SOLID_WALL_TOP = "4";
const std::string BDY_DIRICHLET_CONCENTRATION = "3";
Hermes::vector<std::string> BDY_NATURAL_CONCENTRATION = Hermes::vector<std::string>("2", "1");

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh basemesh;
  H2DReader mloader;
  mloader.load("GAMM-channel.mesh", &basemesh);

  // Initialize the meshes.
  Mesh mesh_flow, mesh_concentration;
  mesh_flow.copy(&basemesh);
  mesh_concentration.copy(&basemesh);

  for(unsigned int i = 0; i < INIT_REF_NUM_CONCENTRATION; i++)
    mesh_concentration.refine_all_elements();

  mesh_concentration.refine_towards_boundary(BDY_DIRICHLET_CONCENTRATION, INIT_REF_NUM_CONCENTRATION_BDY);
  mesh_flow.refine_towards_boundary(BDY_DIRICHLET_CONCENTRATION, INIT_REF_NUM_CONCENTRATION_BDY);

  for(unsigned int i = 0; i < INIT_REF_NUM_FLOW; i++)
    mesh_flow.refine_all_elements();

  // Initialize boundary condition types and spaces with default shapesets.
  // For the concentration.
  EssentialBCs bcs_concentration;
  
  bcs_concentration.add_boundary_condition(new ConcentrationTimedepEssentialBC(BDY_DIRICHLET_CONCENTRATION, CONCENTRATION_EXT, CONCENTRATION_EXT_STARTUP_TIME));
  bcs_concentration.add_boundary_condition(new ConcentrationTimedepEssentialBC(BDY_SOLID_WALL_TOP, 0.0, CONCENTRATION_EXT_STARTUP_TIME));
  
  L2Space space_rho(&mesh_flow, P_INIT_FLOW);
  L2Space space_rho_v_x(&mesh_flow, P_INIT_FLOW);
  L2Space space_rho_v_y(&mesh_flow, P_INIT_FLOW);
  L2Space space_e(&mesh_flow, P_INIT_FLOW);
  // Space for concentration.
  H1Space space_c(&mesh_concentration, &bcs_concentration, P_INIT_CONCENTRATION);

  int ndof = Space::get_num_dofs(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  InitialSolutionEulerDensity prev_rho(&mesh_flow, RHO_EXT);
  InitialSolutionEulerDensityVelX prev_rho_v_x(&mesh_flow, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY prev_rho_v_y(&mesh_flow, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy prev_e(&mesh_flow, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
  InitialSolutionConcentration prev_c(&mesh_concentration, 0.0);

  // Numerical flux.
  OsherSolomonNumericalFlux num_flux(KAPPA);

  // Initialize weak formulation.
  EulerEquationsWeakFormSemiImplicitCoupled wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM,
    BDY_SOLID_WALL_TOP, BDY_INLET, BDY_OUTLET, BDY_NATURAL_CONCENTRATION, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c, EPSILON, (P_INIT_FLOW == 0));
  
  wf.set_time_step(time_step);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c), is_linear);

  // If the FE problem is in fact a FV problem.
  //if(P_INIT == 0) dp.set_fvm();  

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA, RHO_EXT, P_EXT);

  /*
  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));
  ScalarView s5("Concentration", new WinGeom(700, 400, 600, 300));
  */
  
  ScalarView s1("1", new WinGeom(0, 0, 600, 300));
  ScalarView s2("2", new WinGeom(700, 0, 600, 300));
  ScalarView s3("3", new WinGeom(0, 400, 600, 300));
  ScalarView s4("4", new WinGeom(700, 400, 600, 300));
  ScalarView s5("Concentration", new WinGeom(350, 200, 600, 300));

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Set up CFL calculation class.
  CFLCalculation CFL(CFL_NUMBER, KAPPA);

  // Set up Advection-Diffusion-Equation stability calculation class.
  ADEStabilityCalculation ADES(ADVECTION_STABILITY_CONSTANT, DIFFUSION_STABILITY_CONSTANT, EPSILON);

  int iteration = 0; double t = 0;
  for(t = 0.0; t < 100.0; t += time_step) {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    // Set the current time step.
    wf.set_time_step(time_step);
    Space::update_essential_bc_values(&space_c, t);

    // Assemble stiffness matrix and rhs.
    info("Assembling the stiffness matrix and right-hand side vector.");
    dp.assemble(matrix, rhs);

    // Solve the matrix problem.
    info("Solving the matrix problem.");
    scalar* solution_vector = NULL;
    if(solver->solve()) {
      solution_vector = solver->get_solution();
      Solution::vector_to_solutions(solution_vector, Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e, &space_c), Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c));
    }
    else
    error ("Matrix solver failed.\n");

    if(SHOCK_CAPTURING) {
      DiscontinuityDetector discontinuity_detector(Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));

      std::set<int> discontinuous_elements = discontinuity_detector.get_discontinuous_element_ids(DISCONTINUITY_DETECTOR_PARAM);

      FluxLimiter flux_limiter(solution_vector, Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
        &space_rho_v_y, &space_e), Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));

      flux_limiter.limit_according_to_detector(discontinuous_elements);
    }

    util_time_step = time_step;

    CFL.calculate_semi_implicit(Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), &mesh_flow, util_time_step);

    time_step = util_time_step;

    ADES.calculate(Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y), &mesh_concentration, util_time_step);

    if(util_time_step < time_step)
      time_step = util_time_step;

    // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0) {
      // Hermes visualization.
      if(HERMES_VISUALIZATION) {
        /*
        Mach_number.reinit();
        pressure.reinit();
        entropy.reinit();
        pressure_view.show(&pressure);
        entropy_production_view.show(&entropy);
        Mach_number_view.show(&Mach_number);
        s5.show(&prev_c);
        */
        s1.show(&prev_rho);
        s2.show(&prev_rho_v_x);
        s3.show(&prev_rho_v_y);
        s4.show(&prev_e);
        s5.show(&prev_c);
        /*
        s1.save_numbered_screenshot("density%i.bmp", iteration, true);
        s2.save_numbered_screenshot("density_v_x%i.bmp", iteration, true);
        s3.save_numbered_screenshot("density_v_y%i.bmp", iteration, true);
        s4.save_numbered_screenshot("energy%i.bmp", iteration, true);
        s5.save_numbered_screenshot("concentration%i.bmp", iteration, true);
        */
        //s5.wait_for_close();
        
      }
      // Output solution in VTK format.
      if(VTK_VISUALIZATION) {
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
        sprintf(filename, "Concentration-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_c, filename, "Concentration", true);
        sprintf(filename, "Concentration-3D-%i.vtk", iteration - 1);
        lin.save_solution_vtk(&prev_c, filename, "Concentration", true);
 
      }
    }
  }
  
  /*
  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();
  s5.close();
  */
  
  s1.close();
  s2.close();
  s3.close();
  s4.close();
  s5.close();

  return 0;
}
