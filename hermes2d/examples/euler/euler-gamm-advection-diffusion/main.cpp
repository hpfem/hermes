#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

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
// The parameter SETUP_VARIANT has the following meaning:
// 1 - Dirichlet condition (concentration production) on the inlet.
// 2 - Dirichlet condition (concentration production) on the bottom.
// 3 - Dirichlet condition (concentration production) on the top.
const int SETUP_VARIANT = 0;

// If set to true, GAMM channel is used, if false, a simple rectangular channel is.
const bool GAMM_CHANNEL = false;

// Visualization.
const bool HERMES_VISUALIZATION = true;               // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;                  // Set to "true" to enable VTK output.
const unsigned int EVERY_NTH_STEP = 1;                // Set visual output for every nth step.

const int P_INIT_FLOW = 0;                             // Polynomial degree for the Euler equations (for the flow).
const int P_INIT_CONCENTRATION = 1;                    // Polynomial degree for the concentration.
double CFL = 0.8;                                      // CFL value.
double time_step = 1E-4;                               // Time step.
const MatrixSolverType matrix_solver = SOLVER_UMFPACK; // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                       // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

unsigned int INIT_REF_NUM_FLOW = 2;                    // Number of initial uniform mesh refinements of the mesh for the flow.
unsigned int INIT_REF_NUM_CONCENTRATION = 3;           // Number of initial uniform mesh refinements of the mesh for the concentration.
unsigned int INIT_REF_NUM_CONCENTRATION_BDY = 5;       // Number of initial mesh refinements of the mesh for the concentration towards the 
                                                       // part of the boundary where the concentration is prescribed.
// Equation parameters.
const double P_EXT = 2.5;                              // Exterior pressure (dimensionless).
const double RHO_EXT = 1.0;                            // Inlet density (dimensionless).   
const double V1_EXT = 1.25;                            // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;                             // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;                              // Kappa.
const double CONCENTRATION_EXT = 1.0;                  // Concentration on the boundary.

const double EPSILON = 0.01;                           // Diffusivity.
// Boundary markers.
const std::string BDY_INLET = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_BOTTOM = "3";
const std::string BDY_SOLID_WALL_TOP = "4";

// Weak forms.
#include "../forms_explicit.cpp"

// Initial condition.
#include "../constant_initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh basemesh;
  H2DReader mloader;
  if(GAMM_CHANNEL)
    mloader.load("GAMM-channel-4-bnds.mesh", &basemesh);
  else
    mloader.load("channel-4-bnds.mesh", &basemesh);

  // Initialize the meshes.
  Mesh mesh_flow, mesh_concentration;
  mesh_flow.copy(&basemesh);
  mesh_concentration.copy(&basemesh);

  for(unsigned int i = 0; i < INIT_REF_NUM_CONCENTRATION; i++)
    mesh_concentration.refine_all_elements();
  switch(SETUP_VARIANT) {
  case 0:
    mesh_concentration.refine_towards_boundary(BDY_INLET, INIT_REF_NUM_CONCENTRATION_BDY);
    break;
  case 1:
    mesh_concentration.refine_towards_boundary(BDY_SOLID_WALL_BOTTOM, INIT_REF_NUM_CONCENTRATION_BDY);
    break;
  case 2:
    mesh_concentration.refine_towards_boundary(BDY_SOLID_WALL_TOP, INIT_REF_NUM_CONCENTRATION_BDY);
    break;
  }

  for(unsigned int i = 0; i < INIT_REF_NUM_FLOW; i++)
    mesh_flow.refine_all_elements();
  
  // Initialize boundary condition types and spaces with default shapesets.
  // For the flow.
  EssentialBCs bcs_flow;
  
  // For the concentration.
  EssentialBCs bcs_concentration;
  
  switch(SETUP_VARIANT) {
  case 0:
    //bcs_concentration.add_boundary_condition(new NaturalEssentialBC(Hermes::vector<std::string>(BDY_OUTLET, BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP)));
    bcs_concentration.add_boundary_condition(new DefaultEssentialBCConst(BDY_INLET, CONCENTRATION_EXT));
    break;
  case 1:
    //bcs_concentration.add_boundary_condition(new NaturalEssentialBC(Hermes::vector<std::string>(BDY_OUTLET, BDY_INLET, BDY_SOLID_WALL_TOP)));
    bcs_concentration.add_boundary_condition(new DefaultEssentialBCConst(BDY_SOLID_WALL_BOTTOM, CONCENTRATION_EXT));
    break;
  case 2:
    //bcs_concentration.add_boundary_condition(new NaturalEssentialBC(Hermes::vector<std::string>(BDY_OUTLET, BDY_SOLID_WALL_BOTTOM, BDY_INLET)));
    bcs_concentration.add_boundary_condition(new DefaultEssentialBCConst(BDY_SOLID_WALL_TOP, CONCENTRATION_EXT));
    break;
  }

  L2Space space_rho(&mesh_flow, &bcs_flow, P_INIT_FLOW);
  L2Space space_rho_v_x(&mesh_flow, &bcs_flow, P_INIT_FLOW);
  L2Space space_rho_v_y(&mesh_flow, &bcs_flow, P_INIT_FLOW);
  L2Space space_e(&mesh_flow, &bcs_flow, P_INIT_FLOW);
  // Space for concentration.
  H1Space space_c(&mesh_concentration, &bcs_concentration, P_INIT_CONCENTRATION);

  // Initialize solutions, set initial conditions.
  InitialSolutionEulerDensity sln_rho(&mesh_flow, RHO_EXT);
  InitialSolutionEulerDensityVelX sln_rho_v_x(&mesh_flow, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY sln_rho_v_y(&mesh_flow, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy sln_e(&mesh_flow, calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
  InitialSolutionConcentration sln_c(&mesh_concentration, 0.0);

  InitialSolutionEulerDensity prev_rho(&mesh_flow, RHO_EXT);
  InitialSolutionEulerDensityVelX prev_rho_v_x(&mesh_flow, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY prev_rho_v_y(&mesh_flow, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy prev_e(&mesh_flow, calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
  InitialSolutionConcentration prev_c(&mesh_concentration, 0.0);

  // Initialize weak formulation.
  EulerEquationsWeakFormExplicitCoupled wf(SETUP_VARIANT, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP, 
    BDY_INLET, BDY_OUTLET, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, &prev_c, EPSILON);
  wf.set_time_step(time_step);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e, &space_c), is_linear);

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA, RHO_EXT, P_EXT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));

  ScalarView s5("Concentration", new WinGeom(700, 400, 600, 300));
  /*
  ScalarView s1("1", new WinGeom(0, 0, 600, 300));
  ScalarView s2("2", new WinGeom(700, 0, 600, 300));
  ScalarView s3("3", new WinGeom(0, 400, 600, 300));
  ScalarView s4("4", new WinGeom(700, 400, 600, 300));
  ScalarView s5("Concentration", new WinGeom(350, 200, 600, 300));
  */

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  int iteration = 0; double t = 0;
  for(t = 0.0; t < 3.0; t += time_step) {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    bool rhs_only = (iteration == 1 ? false : true);
    // Assemble stiffness matrix and rhs or just rhs.
    if (rhs_only == false) info("Assembling the stiffness matrix and right-hand side vector.");
    else info("Assembling the right-hand side vector (only).");
    dp.assemble(matrix, rhs, rhs_only);
        
    // Solve the matrix problem.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution::vector_to_solutions(solver->get_solution(), Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e, &space_c), Hermes::vector<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e, &sln_c));
    else
    error ("Matrix solver failed.\n");

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&sln_rho);
    prev_rho_v_x.copy(&sln_rho_v_x);
    prev_rho_v_y.copy(&sln_rho_v_y);
    prev_e.copy(&sln_e);
    prev_c.copy(&sln_c);

    // Visualization.
    if((iteration - 1) % EVERY_NTH_STEP == 0) {
      // Hermes visualization.
      if(HERMES_VISUALIZATION) {
        Mach_number.reinit();
        pressure.reinit();
        entropy.reinit();
        pressure_view.show(&pressure);
        entropy_production_view.show(&entropy);
        Mach_number_view.show(&Mach_number);
        s5.show(&prev_c);
        /*
        s1.show(&prev_rho);
        s2.show(&prev_rho_v_x);
        s3.show(&prev_rho_v_y);
        s4.show(&prev_e);
        s5.show(&prev_c);
        */
      }
      // Output solution in VTK format.
      if(VTK_VISUALIZATION) {
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
  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();
  s5.close();

  /*
  s1.close();
  s2.close();
  s3.close();
  s4.close();
  */

  return 0;
}
