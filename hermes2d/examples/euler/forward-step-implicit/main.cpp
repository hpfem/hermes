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
// Domain: forward facing step, see mesh file ffs.mesh
//
// BC: Normal velocity component is zero on solid walls.
//     Full supersonic state prescribed at inlet.
//     Pressure given at outlet, but used only if outlet flow is subsonic/
//
// IC: Constant supersonic state identical to inlet. 
//
// The following parameters can be changed:
// Visualization.
const bool HERMES_VISUALIZATION = true;           // Set to "true" to enable Hermes OpenGL visualization. 
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const unsigned int EVERY_NTH_STEP = 1;            // Set visual output for every nth step.

// Use of preconditioning.
const bool PRECONDITIONING = true;
const double NOX_LINEAR_TOLERANCE = 1e-2;

const int P_INIT = 0;                                   // Initial polynomial degree.                      
const int INIT_REF_NUM = 3;                             // Number of initial uniform mesh refinements.                       
double time_step = 1E-2;                                // Time step.

// Equation parameters.
const double P_EXT = 1.0;         // Exterior pressure (dimensionless).
const double RHO_EXT = 1.4;       // Inlet density (dimensionless).   
const double V1_EXT = 3.0;        // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;                              // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;                               // Kappa.

// Boundary markers.
const std::string BDY_SOLID_WALL = "1";
const std::string BDY_INLET_OUTLET = "2";

// Weak forms.
#include "../forms_implicit.cpp"

// Initial condition.
#include "../initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("ffs.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) 
    mesh.refine_all_elements();

  // Initialize boundary condition types and spaces with default shapesets.
  L2Space space_rho(&mesh, P_INIT);
  L2Space space_rho_v_x(&mesh, P_INIT);
  L2Space space_rho_v_y(&mesh, P_INIT);
  L2Space space_e(&mesh, P_INIT);
  int ndof = Space::get_num_dofs(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  info("ndof: %d", ndof);

  // Initialize solutions, set initial conditions.
  InitialSolutionEulerDensity prev_rho(&mesh, RHO_EXT);
  InitialSolutionEulerDensityVelX prev_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY prev_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy prev_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  // Solutions for the time derivative estimate.
  Solution sln_temp_rho, sln_temp_rho_v_x, sln_temp_rho_v_y, sln_temp_e;

  // Numerical flux.
  OsherSolomonNumericalFlux num_flux(KAPPA); 

  // Initialize weak formulation.
  EulerEquationsWeakFormImplicit wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL, BDY_SOLID_WALL, 
    BDY_INLET_OUTLET, BDY_INLET_OUTLET, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e, PRECONDITIONING);
  wf.set_time_step(time_step);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e));
  
  // If the FE problem is in fact a FV problem.
  if(P_INIT == 0)
    dp.set_fvm();

  // Project the initial solution on the FE space 
  // in order to obtain initial vector for NOX. 
  info("Projecting initial solution on the FE mesh.");
  scalar* coeff_vec = new scalar[Space::get_num_dofs(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e))];
  OGProjection::project_global(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), 
    Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec);

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), KAPPA, RHO_EXT, P_EXT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));

  /*
  ScalarView s1("1", new WinGeom(0, 0, 600, 300));
  ScalarView s2("2", new WinGeom(700, 0, 600, 300));
  ScalarView s3("3", new WinGeom(0, 400, 600, 300));
  ScalarView s4("4", new WinGeom(700, 400, 600, 300));
  */
  
  // Initialize NOX solver.
  NoxSolver solver(&dp);
  solver.set_ls_tolerance(NOX_LINEAR_TOLERANCE);
  solver.disable_abs_resid();
  solver.set_conv_rel_resid(1.00);

  // Select preconditioner.
  RCP<Precond> pc = rcp(new MlPrecond("sa"));
  solver.set_precond(pc);

  int iteration = 0; double t = 0;
  for(t = 0.0; t < 3.0; t += time_step)
  {
    info("---- Time step %d, time %3.5f.", iteration++, t);

    OGProjection::project_global(Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), 
    Hermes::vector<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), coeff_vec);

    info("Assembling by DiscreteProblem, solving by NOX.");
    solver.set_init_sln(coeff_vec);
    if (solver.solve())
      Solution::vector_to_solutions(solver.get_solution(), Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), 
      Hermes::vector<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    else
      error("NOX failed.");
   
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
        /*
        s1.show(&prev_rho);
        s2.show(&prev_rho_v_x);
        s3.show(&prev_rho_v_y);
        s4.show(&prev_e);
        */
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
      }
    }

    info("Number of nonlin iterations: %d (norm of residual: %g)", 
      solver.get_num_iters(), solver.get_residual());
    info("Total number of iterations in linsolver: %d (achieved tolerance in the last step: %g)", 
      solver.get_num_lin_iters(), solver.get_achieved_tol());
  }
  
  pressure_view.close();
  entropy_production_view.close();
  Mach_number_view.close();

  /*
  s1.close();
  s2.close();
  s3.close();
  s4.close();
  */
  
  return 0;
}
