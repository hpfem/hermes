#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This example solves the compressible Euler equations using a basic
// piecewise-constant finite volume method.
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
const int P_INIT = 0;                             // Initial polynomial degree.                      
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.                       
double CFL = 0.8;                                 // CFL value.
double time_step = 1E-4;                                // Time step.
const MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Equation parameters.
const double P_EXT = 1.0;         // Exterior pressure (dimensionless).
const double RHO_EXT = 1.4;       // Inlet density (dimensionless).   
const double V1_EXT = 3.0;        // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;        // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;         // Kappa.

// Boundary markers.
const std::string BDY_SOLID_WALL_BOTTOM = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_TOP = "3";
const std::string BDY_INLET = "4";

// Weak forms.
#include "../forms_explicit.cpp"

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
  InitialSolutionEulerDensity sln_rho(&mesh, RHO_EXT);
  InitialSolutionEulerDensityVelX sln_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY sln_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy sln_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));
  
  InitialSolutionEulerDensity prev_rho(&mesh, RHO_EXT);
  InitialSolutionEulerDensityVelX prev_rho_v_x(&mesh, RHO_EXT * V1_EXT);
  InitialSolutionEulerDensityVelY prev_rho_v_y(&mesh, RHO_EXT * V2_EXT);
  InitialSolutionEulerDensityEnergy prev_e(&mesh, QuantityCalculator::calc_energy(RHO_EXT, RHO_EXT * V1_EXT, RHO_EXT * V2_EXT, P_EXT, KAPPA));

  // Numerical flux.
  OsherSolomonNumericalFlux num_flux(KAPPA); 

  // Initialize weak formulation.
  EulerEquationsWeakFormExplicit wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP, 
    BDY_INLET, BDY_OUTLET, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);

  // Initialize the FE problem.
  bool is_linear = true;
  
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), is_linear);
  
  // If the FE problem is in fact a FV problem.
  if(P_INIT == 0)
    dp.set_fvm();  

  // Filters for visualization of Mach number, pressure and entropy.
  MachNumberFilter Mach_number(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA);
  PressureFilter pressure(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA);
  EntropyFilter entropy(Hermes::vector<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), KAPPA, RHO_EXT, P_EXT);

  ScalarView pressure_view("Pressure", new WinGeom(0, 0, 600, 300));
  ScalarView Mach_number_view("Mach number", new WinGeom(700, 0, 600, 300));
  ScalarView entropy_production_view("Entropy estimate", new WinGeom(0, 400, 600, 300));

  /*
  ScalarView s1("1", new WinGeom(0, 0, 600, 300));
  ScalarView s2("2", new WinGeom(700, 0, 600, 300));
  ScalarView s3("3", new WinGeom(0, 400, 600, 300));
  ScalarView s4("4", new WinGeom(700, 400, 600, 300));
  */
  
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  int iteration = 0; double t = 0;
  for(t = 0.0; t < 10.0; t += time_step) {
    info("---- Time step %d, time %3.5f.", iteration++, t);
    
    bool rhs_only = (iteration == 1 ? false : true);
    // Assemble stiffness matrix and rhs or just rhs.
    if (rhs_only == false) info("Assembling the stiffness matrix and right-hand side vector.");
    else info("Assembling the right-hand side vector (only).");
    // Set the current time step.
    wf.set_time_step(time_step);
    dp.assemble(matrix, rhs, rhs_only);

    // Solve the matrix problem.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution::vector_to_solutions(solver->get_solution(), Hermes::vector<Space *>(&space_rho, &space_rho_v_x, 
      &space_rho_v_y, &space_e), Hermes::vector<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
    else
    error ("Matrix solver failed.\n");

    // Determine the time step according to the CFL condition.
    // Only mean values on an element of each solution component are taken into account.
    double *solution_vector = solver->get_solution();
    double min_condition = 0;
    Element *e;
    for (int _id = 0, _max = mesh.get_max_element_id(); _id < _max; _id++) \
          if (((e) = mesh.get_element_fast(_id))->used) \
            if ((e)->active) {
      AsmList al;
      space_rho.get_element_assembly_list(e, &al);
      double rho = solution_vector[al.dof[0]];
      space_rho_v_x.get_element_assembly_list(e, &al);
      double v1 = solution_vector[al.dof[0]] / rho;
      space_rho_v_y.get_element_assembly_list(e, &al);
      double v2 = solution_vector[al.dof[0]] / rho;
      space_e.get_element_assembly_list(e, &al);
      double energy = solution_vector[al.dof[0]];
      
              double condition = e->get_area() / (std::sqrt(v1*v1 + v2*v2) + QuantityCalculator::calc_sound_speed(rho, rho*v1, rho*v2, energy, KAPPA));
      
      if(condition < min_condition || min_condition == 0.)
        min_condition = condition;
    }
    if(time_step > min_condition)
      time_step = min_condition;
    if(time_step < min_condition * 0.9)
      time_step = min_condition;

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&sln_rho);
    prev_rho_v_x.copy(&sln_rho_v_x);
    prev_rho_v_y.copy(&sln_rho_v_y);
    prev_e.copy(&sln_e);

    // Visualization.
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