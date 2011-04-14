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
const int INIT_REF_NUM = 4;                       // Number of initial uniform mesh refinements.                       
const int P_INIT = 0;                             // Initial polynomial degree.                      
double CFL = 0.8;                                 // CFL value.
double time_step = 1E-4;                                // Time step.
const MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
// Equation parameters.
const double P_EXT = 2.5;         // Exterior pressure (dimensionless).
const double RHO_EXT = 1.0;       // Inlet density (dimensionless).   
const double V1_EXT = 1.25;       // Inlet x-velocity (dimensionless).
const double V2_EXT = 0.0;        // Inlet y-velocity (dimensionless).
const double KAPPA = 1.4;         // Kappa.

// Boundary markers.
const std::string BDY_INLET = "1";
const std::string BDY_OUTLET = "2";
const std::string BDY_SOLID_WALL_BOTTOM = "3";
const std::string BDY_SOLID_WALL_TOP = "4";

// Weak forms.
#include "../../forms_explicit.cpp"

// Initial condition.
#include "../../constant_initial_condition.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../channel.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(2);
//  mesh.refine_towards_boundary(BDY_SOLID_WALL_BOTTOM, 2);

  // Initialize boundary condition types and spaces with default shapesets.
  L2Space space_rho(&mesh, P_INIT);
  L2Space space_rho_v_x(&mesh, P_INIT);
  L2Space space_rho_v_y(&mesh, P_INIT);
  L2Space space_e(&mesh, P_INIT);

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
  EulerEquationsWeakFormExplicitMultiComponentSemiImplicit wf(&num_flux, KAPPA, RHO_EXT, V1_EXT, V2_EXT, P_EXT, BDY_SOLID_WALL_BOTTOM, BDY_SOLID_WALL_TOP, 
    BDY_INLET, BDY_OUTLET, &prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, Hermes::vector<Space*>(&space_rho, &space_rho_v_x, &space_rho_v_y, &space_e), is_linear);
  
  // If the FE problem is in fact a FV problem.
  if(P_INIT == 0)
    dp.set_fvm();
  
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

  for(unsigned int time_iteration = 0; time_iteration < 5; time_iteration++)
  {
    bool rhs_only = (time_iteration == 0 ? false : true);

    if (rhs_only == false) {
      info("Assembling the stiffness matrix and right-hand side vector.");
      dp.assemble(matrix, rhs);
    }
    else {
      info("Assembling the right-hand side vector (only).");
      dp.assemble(NULL, rhs);
    }

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

    // Storing the testing values.
    for(unsigned int j = 0; j < 4; j++)
      for(unsigned int k = j*space_rho.get_num_dofs(); k < (j+1)*space_rho.get_num_dofs(); k++)
        l2_norms[time_iteration][j] += solver->get_solution()[k];    

    point_values[time_iteration][0] = sln_rho_v_x.get_pt_value(0.5, 0.001);
    point_values[time_iteration][1] = sln_rho_v_x.get_pt_value(x, y);
    point_values[time_iteration][2] = sln_rho_v_x.get_pt_value(1.5, 0.001);

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&sln_rho);
    prev_rho_v_x.copy(&sln_rho_v_x);
    prev_rho_v_y.copy(&sln_rho_v_y);
    prev_e.copy(&sln_e);
  }
  bool okay = true;
  switch(P_INIT) {
  case 0:
    if(fabs(l2_norms[0][0] - 16) > 1E-8)
      okay = false;
    if(fabs(l2_norms[0][1] - 20) > 1E-8)
      okay = false;
    if(fabs(l2_norms[0][2]) > 1E-8)
      okay = false;
    if(fabs(l2_norms[0][3] - 112.5) > 1E-8)
      okay = false;
    if(fabs(l2_norms[1][0] - 16) > 1E-8)
      okay = false;
    if(fabs(l2_norms[1][1] - 20) > 1E-8)
      okay = false;
    if(fabs(l2_norms[1][2]) > 1E-8)
      okay = false;
    if(fabs(l2_norms[1][3] - 112.5) > 1E-8)
      okay = false;
    if(fabs(l2_norms[2][0] - 16) > 1E-8)
      okay = false;
    if(fabs(l2_norms[2][1] - 20) > 1E-8)
      okay = false;
    if(fabs(l2_norms[2][2]) > 1E-8)
      okay = false;
    if(fabs(l2_norms[2][3] - 112.5) > 1E-8)
      okay = false;
    if(fabs(l2_norms[3][0] - 16) > 1E-8)
      okay = false;
    if(fabs(l2_norms[3][1] - 20) > 1E-8)
      okay = false;
    if(fabs(l2_norms[3][2]) > 1E-8)
      okay = false;
    if(fabs(l2_norms[3][3] - 112.5) > 1E-8)
      okay = false;
    if(fabs(l2_norms[4][0] - 16) > 1E-8)
      okay = false;
    if(fabs(l2_norms[4][1] - 20) > 1E-8)
      okay = false;
    if(fabs(l2_norms[4][2]) > 1E-8)
      okay = false;
    if(fabs(l2_norms[4][3] - 112.5) > 1E-8)
      okay = false;

    // points
    if(fabs(point_values[0][0] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[0][1] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[0][2] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[1][0] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[1][1] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[1][2] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[2][0] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[2][1] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[2][2] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[3][0] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[3][1] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[3][2] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[4][0] - 1.25)> 1E-8)
      okay = false;
    if(fabs(point_values[4][1] - 1.25) > 1E-8)
      okay = false;
    if(fabs(point_values[4][2] - 1.25) > 1E-8)
      okay = false;
    break;
  }

  for (int i = 0; i<5; i++)
    for (int j=0; j<4; j++)
      info ("l2_norms[%d][%d] = %g", i, j, l2_norms[i][j]);

  for (int i = 0; i<5; i++)
    for (int j=0; j<3; j++)
      info ("point[%d][%d] = %g", i, j, point_values[i][j]);

  if (okay) {      // ndofs was 908 at the time this test was created
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  } 
}
