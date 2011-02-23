#define HERMES_REPORT_INFO
#include "hermes2d.h"

// This test makes sure that example 37-nurbs works correctly.

const char* mesh_file_1 = "domain-1.mesh";          // One control point.
const char* mesh_file_2 = "domain-2.mesh";          // Two control points.
const char* mesh_file_3 = "domain-3.mesh";          // Three control points.

// The following parameters can be also changed:

const int P_INIT = 3;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const int BDY_MARKER = 1;

// Problem parameters.
const double CONST_F = 1.0;  

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Check number of command-line parameters.
  if (argc < 2) {
    error("Not enough parameters.");
  }

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  if (strcasecmp(argv[1], "1") == 0)
    mloader.load(mesh_file_1, &mesh);
  if (strcasecmp(argv[1], "2") == 0)
    mloader.load(mesh_file_2, &mesh);
  if (strcasecmp(argv[1], "3") == 0)
    mloader.load(mesh_file_3, &mesh);

  // Perform initial mesh refinements (optional).
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_MARKER);

  // Enter zero Dirichlet boundary values.
  BCValues bc_values;

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

   // Initialize the solution.
  Solution sln;

  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling the stiffness matrix and right-hand side vector.");
  dp.assemble(matrix, rhs);

  // Solve the linear system and if successful, obtain the solution.
  info("Solving the matrix problem.");
  if(solver->solve())
    Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else
    error ("Matrix solver failed.\n");

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  info("Coordinate ( 0.3, 0.5) value = %lf", sln.get_pt_value(0.3, 0.5));
  info("Coordinate ( 0.7, 0.5) value = %lf", sln.get_pt_value(0.7, 0.5));
  info("Coordinate ( 1.3, 0.5) value = %lf", sln.get_pt_value(1.3, 0.5));
  info("Coordinate ( 1.7, 0.5) value = %lf", sln.get_pt_value(1.7, 0.5));

  double coor_x[4] = {0.3, 0.7, 1.3, 1.7};
  double coor_y = 0.5;

  double value[4] = {0.102569, 0.167907, 0.174203, 0.109630};
  if (strcasecmp(argv[1], "2") == 0)
  {
    value[0] = 0.062896; 
    value[1] = 0.096658;
    value[2] = 0.114445;
    value[3] = 0.081221;
  }
  if (strcasecmp(argv[1], "3") == 0)
  {
    value[0] = 0.048752; 
    value[1] = 0.028585;
    value[2] = 0.028585;
    value[3] = 0.048752;
  }

  bool success = true;
  for (int i = 0; i < 4; i++)
  {
    if (abs(value[i] - sln.get_pt_value(coor_x[i], coor_y)) > 1E-6) success = false;
  }
  if (success) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

