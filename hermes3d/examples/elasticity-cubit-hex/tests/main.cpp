#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "../config.h"
#include <hermes3d.h>

// This test makes sure that the example elasticity-cubit-hex works correctly.

// The following parameters can be changed:
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int P_INIT_X = 2,
          P_INIT_Y = 2,
          P_INIT_Z = 2;                           // Initial polynomial degree of all mesh elements.
bool solution_output = true;                      // Generate output files (if true).
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters. 

const double E  = 200e9; 		          // Young modulus for steel: 200 GPa.
const double nu = 0.3;			          // Poisson ratio. 
const double f_x  = 0;   		          // x-direction load force (N). 
const double f_y  = 1e4;   		          // y-direction load force (N). 
const double f_z  = 0;   		          // z-direction load force (N). 

const double lambda = (E * nu) / ((1 + nu) * (1 - 2*nu));
const double mu = E / (2*(1 + nu));

// Boundary markers.
int bdy_fixed = 1;
int bdy_force = 2;

// Boundary condition types. 
BCType bc_types_x(int marker) 
{
  return (marker == bdy_fixed) ? H3D_BC_ESSENTIAL : H3D_BC_NATURAL;
}

BCType bc_types_y(int marker) 
{
  return (marker == bdy_fixed) ? H3D_BC_ESSENTIAL : H3D_BC_NATURAL;
}

BCType bc_types_z(int marker) 
{
  return (marker == bdy_fixed) ? H3D_BC_ESSENTIAL : H3D_BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values. 
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}

#include "../definitions.cpp"

int main(int argc, char **args) 
{
  // Test variable.
  int success_test = 1;

  // Load the mesh. 
  Mesh mesh;
  ExodusIIReader mloader;
  mloader.load("../brick_with_hole_hex.e", &mesh);

  // Perform initial mesh refinement. 
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

  // Create H1 space with default shapeset for x-displacement component. 
  H1Space xdisp(&mesh, bc_types_x, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
  
  // Create H1 space with default shapeset for y-displacement component. 
  H1Space ydisp(&mesh, bc_types_y, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
  
  // Create H1 space with default shapeset for z-displacement component. 
  H1Space zdisp(&mesh, bc_types_z, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));
  
  // Initialize weak formulation.
  WeakForm wf(3);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0), HERMES_SYM);
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1), HERMES_SYM);
  wf.add_matrix_form(0, 2, callback(bilinear_form_0_2), HERMES_SYM);
  wf.add_vector_form_surf(0, callback(surf_linear_form_x), bdy_force);

  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1), HERMES_SYM);
  wf.add_matrix_form(1, 2, callback(bilinear_form_1_2), HERMES_SYM);
  wf.add_vector_form_surf(1, callback(surf_linear_form_y), bdy_force);

  wf.add_matrix_form(2, 2, callback(bilinear_form_2_2), HERMES_SYM);
  wf.add_vector_form_surf(2, callback(surf_linear_form_z), bdy_force);

  // Initialize discrete problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&xdisp, &ydisp, &zdisp), is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize the preconditioner in the case of SOLVER_AZTECOO.
  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver)->set_solver(iterative_method);
    ((AztecOOSolver*) solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }

  // Assemble stiffness matrix and load vector.
  info("Assembling the linear problem (ndof: %d).", Space::get_num_dofs(Hermes::vector<Space *>(&xdisp, &ydisp, &zdisp)));
  dp.assemble(matrix, rhs);

  // Solve the linear system. If successful, obtain the solution.
  info("Solving the linear problem.");
  Solution xsln(xdisp.get_mesh());
  Solution ysln(ydisp.get_mesh());
  Solution zsln(zdisp.get_mesh());
  if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), 
                      Hermes::vector<Space *>(&xdisp, &ydisp, &zdisp), Hermes::vector<Solution *>(&xsln, &ysln, &zsln));
  else error ("Matrix solver failed.\n");

  // Output all components of the solution.
  if (solution_output) out_fn_vtk(&xsln, &ysln, &zsln, "sln");
  
  double sum = 0;
  for(int i = 0; i < Space::get_num_dofs(Hermes::vector<Space *>(&xdisp, &ydisp, &zdisp)); i++)
    sum += solver->get_solution()[i];

  if (abs(sum - 0.00037) / sum > 0.01)
    success_test = 0;
  // Clean up.
  delete matrix;
  delete rhs;
  delete solver;
   
  if (success_test) {
    info("Success!");
    return ERR_SUCCESS;
  }
  else {
    info("Failure!");
    return ERR_FAILURE;
  }
}
