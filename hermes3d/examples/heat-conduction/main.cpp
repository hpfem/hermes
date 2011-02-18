#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
#include <hermes3d.h>

//  This example shows how to solve a time-dependent PDE discretized
//  in time via the implicit Euler method on a fixed mesh. 
//  You will see how to use the solution from previous time step.
//
//  PDE: stationary heat transfer equation
//  dT/dt - Laplace T = f.
//
//  Domain: (0, 0, 1)x(0, 1, 0)x(1, 0, 0), see the file hexahedron.mesh3d. 
//
//  Known exact solution, see unctions fn() and fndd(). 
//
//  IC:  T = 0.
//  BC:  T = 0. ... Dirichlet,
//
//  Time-stepping: implicit Euler.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 2;			  // Number of initial uniform mesh refinements.
const int P_INIT_X = 2,
          P_INIT_Y = 2,
          P_INIT_Z = 2;                           // Initial polynomial degree of all mesh elements.
const double TAU = 0.05;			  // Time step in seconds. 
bool solution_output = true;                      // Generate output files (if true).
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h)   
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.                                               

// Problem parameters. 
const double FINAL_TIME = 2 * M_PI;		  // Length of time interval in seconds. 

// Global time variable. 
double TIME = TAU;

// Exact solution. 
#include "exact_solution.cpp"

// Boundary condition types.
BCType bc_types(int marker) {
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values. 
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}

#include "forms.cpp"

int main(int argc, char **args) 
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the initial mesh. 
  Mesh mesh;
  H3DReader mesh_loader;
  mesh_loader.load("hexahedron.mesh3d", &mesh);

  // Perform initial mesh refinement. 
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

  // Create H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));

  // Construct initial solution and set it to zero.
  Solution sln_prev(&mesh);
  sln_prev.set_zero();

  // Initialize weak formulation. 
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
  wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY, &sln_prev);

  // Initialize discrete problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

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

  // Output mesh with polynomial orders.
  if (solution_output)
    out_orders_vtk(&space, "order");

  // Time stepping. 
  int nsteps = (int) (FINAL_TIME/TAU + 0.5);
  for (int ts = 0; ts < nsteps;  ts++)
  {
    info("---- Time step %d, time %3.5f.", ts, TIME);
   
    // Assemble the linear problem.
    info("Assembling the linear problem (ndof: %d).", Space::get_num_dofs(&space));

    bool rhsonly = (ts > 0);
    dp.assemble(matrix, rhs, rhsonly);

    // Solve the linear system. If successful, obtain the solution.
    info("Solving the linear problem.");
    Solution sln(space.get_mesh());
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
    else error ("Matrix solver failed.\n");

    // Output solution.
    if (solution_output)
      out_fn_vtk(&sln, "sln", ts);

    // Calculate error wrt. exact solution. 
    info("Calculating exact error.");
    ExactSolution esln(&mesh, fndd);
    double err_exact = h1_error(&sln, &esln) * 100; 
    info("Err. exact: %g%%.", err_exact);

    // Next time step.
    sln_prev = sln;
    TIME += TAU;
  }

  // Time measurement.
  cpu_time.tick();

  // Print timing information.
  info("Solutions and mesh with polynomial orders saved. Total running time: %g s", cpu_time.accumulated());

  // Clean up.
  delete matrix;
  delete rhs;
  delete solver;

  return 0;
}
