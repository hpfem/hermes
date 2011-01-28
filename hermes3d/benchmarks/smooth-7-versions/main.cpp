#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
#include <hermes3d.h>

//  This benchmark is an extension of the H2D benchmarks smooth-iso, smooth-aniso-x
//  and smooth-aniso-y to 3D. Solved is the Poisson equation in a cube and the exact 
//  solutions are sin(x), sin(y), sin(z), sin(x)*sin(y), sin(x)*sin(z), sin(y)*sin(z), 
//  and sin(x)*sin(y)*sin(z). According to this, the executable needs to be called with 
//  'x', 'y', 'z', 'xy', 'xz', 'yz' or 'xyz' as a command-line parameter. 
//
//  PDE: -Laplace u = f.
//
//  Known exact solution, see functions fn() and fndd().
//
//  Domain: cube (0, pi) x (0, pi) x (0, pi), mesh file hex-0-1.mesh3d.
//
//  BC:  Dirichlet and Neumann given by exact solution, see the function bc_types() below.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
int P_INIT_X, P_INIT_Y, P_INIT_Z;                 // Initial polynomial degree of all mesh elements, to be 
                                                  // determined according to the anisotropy. See bc_types().
const double THRESHOLD = 0.3;                     // Error threshold for element refinement of the adapt(...) function 
                                                  // (default) STRATEGY = 0 ... refine elements elements until sqrt(THRESHOLD) 
                                                  // times total error is processed. If more elements have similar errors, 
                                                  // refine all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  // than THRESHOLD times maximum element error.
const double ERR_STOP = 1e-4;                     // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                     // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
bool solution_output = true;                      // Generate output files (if true).
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Exact solution
#define ANISO_X	   1     
#define ANISO_Y	   2
#define ANISO_Z	   4
int ANISO_TYPE;
#include "exact_solution.cpp"

// Boundary condition types. We also assign directional polynomial degrees here. 
BCType bc_types(int marker)
{
  switch (ANISO_TYPE) {
    // u = sin(x), thus we prescribe zero Dirichlet at faces 1 and 2, and zero Neumann elsewhere.
    case ANISO_X: if (marker == 1 || marker == 2) return BC_ESSENTIAL; 
                  else return BC_NATURAL;
    // u = sin(y), thus we prescribe zero Dirichlet at faces 3 and 4, and zero Neumann elsewhere.
    case ANISO_Y: if (marker == 3 || marker == 4) return BC_ESSENTIAL; 
                  else return BC_NATURAL;
    // u = sin(z), thus we prescribe zero Dirichlet at faces 5 and 6, and zero Neumann elsewhere.
    case ANISO_Z: if (marker == 5 || marker == 6) return BC_ESSENTIAL; 
                  else return BC_NATURAL;
    // u = sin(x)*sin(y), thus we prescribe zero Dirichlet at faces 1, 2, 3, 4 and zero Neumann elsewhere.
    case ANISO_X | ANISO_Y: if (marker == 1 || marker == 2 || marker == 3 || marker == 4) return BC_ESSENTIAL; 
                            else return BC_NATURAL;
    // u = sin(x)*sin(z), thus we prescribe zero Dirichlet at faces 1, 2, 5, 6 and zero Neumann elsewhere.
    case ANISO_X | ANISO_Z: if (marker == 1 || marker == 2 || marker == 5 || marker == 6) return BC_ESSENTIAL; 
                            else return BC_NATURAL;
    // u = sin(y)*sin(z), thus we prescribe zero Dirichlet at faces 3, 4, 5, 6 and zero Neumann elsewhere.
    case ANISO_Y | ANISO_Z: if (marker == 3 || marker == 4 || marker == 5 || marker == 6) return BC_ESSENTIAL; 
                            else return BC_NATURAL;
    // u = sin(x)*sin(y)*sin(z), thus we prescribe zero Dirichlet everywhere.
    case ANISO_X | ANISO_Y | ANISO_Z: return BC_ESSENTIAL;
    default: error("Admissible command-line options are x, y, x, xy, xz, yz, xyz.");
  }

  return BC_ESSENTIAL;
}

// Assign the lowest possible directional polynomial degrees so that the problem's NDOF >= 1.
void assign_poly_degrees()
{
  switch (ANISO_TYPE) {
    // u = sin(x), thus we prescribe zero Dirichlet at faces 1 and 2, and zero Neumann elsewhere.
    case ANISO_X: P_INIT_X = 2; P_INIT_Y = 1; P_INIT_Z = 1; break;
    // u = sin(y), thus we prescribe zero Dirichlet at faces 3 and 4, and zero Neumann elsewhere.
    case ANISO_Y: P_INIT_X = 1; P_INIT_Y = 2; P_INIT_Z = 1; break;
    // u = sin(z), thus we prescribe zero Dirichlet at faces 5 and 6, and zero Neumann elsewhere.
    case ANISO_Z: P_INIT_X = 1; P_INIT_Y = 1; P_INIT_Z = 2; break; 
    // u = sin(x)*sin(y), thus we prescribe zero Dirichlet at faces 1, 2, 3, 4 and zero Neumann elsewhere.
    case ANISO_X | ANISO_Y: P_INIT_X = 2; P_INIT_Y = 2; P_INIT_Z = 1; break;
    // u = sin(x)*sin(z), thus we prescribe zero Dirichlet at faces 1, 2, 5, 6 and zero Neumann elsewhere.
    case ANISO_X | ANISO_Z: P_INIT_X = 2; P_INIT_Y = 1; P_INIT_Z = 2; break;
    // u = sin(y)*sin(z), thus we prescribe zero Dirichlet at faces 3, 4, 5, 6 and zero Neumann elsewhere.
    case ANISO_Y | ANISO_Z: P_INIT_X = 1; P_INIT_Y = 2; P_INIT_Z = 2; break;
    // u = sin(x)*sin(y)*sin(z), thus we prescribe zero Dirichlet everywhere.
    case ANISO_X | ANISO_Y | ANISO_Z: P_INIT_X = 2; P_INIT_Y = 2; P_INIT_Z = 2; break;
    default: error("Admissible command-line options are x, y, x, xy, xz, yz, xyz.");
  }
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return fn(x, y, z);
}

// Determine type of anisotropy from the command-line parameter.
int parse_aniso_type(char *str)
{
  int type = 0;
  if (strchr(str, 'x') != NULL) type |= ANISO_X;
  if (strchr(str, 'y') != NULL) type |= ANISO_Y;
  if (strchr(str, 'z') != NULL) type |= ANISO_Z;
  return type;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char **args)
{
  // Check the number of command-line parameters.
  if (argc < 2) {
    info("Use x, y, z, xy, xz, yz, or xyz as a command-line parameter.");
    error("Not enough command-line parameters.");
  }

  // Determine anisotropy type from the command-line parameter.
  ANISO_TYPE = parse_aniso_type(args[1]);

  // Load the mesh.
  Mesh mesh;
  H3DReader mesh_loader;
  mesh_loader.load("hex-0-1.mesh3d", &mesh);

  // Assign the lowest possible directional polynomial degrees so that the problem's NDOF >= 1.
  assign_poly_degrees();

  // Create an H1 space with default shapeset.
  info("Setting directional polynomial degrees %d, %d, %d.", P_INIT_X, P_INIT_Y, P_INIT_Z);
  H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM, HERMES_ANY);
  wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY);

  // Set exact solution.
  ExactSolution exact(&mesh, fndd);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, graph_dof_exact, graph_cpu_exact;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity loop. 
  int as = 1; 
  bool done = false;
  do 
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = construct_refined_space(&space, 1);

    // Initialize discrete problem.
    bool is_linear = true;
    DiscreteProblem dp(&wf, ref_space, is_linear);

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
  
    // Assemble the reference problem.
    info("Assembling on reference mesh (ndof: %d).", Space::get_num_dofs(ref_space));
    dp.assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();

    // Solve the linear system on reference mesh. If successful, obtain the solution.
    info("Solving on reference mesh.");
    Solution ref_sln(ref_space->get_mesh());
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), ref_space, &ref_sln);
    else error ("Matrix solver failed.\n");

    // Time measurement.
    cpu_time.tick();

    // Project the reference solution on the coarse mesh.
    Solution sln(space.get_mesh());
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver);

    // Time measurement.
    cpu_time.tick();

    // Output solution and mesh with polynomial orders.
    if (solution_output) 
    {
      out_fn_vtk(&sln, "sln", as);
      out_orders_vtk(&space, "order", as);
    }

    // Skip the visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate and exact error.");
    Adapt *adaptivity = new Adapt(&space, HERMES_H1_NORM);
    bool solutions_for_adapt = true;
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt) * 100;

    // Calculate exact error.
    solutions_for_adapt = false;
    double err_exact_rel = adaptivity->calc_err_exact(&sln, &exact, solutions_for_adapt) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d.", Space::get_num_dofs(&space), Space::get_num_dofs(ref_space));
    info("err_est_rel: %g%%, err_exact_rel: %g%%.", err_est_rel, err_exact_rel);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");
    graph_dof_exact.add_values(Space::get_num_dofs(&space), err_exact_rel);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel);
    graph_cpu_exact.save("conv_cpu_exact.dat");

    // If err_est_rel is too large, adapt the mesh. 
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      adaptivity->adapt(THRESHOLD);
    }
    if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;

    // Clean up.
    delete ref_space->get_mesh();
    delete ref_space;
    delete matrix;
    delete rhs;
    delete solver;
    delete adaptivity;

    // Increase the counter of performed adaptivity steps.
    as++;
  } while (!done);

  return 0;
}
