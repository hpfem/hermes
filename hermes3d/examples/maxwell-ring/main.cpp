#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
#include <hermes3d.h>

//  This example describes the diffraction
//  of an electromagnetic wave from a re-entrant corner. Convergence graphs are 
//  saved, both exact error and error estimate, and both in terms of DOF and CPU time.
//
//  PDE: time-harmonic Maxwell's equations.
//
//  Known exact solution, see functions exact_sol_val(), exact_sol(), exact().
//
//  Domain: L-shape 3D domain.
//
//  Meshes: "lshape_hex.mesh3d" (hexahedron mesh) See the mesh.load(...) command below.
//
//  BC: perfect conductor on boundary markers 1 and 6 (essential BC),
//      impedance boundary condition on rest of boundary (natural BC).
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int P_INIT_X = 1,
          P_INIT_Y = 1,
          P_INIT_Z = 1;                           // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.7;                     // Error threshold for element refinement of the adapt(...) function 
                                                  // (default) STRATEGY = 0 ... refine elements elements until sqrt(THRESHOLD) 
                                                  // times total error is processed. If more elements have similar errors, 
                                                  // refine all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  // than THRESHOLD times maximum element error.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
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

// Problem parameters.
const double mu_0   = 4*M_PI*1e-7;
const double mu_r   = 1.0;
const double kappa  = 1.0;
const double lambda = 1.0;
const double OMEGA = 2*3.14*1000;
const double SIGMA = 57*10^6;
const double J[3] = {1e+7 *  1 / sqrt((double)2), 0, 1e+7 *   1/ sqrt((double)2)};
//const double J[3] = {1e+7 *  z / sqrt(x^2 + z^2), 0, 1e+7 *  z / sqrt(x^2 + z^2)};


// Weak forms.
#include "forms.cpp"

// Boundary condition types. 
BCType bc_types(int marker)
{
  if (marker == 1 || marker == 6)
    return BC_ESSENTIAL; // perfect conductor
  else
    return BC_NATURAL; // impedance
}

// Essential (Dirichlet) boundary condition values. 
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}

int main(int argc, char **args) 
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh. 
  Mesh mesh;
  ExodusIIReader mloader;
  mloader.load("maxwell.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

  // Create an Hcurl space with default shapeset.
  HcurlSpace space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(biform<double, scalar>, biform<Ord, Ord>, HERMES_SYM);
  wf.add_vector_form(liform<double, scalar>, liform<Ord, Ord>);
   
  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est; 
  
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
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver, HERMES_HCURL_NORM);

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
    info("Calculating error estimate.");
    Adapt *adaptivity = new Adapt(&space, HERMES_HCURL_NORM);
    bool solutions_for_adapt = true;
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt) * 100;

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");

    // Clean up.
    delete ref_space->get_mesh();
    delete ref_space;
    delete matrix;
    delete rhs;
    delete solver;
    delete adaptivity;

    // Increase the counter of performed adaptivity steps.
    as++;

    done = true;

  } while (!done);

  return 0;
}
