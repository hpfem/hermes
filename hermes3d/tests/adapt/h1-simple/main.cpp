#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// Test to verify that hp-adaptivity works correctly.
//
// Starts with linear elements and should find the solution (just using p refinements).
//

// The following parameters can be changed:
int P_INIT_X, P_INIT_Y, P_INIT_Z;                 // Initial polynomial degree of all mesh elements taken from the command line arguments passed.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;                     // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).
const double TOLERANCE = 0.001;		          // error tolerance in percent.
const double THRESHOLD = 0.5;		          // error threshold for element refinement.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// The error should be smaller than this epsilon.
#define EPS								10e-10F

// Problem parameters.
int m = 2;
int n = 2;
int o = 2;



double fnc(double x, double y, double z) {
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
}

template<typename T>
T dfnc(T x, T y, T z) {
	T ddxx = m * (m - 1) * pow(x, m - 2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
	T ddyy = n * (n - 1) * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = o * (o - 1) * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);

	return -(ddxx + ddyy + ddzz);
}

// Exact solution.
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	// u(x, y, z) = pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4)
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 3) + 4 * pow(z, 3);

	return fnc(x, y, z);
}

// Boundary condition types.
BCType bc_types(int marker) 
{
	return BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z) {
	return fnc(x, y, z);
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) {
	return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *data) {
	return int_F_v<Real, Scalar>(n, wt, dfnc, u, e);
}

int main(int argc, char **args) 
{
  // Test variable.
  int success_test = 1;

	if (argc < 5) error("Not enough parameters.");

  // Load the mesh.
	Mesh mesh;
	H3DReader mloader;
	if (!mloader.load(args[1], &mesh)) error("Loading mesh file '%s'.", args[1]);
  
  // Initialize the space according to the
  // command-line parameters passed.
  sscanf(args[2], "%d", &P_INIT_X);
	sscanf(args[3], "%d", &P_INIT_Y);
	sscanf(args[4], "%d", &P_INIT_Z);
	Ord3 order(P_INIT_X, P_INIT_Y, P_INIT_Z);
  H1Space space(&mesh, bc_types, essential_bc_values, order);

  // Initialize the weak formulation.
	WeakForm wf;
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM, HERMES_ANY);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY);

	// Time measurement.
	TimePeriod cpu_time;
	cpu_time.tick();

	// Adaptivity loop.
  int as = 1; 
	bool done = false;
	do {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
  	Space* ref_space = construct_refined_space(&space, 1);
  
    out_orders_vtk(ref_space, "space", as);
	
	  // Initialize the FE problem.
	  bool is_linear = true;
	  DiscreteProblem lp(&wf, ref_space, is_linear);
		
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
    lp.assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();

    // Solve the linear system on reference mesh. If successful, obtain the solution.
    info("Solving on reference mesh.");
    Solution ref_sln(ref_space->get_mesh());
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), ref_space, &ref_sln);
    else {
		  info("Matrix solver failed.");
		  success_test = 0;
	  }
    
    // Time measurement.
    cpu_time.tick();

    // Project the reference solution on the coarse mesh.
    Solution sln(space.get_mesh());
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver);

    // Time measurement.
    cpu_time.tick();

	  // Calculate element errors and total error estimate.
    info("Calculating error estimate and exact error.");
    Adapt *adaptivity = new Adapt(&space, HERMES_H1_NORM);
    bool solutions_for_adapt = true;
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d.", Space::get_num_dofs(&space), Space::get_num_dofs(ref_space));
    info("err_est_rel: %g%%.", err_est_rel);

	  // If err_est_rel is too large, adapt the mesh. 
    if (err_est_rel < ERR_STOP) {
		  done = true;
      ExactSolution ex_sln(&mesh, exact_solution);
		  
      // Calculate exact error.
      info("Calculating exact error.");
      Adapt *adaptivity = new Adapt(&space, HERMES_H1_NORM);
      bool solutions_for_adapt = false;
      double err_exact = adaptivity->calc_err_exact(&sln, &ex_sln, solutions_for_adapt, HERMES_TOTAL_ERROR_ABS);

      if (err_exact > EPS)
		    // Calculated solution is not precise enough.
		    success_test = 0;
	 
      break;
	  }	
    else {
      info("Adapting coarse mesh.");
      adaptivity->adapt(THRESHOLD);
    }

    // If we reached the maximum allowed number of degrees of freedom, set the return flag to failure.
    if (Space::get_num_dofs(&space) >= NDOF_STOP)
    {
      done = true;
      success_test = 0;
    }

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
  
  if (success_test) {
    info("Success!");
    return ERR_SUCCESS;
  }
  else {
    info("Failure!");
    return ERR_FAILURE;
  }
}

