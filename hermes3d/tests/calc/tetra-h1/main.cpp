#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// The following parameters can be changed:
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// The error should be smaller than this epsilon.
#define EPS								10e-10F

// Problem parameters.
#define la0(x,y,z) (((y) + 1) / 2)
#define la1(x,y,z) (-(1 + (x) + (y) + (z)) / 2)
#define la2(x,y,z) (((x) + 1) / 2)
#define la3(x,y,z) (((z) + 1) / 2)
#define ph0 (-2.0 * 1.22474487139158904909864203735)

double fnc(double x, double y, double z) {
	return la0(x, y, z) * la1(x, y, z) * la2(x, y, z) * la3(x, y, z) * ph0 * ph0 * ph0;
}

// Exact solution.
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = 0.75 * sqrt(3.0 / 2.0) * (((x + 1) * (y + 1) * (z + 1)) - ((y + 1) * (- z - y - x - 1) * (z + 1)));
	dy = 0.75 * sqrt(3.0 / 2.0) * (((x + 1) * (y + 1) * (z + 1)) - ((x + 1) * (- z - y - x - 1) * (z + 1)));
	dz = 0.75 * sqrt(3.0 / 2.0) * (((x + 1) * (y + 1) * (z + 1)) - ((x + 1) * (y + 1) * (- z - y - x - 1)));

	return fnc(x, y, z);
}

// Boundary condition types.
BCType bc_types(int marker) 
{
	return BC_ESSENTIAL;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *user_data) {
	return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename T>
T dfnc(T x, T y, T z) {
	return
		- 1.5 * sqrt(3.0 / 2.0) * (y + 1) * (z + 1)
		- 1.5 * sqrt(3.0 / 2.0) * (x + 1) * (z + 1)
		- 1.5 * sqrt(3.0 / 2.0) * (x + 1) * (y + 1);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *user_data) {
	return int_F_v<Real, Scalar>(n, wt, dfnc, u, e);
}

int main(int argc, char **args) 
{
  // Test variable.
  int success_test = 1;

	if (argc < 3) error("Not enough parameters.");

  // Load the mesh.
	Mesh mesh;
  H3DReader mloader;
  if (!mloader.load(args[1], &mesh)) error("Loading mesh file '%s'.", args[1]);

  // Initialize the space according to the
  // command-line parameters passed.
	int o;
	sscanf(args[2], "%d", &o);
	H1Space space(&mesh, bc_types, NULL, o);

  // Initialize the weak formulation.
	WeakForm wf;
	wf.add_matrix_form(callback(bilinear_form), HERMES_SYM);
	wf.add_vector_form(callback(linear_form));

  // Initialize the FE problem.
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
  
  // Assemble the linear problem.
  info("Assembling (ndof: %d).", Space::get_num_dofs(&space));
  dp.assemble(matrix, rhs);
    
  // Solve the linear system. If successful, obtain the solution.
  info("Solving.");
		Solution sln(&mesh);
  if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
  else error ("Matrix solver failed.\n");

	ExactSolution ex_sln(&mesh, exact_solution);

  // Calculate exact error.
  info("Calculating exact error.");
  Adapt *adaptivity = new Adapt(&space, HERMES_H1_NORM);
  bool solutions_for_adapt = false;
  double err_exact = adaptivity->calc_err_exact(&sln, &ex_sln, solutions_for_adapt, HERMES_TOTAL_ERROR_ABS);

  if (err_exact > EPS)
		// Calculated solution is not precise enough.
		success_test = 0;

  // Clean up.
  delete matrix;
  delete rhs;
  delete solver;
  delete adaptivity;
  
  if (success_test) {
    info("Success!");
    return ERR_SUCCESS;
  }
	else {
    info("Failure!");
    return ERR_FAILURE;
	}
}

