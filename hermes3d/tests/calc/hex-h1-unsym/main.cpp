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

// Error should be smaller than this epsilon.
#define EPS								10e-10F

// Problem parameters.
int m = 2;
int n = 2;
int o = 2;

double fnc(double x, double y, double z)
{
	_F_
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
}

// Exact solution.
double exact_sln_fn(double x, double y, double z, double &dx, double &dy, double &dz)
{
	_F_
	dx = m * pow(x, m - 1) * pow(y, n) * pow(z, o) - 3 * pow(x, 2) * z + 2 * x * pow(y, 3);
	dy = n * pow(x, m) * pow(y, n - 1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o - 1) + 4 * pow(z, 3) - pow(x, 3);

	return fnc(x, y, z);
}

// Boundary condition types.
BCType bc_types(int marker)
{
	return BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	_F_
	return fnc(x, y, z);
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                    ExtData<Scalar> *data)
{
	_F_
	return
		int_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) +
		int_dudx_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename T>
T f(T x, T y, T z)
{
	_F_
	T ddxx = (m - 1) * m * pow(x, m - 2) * pow(y, n) * pow(z, o) - 6 * x * z + 2 * pow(y, 3);
	T ddyy = (n - 1) * n * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = (o - 1) * o * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);
	T dx = m * pow(x, m - 1) * pow(y, n) * pow(z, o) - 3 * pow(x, 2) * z + 2 * x * pow(y, 3);

	return -(ddxx + ddyy + ddzz) + dx;
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *data)
{
	_F_
	return int_F_v<Real, Scalar>(n, wt, f, u, e);
}

int main(int argc, char **args)
{
  // Test variable.
  int success_test = 1;

	if (argc < 2) error("Not enough parameters.");

  // Load the mesh.
	Mesh mesh;
  H3DReader mloader;
  if (!mloader.load(args[1], &mesh)) error("Loading mesh file '%s'.", args[1]);

	// Initialize the space.
	Ord3 o(4, 4, 4);
	H1Space space(&mesh, bc_types, essential_bc_values, o);

  // Initialize the weak formulation.
	WeakForm wf;
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_NONSYM);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>);

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

  ExactSolution ex_sln(&mesh, exact_sln_fn);
  
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
