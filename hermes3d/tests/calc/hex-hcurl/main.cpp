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
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// The error should be smaller than this epsilon.
#define EPS								10e-10F

// Problem parameters.
const double alpha = 1.0;

// Exact solution.
scalar3 exact_solution(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz) {
	dx[0] = (1 - y*y)*(1 - z*z)*(z - 6*x*x);
	dx[1] = 2*(1 - x*x)*(1 - z*z) - 2*x*(1 - z*z)*(y*y*y + 2*x);
	dx[2] = -2*x*(1 - y*y)*(z*z - 3*x*y*z) - 3*y*z*(1 - x*x)*(1 - y*y);

	dy[0] = -2*y*(1 - z*z)*(1 - 2*x*x*x + x*z);
	dy[1] = 3*y*y*(1 - x*x)*(1 - z*z);
	dy[2] = -2*y*(1 - x*x)*(z*z - 3*x*y*z) - 3*x*z*(1 - x*x)*(1 - y*y);

	dz[0] = x*(1 - y*y)*(1 - z*z) - 2*z*(1 - y*y)*(1 - 2*x*x*x + x*z);
	dz[1] = -2*z*(1 - x*x)*(y*y*y + 2*x);
	dz[2] = (1 - x*x)*(1 - y*y)*(2*z - 3*x*y);

	static scalar3 val(0.0, 0.0, 0.0);
	val[0] = (1-y*y) * (1-z*z) * (x*z - 2*x*x*x + 1);
	val[1] = (1-x*x) * (1-z*z) * (y*y*y + 2*x);
	val[2] = (1-x*x) * (1-y*y) * (z*z - 3*x*y*z);

	return val;
}

template<typename S, typename T>
void exact_sln(S x, S y, S z, T (&val)[3], T (&dx)[3], T (&dy)[3], T (&dz)[3]) {
	val[0] = (1-y*y) * (1-z*z) * (x*z - 2*x*x*x + 1);
	val[1] = (1-x*x) * (1-z*z) * (y*y*y + 2*x);
	val[2] = (1-x*x) * (1-y*y) * (z*z - 3*x*y*z);

	dx[0] = (1 - y*y)*(1 - z*z)*(z - 6*x*x);
	dx[1] = 2*(1 - x*x)*(1 - z*z) - 2*x*(1 - z*z)*(y*y*y + 2*x);
	dx[2] = -2*x*(1 - y*y)*(z*z - 3*x*y*z) - 3*y*z*(1 - x*x)*(1 - y*y);

	dy[0] = -2*y*(1 - z*z)*(1 - 2*x*x*x + x*z);
	dy[1] = 3*y*y*(1 - x*x)*(1 - z*z);
	dy[2] = -2*y*(1 - x*x)*(z*z - 3*x*y*z) - 3*x*z*(1 - x*x)*(1 - y*y);

	dz[0] = x*(1 - y*y)*(1 - z*z) - 2*z*(1 - y*y)*(1 - 2*x*x*x + x*z);
	dz[1] = -2*z*(1 - x*x)*(y*y*y + 2*x);
	dz[2] = (1 - x*x)*(1 - y*y)*(2*z - 3*x*y);
}

template<typename S, typename T>
void f(S x, S y, S z, T (&val)[3]) {
	T ev[3], dx[3], dy[3], dz[3];
	exact_sln(x, y, z, ev, dx, dy, dz);

	T curlpart[3] = {
		2*(1 - y*y)*(1 - 2*x*x*x + x*z) + 2*(1 - z*z)*(1 - 2*x*x*x + x*z) - 6*x*y*y*(1 - z*z) - 3*y*(1 - x*x)*(1 - y*y) - 2*x*(1 - y*y)*(2*z - 3*x*y) + 4*x*z*(1 - y*y),
		2*(1 - x*x)*(y*y*y + 2*x) + 2*(1 - z*z)*(y*y*y + 2*x) + 8*x*(1 - z*z) - 3*x*(1 - x*x)*(1 - y*y) - 2*y*(1 - x*x)*(2*z - 3*x*y) - 2*y*(1 - z*z)*(z - 6*x*x),
		(1 - y*y)*(1 - z*z) + 2*(1 - x*x)*(z*z - 3*x*y*z) + 2*(1 - y*y)*(z*z - 3*x*y*z) - 6*z*y*y*(1 - x*x) - 2*z*(1 - y*y)*(z - 6*x*x) - 12*x*y*z*(1 - x*x) - 12*x*y*z*(1 - y*y)
	};

	val[0] = curlpart[0] - alpha * ev[0];
	val[1] = curlpart[1] - alpha * ev[1];
	val[2] = curlpart[2] - alpha * ev[2];
}

// Boundary condition types.
BCType bc_types(int marker) 
{
	return H3D_BC_ESSENTIAL;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) {
	return
		hcurl_int_curl_u_curl_v<Real, Scalar>(n, wt, u, v, e) -
		alpha * hcurl_int_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *data) {
	return hcurl_int_F_v<Real, Scalar>(n, wt, f<Real, Scalar>, u, e);
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
	Ord3 order(o, o, o);
	HcurlSpace space(&mesh, bc_types, NULL, order);

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
  Adapt *adaptivity = new Adapt(&space, HERMES_HCURL_NORM);
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

