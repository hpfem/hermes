#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
#include <math.h>
#include <hermes3d.h>

const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).

#if defined GMSH
	GmshOutputEngine output(stdout);
#elif defined VTK
	VtkOutputEngine output(stdout, 1);
#endif


int m, n, o;

// error should be smaller than this epsilon
#define EPSILON								10e-10F

double fnc(double x, double y, double z)
{
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
}

template<typename T>
T dfnc(T x, T y, T z)
{
	T ddxx = m * (m - 1) * pow(x, m - 2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
	T ddyy = n * (n - 1) * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = o * (o - 1) * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);

	return -(ddxx + ddyy + ddzz);
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz)
{
	// u(x, y, z) = x^m * y^n * z^o + x^2 * y^3 - x^3 * z + z^4
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 3) + 4 * pow(z, 3);

	return fnc(x, y, z);
}

// exact solution for vector functions
double exact_solution0(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = dy = dz = 0;
	return 1.0;
}

double exact_solution1(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = dy = dz = 0;
	return 0.0;
}

double exact_solution2(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = dy = dz = 0;
	return 0.0;
}

//
scalar3 exact_vec_solution(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz)
{
  static scalar3 val(0.0, 0.0, 0.0);

	dx[0] = dx[1] = dx[2] = 0;
	dy[0] = dy[1] = dy[2] = 0;
	dz[0] = dz[1] = dz[2] = 0;

	val[0] = 0.0;
	val[1] = 1.0;
	val[2] = 0.0;

	return val;
}

// Boundary condition types.
BCType bc_types(int marker)
{
	return H3D_BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return fnc(x, y, z);
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int np, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                    ExtData<Scalar> *data)
{
	return int_grad_u_grad_v<Real, Scalar>(np, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *data)
{
	return int_F_v<Real, Scalar>(n, wt, dfnc, u, e);
}

#if defined WITH_UMFPACK
#define StiffMatrix		UMFPackMatrix
#elif defined WITH_PETSC
#define StiffMatrix		PetscMatrix
#elif defined WITH_MUMPS
#define StiffMatrix		MumpsMatrix
#elif defined WITH_TRILINOS
#define StiffMatrix		EpetraMatrix
#endif

void test_mat(Mesh *mesh, StiffMatrix &mat)
{
#if defined WITH_UMFPACK
MatrixSolverType matrix_solver = SOLVER_UMFPACK;
#elif defined WITH_PETSC
MatrixSolverType matrix_solver = SOLVER_PETSC;
#elif defined WITH_MUMPS
MatrixSolverType matrix_solver = SOLVER_MUMPS;
#elif defined WITH_TRILINOS 
MatrixSolverType matrix_solver = SOLVER_AZTECOO;
#endif

	m = n = o = 2;
	int mx = maxn(4, m, n, o, 4);
	Ord3 order(mx, mx, mx);

	// Create an H1 space with default shapeset.
	H1Space space(mesh, bc_types, essential_bc_values, order);

	WeakForm wf(1);
	wf.add_matrix_form(0, 0, bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
	wf.add_vector_form(0, linear_form<double, scalar>, linear_form<Ord, Ord>);

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

	// Assemble stiffness matrix and load vector.
	info("Assembling the linear problem (ndof: %d).", Space::get_num_dofs(&space));
	dp.assemble(matrix, rhs);

	// Solve the linear system. If successful, obtain the solution.
	info("Solving the linear problem.");
	Solution sln(space.get_mesh());
	if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
	else error ("Matrix solver failed.\n");

}

void test_mm(Mesh *mesh)
{
	_F_
	// testing a visualization of a vector-valued solution, where each component is on
	// a different mesh
	Mesh mesh0, mesh1, mesh2;
	mesh0.copy(*mesh);
	mesh1.copy(*mesh);
	mesh2.copy(*mesh);

	mesh0.refine_element(1, H3D_REFT_HEX_X);
	mesh1.refine_element(1, H3D_REFT_HEX_Y);
	mesh2.refine_element(1, H3D_REFT_HEX_Z);

	ExactSolution ex_sln0(&mesh0, exact_solution0);
	ExactSolution ex_sln1(&mesh1, exact_solution1);
	ExactSolution ex_sln2(&mesh2, exact_solution2);
	output.out(&ex_sln0, &ex_sln1, &ex_sln2, "U");
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args)
{
	if (argc < 3) error("Not enough parameters.");

	char *type = args[1];

	Mesh mesh;
	H3DReader mesh_loader;
	if (!mesh_loader.load(args[2], &mesh)) error("Loading mesh file '%s'\n", args[2]);

	if (strcmp(type, "sln") == 0) {
		// Testing on Exact solution which always gives the same value (values from Solution may differ by epsilon)
		ExactSolution ex_sln(&mesh, exact_solution);
		output.out(&ex_sln, "U");
	}
	else if (strcmp(type, "vec-sln") == 0) {
		// Testing on Exact solution which always gives the same value (values from Solution may differ by epsilon)
		ExactSolution ex_sln(&mesh, exact_vec_solution);
		output.out(&ex_sln, "U");
	}
	else if (strcmp(type, "3sln") == 0) {
		// Testing on Exact solution which always gives the same value (values from Solution may differ by epsilon)
		ExactSolution ex_sln0(&mesh, exact_solution0);
		ExactSolution ex_sln1(&mesh, exact_solution1);
		ExactSolution ex_sln2(&mesh, exact_solution2);
		output.out(&ex_sln0, &ex_sln1, &ex_sln2, "U");
	}
	else if (strcmp(type, "ord") == 0) {

		Ord3 order;
		if (mesh.elements[1]->get_mode() == HERMES_MODE_HEX)
			order = Ord3(2, 3, 4);
		else if (mesh.elements[1]->get_mode() == HERMES_MODE_TET)
			order = Ord3(3);
		else
			error(HERMES_ERR_NOT_IMPLEMENTED);

		H1Space space(&mesh, bc_types, essential_bc_values, order);

#if defined GMSH
		output.out_orders_gmsh(&space, "orders_gmsh");
#elif defined VTK
		output.out_orders_vtk(&space, "orders_vtk");
#endif
	}
	else if (strcmp(type, "bc") == 0) {

#if defined GMSH
		output.out_bc_gmsh(&mesh);
#elif defined VTK
		output.out_bc_vtk(&mesh);
#endif
	}
	else if (strcmp(type, "mat") == 0) {
		StiffMatrix mat;
		test_mat(&mesh, mat);
		output.out(&mat);
	}
	else if (strcmp(type, "mm") == 0) {
		test_mm(&mesh);
	}

	return 0;
}
