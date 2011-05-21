#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// PDE: curl(curl(u)) = f.	
//
// BC:	u x n = 0 on Gamma_P
//  	curl(u) x n - i (n x u) x n = g on Gamma_I.

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
cplx img(0, 1);

// Exact solution.
scalar3 exact_solution(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz) {
	dx[0] = 3.*x*x*y*y - 3.*y*y*y*z;
	dy[0] = 2.*x*x*x*y - 9.*x*y*y*z;
	dz[0] = -3.*x*y*y*y;

	dx[1] = 3.*x*x*y*y*y*z*z*z + 4.*x*y*z;
	dy[1] = 3.*x*x*x*y*y*z*z*z + 2.*x*x*z;
	dz[1] = 3.*x*x*x*y*y*y*z*z + 2.*x*x*y;

	dx[2] = -12.*x*x;
	dy[2] = z*z*z;
	dz[2] = 3.*y*z*z;

	static scalar3 val(0.0, 0.0, 0.0);
	val[0] = x*x*x*y*y - 3.*x*y*y*y*z;
	val[1] = x*x*x*y*y*y*z*z*z + 2.*x*x*y*z;
	val[2] = y*z*z*z - 4.*x*x*x;
	return val;
}

template<typename S, typename T>
void exact_sln(S x, S y, S z, T (&fn)[3], T (&dx)[3], T (&dy)[3], T (&dz)[3]) {
	fn[0] = x*x*x*y*y - 3.*x*y*y*y*z;
	fn[1] = x*x*x*y*y*y*z*z*z + 2.*x*x*y*z;
	fn[2] = y*z*z*z - 4.*x*x*x;

	dx[0] = 3.*x*x*y*y - 3.*y*y*y*z;
	dx[1] = 3.*x*x*y*y*y*z*z*z + 4.*x*y*z;
	dx[2] = -12.*x*x;

	dy[0] = 2.*x*x*x*y - 9.*x*y*y*z;
	dy[1] = 3.*x*x*x*y*y*z*z*z + 2.*x*x*z;
	dy[2] = z*z*z;

	dz[0] = -3.*x*y*y*y;
	dz[1] = 3.*x*x*x*y*y*y*z*z + 2.*x*x*y;
	dz[2] = 3.*y*z*z;
}

template<typename Real, typename Scalar>
void f(Real x, Real y, Real z, Scalar (&val)[3]) {
	Scalar ev[3], dx[3], dy[3], dz[3];
	exact_sln(x, y, z, ev, dx, dy, dz);

	val[0] = 4*x*z + 18*x*y*z - 2*x*x*x + 9*x*x*y*y*z*z*z - ev[0];
	val[1] = -4*y*z + 3*z*z - 9*z*y*y + 6*y*x*x - 6*x*y*y*y*z*z*z - 6*z*x*x*x*y*y*y - ev[1];
	val[2] = 24*x + 2*x*x + 9*x*x*x*y*y*z*z - 3*y*y*y - ev[2];
}

/*
// TODO: this could be written in a much simpler way. Just use curl of exact solution
// and cross product defined in Scalar3D...
scalar3 &bc_values(int ess_bdy_marker, double x, double y, double z) {
	static scalar bc[3] = { 0., 0., 0. };

	switch (marker) {
		case 1:
			bc[1] = -4*x*y*z - 9*x*z*y*y + 2*y*x*x*x - 3*x*x*y*y*y*z*z*z;
			bc[2] = 12*x*x - 3*x*y*y*y;
			break;

		case 2:
			bc[1] = 4*x*y*z + 9*x*z*y*y - 2*y*x*x*x + 3*x*x*y*y*y*z*z*z;
			bc[2] = -12*x*x + 3*x*y*y*y;
			break;

		case 3:
			bc[0] = 4*x*y*z + 9*x*z*y*y - 2*y*x*x*x + 3*x*x*y*y*y*z*z*z;
			bc[2] = 2*y*x*x + 3*x*x*x*y*y*y*z*z - z*z*z;
			break;

		case 4:
			bc[0] = -4*x*y*z - 9*x*z*y*y + 2*y*x*x*x - 3*x*x*y*y*y*z*z*z;
			bc[2] = -2*y*x*x - 3*x*x*x*y*y*y*z*z + z*z*z;
			break;

		case 5:
			bc[0] = -12*x*x + 3*x*y*y*y;
			bc[1] = -2*y*x*x - 3*x*x*x*y*y*y*z*z + z*z*z;
			break;

		case 6:
			bc[0] = 12*x*x - 3*x*y*y*y;
			bc[1] = 2*y*x*x + 3*x*x*x*y*y*y*z*z - z*z*z;
			break;

		default:
			EXIT(H3D_ERR_FACE_INDEX_OUT_OF_RANGE);
	}

	switch (marker) {
		case 1:
		case 2:
			bc[1] -= img * (2*y*z*x*x + x*x*x*y*y*y*z*z*z);
			bc[2] -= img * (-4*x*x*x + y*z*z*z);
			break;

		case 3:
		case 4:
			bc[0] -= img * (x*x*x*y*y - 3*x*z*y*y*y);
			bc[2] -= img * (-4*x*x*x + y*z*z*z);
			break;

		case 5:
		case 6:
			bc[0] -= img * (x*x*x*y*y - 3*x*z*y*y*y);
			bc[1] -= img * (2*y*z*x*x + x*x*x*y*y*y*z*z*z);
			break;

		default:
			EXIT(H3D_ERR_FACE_INDEX_OUT_OF_RANGE);
	}

	return bc;
}
*/

// Boundary condition types.
BCType bc_types(int marker) 
{
	return H3D_BC_NATURAL;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) {
	return
		hcurl_int_curl_u_curl_v<Real, Scalar>(n, wt, u, v, e) -
		hcurl_int_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) {
	return hcurl_int_F_v<Real, Scalar>(n, wt, f, v, e);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) {
	return -img * hcurl_int_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) {
	Scalar result = 0.0;
	for (int i = 0; i < n; i++) {
		Scalar ev[3], dx[3], dy[3], dz[3];
		exact_sln(e->x[i], e->y[i], e->z[i], ev, dx, dy, dz);

		Scalar curl_e[3];
		calc_curl(dx, dy, dz, curl_e);
		Scalar tpe[3];
		calc_tan_proj(e->nx[i], e->ny[i], e->nz[i], ev, tpe);

		Scalar g[3] = {
			(e->nz[i] * curl_e[1] - e->ny[i] * curl_e[2]) - img * tpe[0],
			(e->nx[i] * curl_e[2] - e->nz[i] * curl_e[0]) - img * tpe[1],
			(e->ny[i] * curl_e[0] - e->nx[i] * curl_e[1]) - img * tpe[2],
		};
		result += wt[i] * (v->val0[i] * g[0] + v->val1[i] * g[1] + v->val2[i] * g[2]);
	}
	return result;
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
  wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_NONSYM);
  wf.add_matrix_form_surf(bilinear_form_surf<double, scalar>, bilinear_form_surf<Ord, Ord>);
  wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>);
  wf.add_vector_form_surf(linear_form_surf<double, scalar>, linear_form_surf<Ord, Ord>);

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
printf("err_exact = %lf", err_exact);
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

