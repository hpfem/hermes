#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// Solving a linear problem using JFNK method
//
// Dirichlet BC:
//
//   -\Delta u = f in \Omega
//           u = g on d\Omega
//
// Neumann BC:
//
//   -\Delta u + u = f in \Omega
//           du/dn = g on d\Omega
//
// Newton BC:
//
//   -\Delta u = f in \Omega
//   du/dn + u = g on d\Omega
//
//   u(x,y,z) = x^2 + y^2 + z^2
//

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

#define grad_grad(u, v) (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->dz[i] * v->dz[i])

// Problem parameters.
template<typename T>
T fnc(T x, T y, T z)
{
	return x*x + y*y + z*z;
}

// Exact solution.
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return fnc(x, y, z);
}

#ifdef LIN_DIRICHLET

// case with dirichlet BC

template<typename T>
T dfnc(T x, T y, T z)
{
	return -6.0;
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
Scalar form_0(int n, double *wt, Func<Scalar> *u[], Func<Real> *vi, Geom<Real> *e,
             ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(u[0], vi) - dfnc(e->x[i], e->y[i], e->z[i]) * vi->val[i]);
	return res;
}

// Preconditioning.
template<typename Real, typename Scalar>
Scalar precond_0_0(int n, double *wt, Func<Real> *u[0], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e,
                  ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(vi, vj));
	return res;
}

#elif defined LIN_NEUMANN

// Case with neumann BC.

template<typename T>
T dfnc(T x, T y, T z)
{
	T ddxx = 2;
	T ddyy = 2;
	T ddzz = 2;

	return -(ddxx + ddyy + ddzz) + fnc(x, y, z);
}

// Boundary condition types.
BCType bc_types(int marker)
{
	return H3D_BC_NATURAL;
}

template<typename Real, typename Scalar>
Scalar form_0(int n, double *wt, Func<Real> *u[0], Func<Real> *vi, Geom<Real> *e,
             ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(u[0], vi)
				+ u[0]->val[i] * vi->val[i]
				- dfnc(e->x[i], e->y[i], e->z[i]) * vi->val[i]);
	return res;
}

template<typename Real, typename Scalar>
Scalar form_0_surf(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e,
                  ExtData<Scalar> *data)
{
	Scalar result = 0;
	for (int i = 0; i < n; i++) {
		Scalar dx = 2 * e->x[i];
		Scalar dy = 2 * e->y[i];
		Scalar dz = 2 * e->z[i];

		result += wt[i] * (vi->val[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i]));
	}
	return -result;
}

// precond
template<typename Real, typename Scalar>
Scalar precond_0_0(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e,
                  ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(vi, vj));
	return res;
}

#elif defined LIN_NEWTON

// Case with newton BC.

template<typename T>
T dfnc(T x, T y, T z)
{
	T ddxx = 2;
	T ddyy = 2;
	T ddzz = 2;

	return -(ddxx + ddyy + ddzz);
}

// Boundary condition types.
BCType bc_types(int marker)
{
	return H3D_BC_NATURAL;
}

template<typename Real, typename Scalar>
Scalar form_0(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e,
             ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(u[0], vi)
				- dfnc(e->x[i], e->y[i], e->z[i]) * vi->val[i]);
	return res;
}

template<typename Real, typename Scalar>
Scalar form_0_surf(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e,
                  ExtData<Scalar> *data)
{
	Scalar res = 0;
	for (int i = 0; i < n; i++) {
		Scalar dx = 2 * e->x[i];
		Scalar dy = 2 * e->y[i];
		Scalar dz = 2 * e->z[i];

		res += wt[i] * (u[0]->val[i] * vi->val[i]
				- (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc(e->x[i], e->y[i], e->z[i])) * vi->val[i]);
	}
	return res;
}

#elif defined NLN_DIRICHLET

// Nonlinear case with dirichlet BC.

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

template<typename T>
inline T f(T x, T y, T z)
{
	return -10 * (x*x + y*y + z*z);
}


template<typename Real, typename Scalar>
Scalar form_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *vi, Geom<Real> *e,
             ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (u_ext[0]->val[i] * grad_grad(u_ext[0], vi) - f(e->x[i], e->y[i], e->z[i]) * vi->val[i]);
	return res;
}

// precond
template<typename Real, typename Scalar>
Scalar precond_0_0(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e,
                  ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(vi, vj));
	return res;
}

#endif



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
	int mx = 2;
	Ord3 order(mx, mx, mx);
	H1Space space(&mesh, bc_types, NULL, order);
#if defined LIN_DIRICHLET || defined NLN_DIRICHLET
	space.set_essential_bc_values(essential_bc_values);
#endif
	// Initialize the weak formulation.
	WeakForm wf;
	wf.add_vector_form(form_0<double, scalar>, form_0<Ord, Ord>);
#if defined LIN_NEUMANN || defined LIN_NEWTON
	wf.add_vector_form_surf(form_0_surf<double, scalar>, form_0_surf<Ord, Ord>);
#endif
#if defined LIN_DIRICHLET || defined NLN_DIRICHLET
	// preconditioner
	wf.add_matrix_form(precond_0_0<double, scalar>, precond_0_0<Ord, Ord>, HERMES_SYM);
#endif

	// Initialize the FE problem.
	DiscreteProblem fep(&wf, &space);

#if defined LIN_DIRICHLET || defined NLN_DIRICHLET
	// use ML preconditioner to speed-up things
	MlPrecond pc("sa");
	pc.set_param("max levels", 6);
	pc.set_param("increasing or decreasing", "decreasing");
	pc.set_param("aggregation: type", "MIS");
	pc.set_param("coarse: type", "Amesos-KLU");
#endif

	NoxSolver solver(&fep);
#if defined LIN_DIRICHLET || defined NLN_DIRICHLET
//	solver.set_precond(&pc);
#endif

	info("Solving.");
	Solution sln(&mesh);
	if(solver.solve()) Solution::vector_to_solution(solver.get_solution(), &space, &sln);
  else error ("Matrix solver failed.\n");
	

		Solution ex_sln(&mesh);
		ex_sln.set_exact(exact_solution);

		// Calculate exact error.
  info("Calculating exact error.");
  Adapt *adaptivity = new Adapt(&space, HERMES_H1_NORM);
  bool solutions_for_adapt = false;
  double err_exact = adaptivity->calc_err_exact(&sln, &ex_sln, solutions_for_adapt, HERMES_TOTAL_ERROR_ABS);

  if (err_exact > EPS)
		// Calculated solution is not precise enough.
		success_test = 0;

  if (success_test) {
    info("Success!");
    return ERR_SUCCESS;
  }
  else {
    info("Failure!");
    return ERR_FAILURE;
  }
}

