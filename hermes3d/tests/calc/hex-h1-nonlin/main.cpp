#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

//  Testing nonlinear solver using Newton's method.

// The following parameters influence the projection system solving:
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


#if defined LINEAR

//  Solving a linear problem using nonlinear solver with Newton's method.
//
//  PDE: stationary heat transfer with nonlinear thermal conductivity.
//  -\Laplace u = f.
//
//  Exact solution: u(x,y,z) = x^2 + y^2 + z^2.

double fnc(double x, double y, double z)
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

// Boundary condition types.
BCType bc_types(int marker)
{
	return BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return fnc(x, y, z);
}

template<typename Real, typename Scalar>
Scalar jacobi_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e,
                  ExtData<Scalar> *data)
{
	return int_grad_u_grad_v<Real, Scalar>(n, wt, vi, vj, e);
}

template<typename Real, typename Scalar>
Scalar resid_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *vi, Geom<Real> *e,
                 ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (grad_grad(u_ext[0], vi) + 6.0 * vi->val[i]);
	return res;
}

#elif defined NONLIN1

//  Solving a nonlinear problem using Newton's method.
//
//  PDE: stationary heat transfer with nonlinear thermal conductivity.
//  - div[lambda(u) grad u] = f,
//                lambda(u) = 1 + u^2.
//
//  BC:  T = 100 on the left, top and bottom edges,
//       dT/dn = 0 on the face.
//
//  Exact solution: u(x,y,z) = 100.

// Thermal conductivity (temperature-dependent).
// For any u, this function has to be  positive in the entire domain!
template<typename T>
inline T lambda(T temp)  { return 10 + 0.1 * pow(temp, 2); }

// Derivate of lambda wrt temperature.
template<typename T>
inline T dlambda(T temp) { return 0.2 * temp; }

const int marker_right = 2;

// Boundary condition types.
BCType bc_types(int marker)
{
	if (marker == marker_right) return BC_NATURAL;
	else return BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return 100.0;
}

template<typename Real, typename Scalar>
Scalar jacobi_form(int n, double *wt, Func<Scalar> *itr[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                  ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] *
			(dlambda(itr[0]->val[i]) * u->val[i] * grad_grad(itr[0], v) +
			  lambda(itr[0]->val[i]) * grad_grad(u, v));
	return res;
}

template<typename Real, typename Scalar>
Scalar resid_form(int n, double *wt, Func<Scalar> *itr[], Func<Real> *v, Geom<Real> *e,
                 ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (lambda(itr[0]->val[i]) * grad_grad(itr[0], v));
	return res;
}

#elif defined NONLIN2

//  Solving a nonlinear problem using Newton's method.
//
//  PDE: stationary heat transfer with nonlinear thermal conductivity.
//  - div [u grad u] = f.
//
//  BC:  u = x^2 + y^2 + z^2 on dOmega.
//
//  Exact solution: u(x,y,z) = x^2 + y^2 + z^2.

double fnc(double x, double y, double z)
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

// Boundary condition types.
BCType bc_types(int marker)
{
	return BC_ESSENTIAL;
}

// Dirichlet boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return fnc(x, y, z);
}

template<typename T>
inline T f(T x, T y, T z)
{
	return -10.0 * (x*x + y*y + z*z);
}


template<typename Real, typename Scalar>
Scalar jacobi_form(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Func<Real> *vj, Geom<Real> *e,
                  ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (vj->val[i] * grad_grad(u[0], vi) + u[0]->val[i] * grad_grad(vj, vi));
	return res;
}


template<typename Real, typename Scalar>
Scalar resid_form(int n, double *wt, Func<Real> *u[], Func<Real> *vi, Geom<Real> *e,
                 ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (u[0]->val[i] * grad_grad(u[0], vi) - f(e->x[i], e->y[i], e->z[i]) * vi->val[i]);
	return res;
}

template<typename Real, typename Scalar>
Scalar biproj_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e,
                  ExtData<Scalar> *data)
{
	Scalar res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (u->val[i] * v->val[i] + grad_grad(u, v));
	return res;
}

template<typename Real, typename Scalar>
Scalar liproj_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data)
{
	Scalar res = 0.0;
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
#if defined NONLIN1
	Ord3 order(1, 1, 1);
#else
	Ord3 order(2, 2, 2);
#endif
	H1Space space(&mesh, bc_types, essential_bc_values, order);

#if defined NONLIN2
	// Do L2 projection of zero function.
	WeakForm proj_wf;
	proj_wf.add_matrix_form(biproj_form<double, scalar>, biproj_form<Ord, Ord>, HERMES_SYM);
	proj_wf.add_vector_form(liproj_form<double, scalar>, liproj_form<Ord, Ord>);

	bool is_linear = true;
	DiscreteProblem lp(&proj_wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver_proj = create_linear_solver(matrix_solver, matrix, rhs);
  
  // Initialize the preconditioner in the case of SOLVER_AZTECOO.
  if (matrix_solver == SOLVER_AZTECOO) 
  {
    ((AztecOOSolver*) solver_proj)->set_solver(iterative_method);
    ((AztecOOSolver*) solver_proj)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }
  
  // Assemble the linear problem.
  info("Assembling (ndof: %d).", Space::get_num_dofs(&space));
  lp.assemble(matrix, rhs);
    
  // Solve the linear system.
  info("Solving.");
  if(!solver_proj->solve()) error ("Matrix solver failed.\n");

  delete matrix;
  delete rhs;
#endif

	// Initialize the weak formulation.
	WeakForm wf(1);
	wf.add_matrix_form(0, 0, jacobi_form<double, scalar>, jacobi_form<Ord, Ord>, HERMES_NONSYM);
	wf.add_vector_form(0, resid_form<double, scalar>, resid_form<Ord, Ord>);

	// Initialize the FE problem.
#if defined NONLIN2
  is_linear = false;
#else
  bool is_linear = false;
#endif
	DiscreteProblem dp(&wf, &space, is_linear);

	NoxSolver solver(&dp);

#if defined NONLIN2
solver.set_init_sln(solver_proj->get_solution());
delete solver_proj;
#endif

solver.set_conv_iters(10);
	info("Solving.");
	Solution sln(&mesh);
	if(solver.solve()) Solution::vector_to_solution(solver.get_solution(), &space, &sln);
  else error ("Matrix solver failed.\n");

		Solution ex_sln(&mesh);
#ifdef NONLIN1
		ex_sln.set_const(100.0);
#else
		ex_sln.set_exact(exact_solution);
#endif
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

