#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// This tets makes sure that tetrahedra work with Neumann BC.
//
// PDE: -Laplace u + u = f
//
// BC: Neumann, given by exact solution.

// The following parameters can be changed:
int P_INIT;                                       // Polynomial degree of the element, to be read from command-line.
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

// Exact solution.
template<typename T>
T fnc(T x, T y, T z) 
{
  return x*x + y*y + z*z;
}

// Right-hand side.
template<typename T>
T dfnc(T x, T y, T z) 
{
  return -6.0 + fnc(x, y, z);
}

// Exact solution with derivatives.
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
  return H3D_BC_NATURAL;
}

// Weak forms.
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *user_data) 
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e) + int_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *user_data) 
{
  return int_F_v<Real, Scalar>(n, wt, dfnc, u, e);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int np, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *user_data) 
{
  //printf("marker: %d\n", e->marker); // debug
  Scalar result = 0;
  Scalar dx, dy, dz;
  for (int i = 0; i < np; i++) {
    dx = 2. * e->x[i];
    dy = 2. * e->y[i];
    dz = 2. * e->z[i];
    //printf("normal vector: %g %g %g\n", e->nx[i], e->ny[i], e->nz[i]); // debug - print normal vector to face
    result += wt[i] * (u->val[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i]));
    //result += wt[i];  // debug - this is to check whether weights on a face add up to the area of that face.
  }
  //printf("result: %g\n", result); // debug
  return result;
}

Ord linear_form_surf_ord(int np, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Geom<Ord> *e, ExtData<Ord> *user_data) 
{
  return Ord(2*P_INIT);
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
  sscanf(args[2], "%d", &P_INIT);
  H1Space space(&mesh, bc_types, NULL, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
  wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>);
  wf.add_vector_form_surf(linear_form_surf<double, scalar>, linear_form_surf_ord);

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

