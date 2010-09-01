#include "hermes2d.h"

// This test makes sure that example 07-general works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

const int P_INIT = 2;             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;       // Number of initial uniform refinements
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters
double a_11(double x, double y) {
  if (y > 0) return 1 + x*x + y*y;
  else return 1;
}

double a_22(double x, double y) {
  if (y > 0) return 1;
  else return 1 + x*x + y*y;
}

double a_12(double x, double y) {
  return 1;
}

double a_21(double x, double y) {
  return 1;
}

double a_1(double x, double y) {
  return 0.0;
}

double a_2(double x, double y) {
  return 0.0;
}

double a_0(double x, double y) {
  return 0.0;
}

double rhs(double x, double y) {
  return 1 + x*x + y*y;
}

double g_D(double x, double y) {
  return -cos(M_PI*x);
}

double g_N(double x, double y) {
  return 0;
}

/********** Boundary conditions ***********/

// Boundary condition types
BCType bc_types(int marker)
{
  if (marker == 1) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Dirichlet boundary condition values
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return g_D(x, y);
}

// (Volumetric) bilinear form
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    double x = e->x[i];
    double y = e->y[i];
    result += (a_11(x, y)*u->dx[i]*v->dx[i] +
               a_12(x, y)*u->dy[i]*v->dx[i] +
               a_21(x, y)*u->dx[i]*v->dy[i] +
               a_22(x, y)*u->dy[i]*v->dy[i] +
               a_1(x, y)*u->dx[i]*v->val[i] +
               a_2(x, y)*u->dy[i]*v->val[i] +
               a_0(x, y)*u->val[i]*v->val[i]) * wt[i];
  }
  return result;
}

// Integration order for the bilinear form
template<typename Real, typename Scalar>
Scalar bilinear_form_ord(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                         Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return u->val[0] * v->val[0] * e->x[0] * e->x[0]; // returning the sum of the degrees of the basis
                                                    // and test function plus two
}

// Surface linear form (natural boundary conditions)
template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, g_N, v, e);
}

// Integration order for surface linear form
template<typename Real, typename Scalar>
Scalar linear_form_surf_ord(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return v->val[0] * e->x[0] * e->x[0];  // returning the polynomial degree of the test function plus two
}

// Volumetric linear form (right-hand side)
template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}

// Integration order for the volumetric linear form
template<typename Real, typename Scalar>
Scalar linear_form_ord(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return v->val[0] * e->x[0] * e->x[0];  // returning the polynomial degree of the test function plus two
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);
  for (int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space.
  H1Space* space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form, bilinear_form_ord, H2D_SYM);
  wf.add_vector_form(linear_form, linear_form_ord);
  wf.add_vector_form_surf(linear_form_surf, linear_form_surf_ord, 2);

  // Testing n_dof and correctness of solution vector
  // for p_init = 1, 2, ..., 10
  int success = 1;
  Solution sln;
  for (int p_init = 1; p_init <= 10; p_init++) {

    printf("********* p_init = %d *********\n", p_init);
    space->set_uniform_order(p_init);

    // Initialize the linear problem.
    LinearProblem lp(&wf, space);

    // Select matrix solver.
    Matrix* mat; Vector* rhs; CommonSolver* solver;
    init_matrix_solver(matrix_solver, get_num_dofs(space), mat, rhs, solver);

    // Assemble stiffness matrix and rhs.
    bool rhsonly = false;
    lp.assemble(mat, rhs, rhsonly);

    // Solve the matrix problem.
    if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

    int ndof = get_num_dofs(space);
    printf("ndof = %d\n", ndof);
    double sum = 0;
    for (int i=0; i < ndof; i++) sum += rhs->get(i);
    printf("coefficient sum = %g\n", sum);

    // Actual test. The values of 'sum' depend on the
    // current shapeset. If you change the shapeset,
    // you need to correct these numbers.
    if (p_init == 1 && fabs(sum - 1.72173) > 1e-2) success = 0;
    if (p_init == 2 && fabs(sum - 0.639908) > 1e-2) success = 0;
    if (p_init == 3 && fabs(sum - 0.826367) > 1e-2) success = 0;
    if (p_init == 4 && fabs(sum - 0.629395) > 1e-2) success = 0;
    if (p_init == 5 && fabs(sum - 0.574235) > 1e-2) success = 0;
    if (p_init == 6 && fabs(sum - 0.62792) > 1e-2) success = 0;
    if (p_init == 7 && fabs(sum - 0.701982) > 1e-2) success = 0;
    if (p_init == 8 && fabs(sum - 0.7982) > 1e-2) success = 0;
    if (p_init == 9 && fabs(sum - 0.895069) > 1e-2) success = 0;
    if (p_init == 10 && fabs(sum - 1.03031) > 1e-2) success = 0;
  }

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  if (success == 1) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
