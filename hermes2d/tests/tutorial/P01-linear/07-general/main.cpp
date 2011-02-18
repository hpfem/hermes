#include "hermes2d.h"

// This test makes sure that example 07-general works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

const int P_INIT = 2;             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

const char* iterative_method = "cg";              // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).

// Boundary markers.
const std::string BDY_HORIZONTAL = "Boundary horizontal";
const std::string BDY_VERTICAL = "Boundary vertical";

// Problem parameters.
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

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(double x, double y)
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
Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                         Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(u->val[0] * v->val[0] * e->x[0] * e->x[0]); // returning the sum of the degrees of the basis
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
Ord linear_form_surf_ord(int n, double *wt, Func<Ord> *u_ext[], 
Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(v->val[0] * e->x[0] * e->x[0]);  // returning the polynomial degree of the test function plus two
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

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_HORIZONTAL);
  bc_types.add_bc_neumann(BDY_VERTICAL);

  BCValues bc_values;
  bc_values.add_function(BDY_HORIZONTAL, essential_bc_values);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form, bilinear_form_ord, HERMES_SYM);
  wf.add_vector_form(linear_form, linear_form_ord);
  wf.add_vector_form_surf(linear_form_surf, linear_form_surf_ord, BDY_VERTICAL);

  // Testing n_dof and correctness of solution vector
  // for p_init = 1, 2, ..., 10
  int success = 1;
  Solution sln;
  for (int p_init = 1; p_init <= 10; p_init++) {

    printf("********* p_init = %d *********\n", p_init);
    space.set_uniform_order(p_init);

    // Initialize the FE problem.
    bool is_linear = true;
    DiscreteProblem dp(&wf, &space, is_linear);

    // Set up the solver, matrix, and rhs according to the solver selection. 
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
    
    // Initialize the solution.
    Solution sln;

    // Assemble the stiffness matrix and right-hand side vector.  
    info("Assembling the stiffness matrix and right-hand side vector.");
    bool rhsonly = false;
    dp.assemble(matrix, rhs, rhsonly);

    // Solve the linear system and if successful, obtain the solution.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution::vector_to_solution(solver->get_solution(), &space, &sln);
    else
      error ("Matrix solver failed.\n");
  
    int ndof = Space::get_num_dofs(&space);
    printf("ndof = %d\n", ndof);
    double sum = 0;
    for (int i=0; i < ndof; i++) sum += solver->get_solution()[i];
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

  if (success == 1) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
