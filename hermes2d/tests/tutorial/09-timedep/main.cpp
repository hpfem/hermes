#include "hermes2d.h"

// This test makes sure that example 09-timedep works correctly.
// CAUTION: This test will fail when any changes to the shapeset
// are made, but it is easy to fix (see below).

const int P_INIT = 3;            // initial polynomial degree in elements
const int INIT_REF_NUM = 0;      // number of initial uniform refinements
const double TAU = 200.0;        // time step in seconds
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem constants
const double T_INIT = 10;        // temperature of the ground (also initial temperature)
const double ALPHA = 10;         // heat flux coefficient for Newton's boundary condition
const double LAMBDA = 1e5;       // thermal conductivity of the material
const double HEATCAP = 1e6;      // heat capacity
const double RHO = 3000;         // material density
const double FINAL_TIME = 2100;  // length of time interval (24 hours) in seconds

// Global variable.
double TIME = 0;

// Time-dependent exterior temperature.
double temp_ext(double t) {
  return T_INIT + 10. * sin(2*M_PI*t/FINAL_TIME);
}

// Boundary markers.
int marker_ground = 1;
int marker_air = 2;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == marker_ground) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Function values for Dirichlet boundary markers.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return T_INIT;
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, u, v) / TAU +
         LAMBDA * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return HEATCAP * RHO * int_u_v<Real, Scalar>(n, wt, ext->fn[0], v) / TAU;
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return LAMBDA * ALPHA * int_u_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], 
Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return LAMBDA * ALPHA * temp_ext(TIME) * int_v<Real, Scalar>(n, wt, v);
}


int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("cathedral.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(2, 5);

  // Initialize an H1 space.
  H1Space* space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = get_num_dofs(space);
  printf("ndof = %d\n", ndof);

  // Set initial condition.
  Solution tsln;
  tsln.set_const(&mesh, T_INIT);

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, double>, bilinear_form<Ord, Ord>);
  wf.add_matrix_form_surf(bilinear_form_surf<double, double>, bilinear_form_surf<Ord, Ord>, marker_air);
  wf.add_vector_form(linear_form<double, double>, linear_form<Ord, Ord>, H2D_ANY, &tsln);
  wf.add_vector_form_surf(linear_form_surf<double, double>, linear_form_surf<Ord, Ord>, marker_air);

  // Initialize linear system.
  LinearProblem lp(&wf, space);

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, ndof, mat, rhs, solver);

  // time stepping
  int nsteps = (int)(FINAL_TIME/TAU + 0.5);
  bool rhsonly = false;
  for(int n = 1; n <= nsteps; n++)
  {
    info("---- Time step %d, time %3.5f, ext_temp %g", ts, TIME, temp_ext(TIME));

    // Assemble stiffness matrix and rhs.
    lp.assemble(mat, rhs, rhsonly);
    rhsonly = true;

    // Solve the matrix problem.
    if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");

    // Update tsln.
    tsln.set_fe_solution(space, rhs);

    // Update the time variable.
    TIME += TAU;
  }

  double sum = 0;
  for (int i=0; i < ndof; i++) sum += rhs->get(i);
  printf("coefficient sum = %g\n", sum);

  // Actual test. The value of 'sum' depend on the
  // current shapeset. If you change the shapeset,
  // you need to correct this number.
  int success = 1;
  if (fabs(sum - 4737.73) > 1e-1) success = 0;

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

