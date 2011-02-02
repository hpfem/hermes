#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This example solves a simple linear wave equation by converting it 
// into a system of two first-order equations in time. Time discretization
// is done using the implicit Euler method.
//
// PDE: \frac{\partial^2 u}{\partial t^2} - \Delta u = 0,
// converted into
//
//      \frac{\partial u}{\partial t} - v = 0,
//      \frac{\partial v}{\partial t} - \Delta u = 0.
//
// BC:  u = 0 on the boundary,
//      v = 0 on the boundary (u = 0 => \partial u / \partial t = 0).
//
// IC:  smooth peak for u, zero for v.
//
// The following parameters can be changed:

const int P_INIT = 6;                              // Initial polynomial degree of all elements.
const double TAU = 0.01;                           // Time step.
const double T_FINAL = 2.15;                       // Final time.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                   // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const int BDY = 1;

// Problem parameters.
const double C_SQUARED = 100;                      // Square of wave speed.                     

// Initial condition for u.
double init_cond_u(double x, double y, double& dx, double& dy) {
  dx = exp(-x*x - y*y) * (-2*x);
  dy = exp(-x*x - y*y) * (-2*y);
  return exp(-x*x - y*y);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Convert to quadrilaterals.
  mesh.convert_triangles_to_quads();

  // Refine once towards vertex #4.
  mesh.refine_towards_vertex(4, 1);

  // Refine towards boundary.
  mesh.refine_towards_boundary(1, 1);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY);

  // Enter Dirichlet boundary values;
  BCValues bc_values;
  bc_values.add_zero(BDY);

  // Create x- and y- displacement space using the default H1 shapeset.
  H1Space u_space(&mesh, &bc_types, &bc_values, P_INIT);
  H1Space v_space(&mesh, &bc_types, &bc_values, P_INIT);
  info("ndof = %d.", Space::get_num_dofs(Hermes::vector<Space *>(&u_space, &v_space)));

  // Initialize solutions.
  Solution u_sln(&mesh, init_cond_u);
  Solution v_sln(&mesh, 0.0);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0));  
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1));  
  wf.add_matrix_form(1, 0, callback(bilinear_form_1_0));  
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1));  
  wf.add_vector_form(0, callback(linear_form_0), HERMES_ANY, &u_sln);  
  wf.add_vector_form(1, callback(linear_form_1), HERMES_ANY, &v_sln);  

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, Hermes::vector<Space *>(&u_space, &v_space), is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);

  // Initialize views.
  ScalarView u_view("Solution u", new WinGeom(0, 0, 500, 400));
  //u_view.show_mesh(false);
  u_view.fix_scale_width(50);
  ScalarView v_view("Solution v", new WinGeom(510, 0, 500, 400));
  //v_view.show_mesh(false);
  v_view.fix_scale_width(50);

  // Time stepping loop.
  double current_time = TAU; int ts = 1;
  bool rhs_only = false;
  do 
  {
    info("---- Time step %d, t = %g s.", ts, current_time); ts++;

    // First time assemble both the stiffness matrix and right-hand side vector,
    // then just the right-hand side vector.
    if (rhs_only == false) info("Assembling the stiffness matrix and right-hand side vector.");
    else info("Assembling the right-hand side vector (only).");
    dp.assemble(matrix, rhs, rhs_only);
    rhs_only = true;

    // Solve the linear system and if successful, obtain the solutions.
    info("Solving the matrix problem.");
    if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), 
                                                      Hermes::vector<Space *>(&u_space, &v_space), 
                                                      Hermes::vector<Solution *>(&u_sln, &v_sln));
    else error ("Matrix solver failed.\n");
  
    // Visualize the solutions.
    char title[100];
    sprintf(title, "Solution u, t = %g", current_time);
    u_view.set_title(title);
    u_view.show(&u_sln);
    sprintf(title, "Solution v, t = %g", current_time);
    v_view.set_title(title);
    v_view.show(&v_sln);

    // Update time.
    current_time += TAU;

  } while (current_time < T_FINAL);

  // Wait for the view to be closed.
  View::wait();

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  return 0;
}

