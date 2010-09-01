#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This is a long version of example "waveguide": function solve_linear_adapt() is not used.
//
// PDE: time-harmonic Maxwell's equations;
//      there is circular load in the middle of the large cavity, whose permittivity
//      is different from the rest of the domain.
//
// Domain: square cavity with another small square cavity attached from outside
//         on the right.
//
// Meshes: you can either use "oven_load_circle.mesh" containing curved elements
//         aligned with the circular load, or "oven_load_square.mesh" which is not
//         aligned.
//
// BC: perfect conductor on the boundary except for the right-most edge of the small
//     cavity, where a harmonic surface current is prescribed
//
// The following parameters can be changed:

const int INIT_REF_NUM = 0;              // Number of initial uniform mesh refinements.
const int P_INIT = 2;                    // Initial polynomial degree. NOTE: The meaning is different from
                                         // standard continuous elements in the space H1. Here, P_INIT refers
                                         // to the maximum poly order of the tangential component, and polynomials
                                         // of degree P_INIT + 1 are present in element interiors. P_INIT = 0
                                         // is for Whitney elements.
const bool ALIGN_MESH = true;           // if ALIGN_MESH == true, curvilinear elements aligned with the
                                         // circular load are used, otherwise one uses a non-aligned mesh.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See the User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 2.0;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows
                                         // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
const double e_0   = 8.8541878176 * 1e-12;
const double mu_0   = 1.256 * 1e-6;
const double e_r = 1.0;
const double mu_r = 1.0;
const double rho = 3820.0;
const double Cp = 7.531000;
const double freq = 1.0*2450000000.0;
const double omega = 2 * M_PI * freq;
const double c = 1 / sqrt(e_0 * mu_0);
const double kappa  = 2 * M_PI * freq * sqrt(e_0 * mu_0);
const double J = 0.0000033333;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == 2) return BC_ESSENTIAL; // perfect conductor BC
  else return BC_NATURAL;               // impedance BC
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

// Geometry of the load.
bool in_load(double x, double y)
{
  double cx = -0.152994121;
  double cy =  0.030598824;
  double r = 0.043273273;
  if (sqr(cx - x) + sqr(cy - y) < sqr(r)) return true;
  else return false;
}

// Gamma as a function of x, y.
double gam(int marker, double x, double y)
{
  if (ALIGN_MESH && marker == 1) return 0.03;
  if (!ALIGN_MESH && in_load(x,y)) {
    double cx = -0.152994121;  double cy =  0.030598824;
    double r = sqrt(sqr(cx - x) + sqr(cy - y));
    return (0.03 + 1)/2.0 - (0.03 - 1) * atan(10.0*(r -  0.043273273)) / M_PI;
  }
  return 0.0;
}
double gam(int marker, Ord x, Ord y)
{  return 0.0; }

// Relative permittivity as a function of x, y.
double er(int marker, double x, double y)
{
  if (ALIGN_MESH && marker == 1) return 7.5;
  if (!ALIGN_MESH && in_load(x,y)) {
    double cx = -0.152994121;  double cy =  0.030598824;
    double r = sqrt(sqr(cx - x) + sqr(cy - y));
    return (7.5 + 1)/2.0 - (7.5 - 1) * atan(10.0*(r -  0.043273273)) / M_PI;
  }
  return 1.0;
}
double er(int marker, Ord x, Ord y)
{  return 1.0; }

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  if (ALIGN_MESH) mloader.load("oven_load_circle.mesh", &mesh);
  else mloader.load("oven_load_square.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Create an Hcurl space.
  HcurlSpace space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form_surf(callback(linear_form_surf));

  // Initialize views.
  VectorView eview("Electric field", 0, 0, 580, 400);
  OrderView ord("Order", 590, 0, 550, 400);

  /*
  // View the basis functions.
  VectorBaseView bview;
  vbview.show(&space);
  */

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Initialize refinements selector.
  HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  bool is_complex = true; 
  init_matrix_solver(matrix_solver, space.get_num_dofs(), mat, rhs, solver, is_complex);

  // Adaptivity loop:
  Solution sln, ref_sln;
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    info("Solving on reference mesh.");

    // Construct the globally refined reference mesh.
    Mesh ref_mesh;
    ref_mesh.copy(&mesh);
    ref_mesh.refine_all_elements();

    // Setup space for the reference solution.
    Space *ref_space = space.dup(&ref_mesh);
    int order_increase = 1;
    ref_space->copy_orders(&space, order_increase);
 
    // Solve the reference problem.
    solve_linear(ref_space, &wf, matrix_solver, &ref_sln, NULL, is_complex);

    // Project the reference solution on the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    // NULL means that we do not want to know the resulting coefficient vector.
    project_global(&space, H2D_HCURL_NORM, &ref_sln, &sln, NULL, is_complex);

    // Time measurement.
    cpu_time.tick();

    // Show real part of the solution.
    AbsFilter abs(&sln);
    eview.set_min_max_range(0, 4e3);
    eview.show(&abs);
    ord.show(&space);

    // Skip visualization time. 
    cpu_time.tick(HERMES_SKIP);

    // Calculate error estimate wrt. reference solution.
    info("Calculating error.");
    Adapt hp(&space, H2D_HCURL_NORM);
    hp.set_solutions(&sln, &ref_sln);
    hp.set_error_form(callback(hcurl_form_kappa));  // Adaptivity is done using an "energetic norm".
    double err_est_adapt = hp.calc_elem_errors() * 100;
    double err_est_hcurl = calc_abs_error(&sln, &ref_sln, H2D_HCURL_NORM);
    double norm_est = calc_norm(&ref_sln, H2D_HCURL_NORM); 
    double err_est_hcurl_rel = err_est_hcurl / norm_est * 100;

    // Report results.
    info("ndof: %d, ref_ndof: %d, err_est: %g%%", 
         space.get_num_dofs(), ref_space->get_num_dofs(), err_est_hcurl_rel);

    // Add entries to DOF convergence graphs.
    graph_dof_est.add_values(space.get_num_dofs(), err_est_hcurl_rel);
    graph_dof_est.save("conv_dof_est.dat");

    // Add entries to CPU convergence graphs.
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_hcurl_rel);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est_adapt too large, adapt the mesh.
    if (err_est_hcurl_rel < ERR_STOP) done = true;
    else {
      info("Adapting coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      if (space.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the fine mesh solution - the final result.
  eview.set_title("Final solution");
  eview.show(&ref_sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

