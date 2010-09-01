#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "function.h"

using namespace RefinementSelectors;

// This is a long version of example 17-newton-elliptic-adapt: function solve_newton_adapt() is not used.

const int P_INIT = 1;                             // Initial polynomial degree.
const int INIT_GLOB_REF_NUM = 1;                  // Number of initial uniform mesh refinements.
const int INIT_BDY_REF_NUM = 3;                   // Number of initial refinements towards boundary.

const double THRESHOLD = 0.2;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See the User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
const double NEWTON_TOL_COARSE = 1e-4;            // Stopping criterion for the Newton's method on coarse mesh.
const double NEWTON_TOL_FINE = 1e-4;              // Stopping criterion for the Newton's method on fine mesh.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Thermal conductivity (temperature-dependent).
// Note: for any u, this function has to be positive.
template<typename Real>
Real lam(Real u)
{
  return 1 + pow(u, 4);
}

// Derivative of the thermal conductivity with respect to 'u'.
template<typename Real>
Real dlam_du(Real u) {
  return 4*pow(u, 3);
}

// This function is used to define Dirichlet boundary conditions.
double dir_lift(double x, double y, double& dx, double& dy) {
  dx = (y+10)/10.;
  dy = (x+10)/10.;
  return (x+10)*(y+10)/100.;
}

// Initial condition. It will be projected on the FE mesh 
// to obtain initial coefficient vector for the Newton's method.
scalar init_cond(double x, double y, double& dx, double& dy)
{
  // Using the Dirichlet lift elevated by two
  double val = dir_lift(x, y, dx, dy) + 2;
  return val;
}

// Boundary condition types.
BCType bc_types(int marker) { return BC_ESSENTIAL;}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  double dx, dy;
  return dir_lift(x, y, dx, dy);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(1, INIT_BDY_REF_NUM);

  // Create an H1 space with default shapeset.
  H1Space* space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(jac), H2D_UNSYM, H2D_ANY);
  wf.add_vector_form(callback(res), H2D_ANY);

  // Create a selector which will select optimal candidate.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView s_view("Solution", new WinGeom(0, 0, 440, 350));
  s_view.show_mesh(false);
  OrderView o_view("Mesh", new WinGeom(450, 0, 400, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Project the initial condition on the FE space to obtain initial 
  // coefficient vector for the Newton's method.
  info("Projecting initial condition to obtain initial vector on the coarse mesh.");
  // The NULL pointer means that we do not want the projection result as a Solution.
  Vector* coeff_vec = new AVector();
  Solution* sln_tmp = new Solution(&mesh, init_cond);
  project_global(space, H2D_H1_NORM, sln_tmp, NULL, coeff_vec); 
  delete sln_tmp;

  // Newton's loop on the coarse mesh.
  info("Solving on coarse mesh:");
  Solution* sln = new Solution();
  Solution* ref_sln = new Solution();
  bool verbose = true;
  if (!solve_newton(space, &wf, coeff_vec, matrix_solver, 
		    NEWTON_TOL_COARSE, NEWTON_MAX_ITER, verbose)) {
    error("Newton's method did not converge.");
  }

  // Store the result in sln.
  sln->set_fe_solution(space, coeff_vec);

  // Adaptivity loop:
  bool done = false; int as = 1;
  double err_est;
  do {
    info("---- Adaptivity step %d:", as);

    // Time measurement..
    cpu_time.tick();

    // View the coarse mesh solution.
    s_view.show(sln);
    o_view.show(space);

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Construct globally refined reference mesh
    // and setup reference space.
    Space* ref_space;
    Mesh *ref_mesh = new Mesh();
    ref_mesh->copy(space->get_mesh());
    ref_mesh->refine_all_elements();
    ref_space = space->dup(ref_mesh);
    int order_increase = 1;
    ref_space->copy_orders(space, order_increase);

    // Calculate initial coefficient vector for Newton on the fine mesh.
    if (as == 1) {
      info("Projecting coarse mesh solution to obtain initial vector on new fine mesh.");
      // The NULL means that we do not want the result as a Solution.
      project_global(ref_space, H2D_H1_NORM, sln, NULL, coeff_vec);
    }
    else {
      info("Projecting previous fine mesh solution to obtain initial vector on new fine mesh.");
      // The NULL means that we do not want the result as a Solution.
      project_global(ref_space, H2D_H1_NORM, ref_sln, NULL, coeff_vec);
    }

    // Newton's method on fine mesh
    info("Solving on fine mesh.");
    if (!solve_newton(ref_space, &wf, coeff_vec, matrix_solver, NEWTON_TOL_FINE, NEWTON_MAX_ITER, verbose))
      error("Newton's method did not converge.");

    // Store the result in ref_sln.
    ref_sln->set_fe_solution(ref_space, coeff_vec);

    // Calculate element errors.
    info("Calculating error (est).");
    Adapt hp(space, H2D_H1_NORM);
    // Pass coarse mesh and reference solutions for error estimation.
    hp.set_solutions(sln, ref_sln);
    double err_est_rel_total = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100.;

    // Report results.
    info("ndof: %d, ref_ndof: %d, err_est_rel: %g%%", 
         get_num_dofs(space), get_num_dofs(ref_space), err_est_rel_total);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(get_num_dofs(space), err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < ERR_STOP) done = true;
    else {
      if (verbose) info("Adapting the coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (get_num_dofs(space) >= NDOF_STOP) {
        done = true;
        break;
      }
      
      // Project last fine mesh solution on the new coarse mesh
      // to obtain new coars emesh solution.
      info("Projecting reference solution on new coarse mesh for error calculation.");
      // NULL means that we do not want to know the resulting coefficient vector.
      project_global(space, H2D_H1_NORM, ref_sln, sln, NULL); 
    }

    // Free the reference space and mesh.
    ref_space->free();
    ref_mesh->free();

    as++;
  }
  while (done == false);
  if (verbose) info("Total running time: %g s", cpu_time.accumulated());

  delete coeff_vec;
  
  // Wait for keyboard or mouse input.
  View::wait();
  return true;
}


