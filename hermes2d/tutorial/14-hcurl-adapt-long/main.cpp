#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This is a long version of example 14-hcurl-adapt: function solve_linear_adapt() is not used.

const int P_INIT = 2;                             // Initial polynomial degree. NOTE: The meaning is different from
                                                  // standard continuous elements in the space H1. Here, P_INIT refers
                                                  // to the maximum poly order of the tangential component, and polynomials
                                                  // of degree P_INIT + 1 are present in element interiors. P_INIT = 0
                                                  // is for Whitney elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
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
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
const double mu_r   = 1.0;
const double kappa  = 1.0;
const double lambda = 1.0;

// Bessel functions, exact solution, and weak forms.
#include "forms.cpp"

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == 1 || marker == 6)
    return BC_ESSENTIAL; // perfect conductor
  else
    return BC_NATURAL;   // impedance
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("lshape3q.mesh", &mesh);    // quadrilaterals
  //mloader.load("lshape3t.mesh", &mesh);  // triangles

  // Perform initial mesh refinemets.
  for (int i=0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Create an Hcurl space with default shapeset.
  HcurlSpace* space = new HcurlSpace(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form), H2D_SYM);
  wf.add_matrix_form_surf(callback(bilinear_form_surf));
  wf.add_vector_form_surf(linear_form_surf, linear_form_surf_ord);

  // Initialize refinement selector.
  HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Initialize matrix solver.
  bool is_complex = true;
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, get_num_dofs(space), mat, rhs, solver, is_complex);

  // Initialize views.
  VectorView v_view("Solution (magnitude)", new WinGeom(0, 0, 440, 350));
  OrderView  o_view("Mesh", new WinGeom(450, 0, 400, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Adaptivity loop.
  Solution *sln = new Solution();
  Solution *ref_sln = new Solution();
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    info("Solving on reference mesh.");

    // Construct globally refined reference mesh
    // and setup reference space.
    Mesh *ref_mesh = new Mesh();
    ref_mesh->copy(space->get_mesh());
    ref_mesh->refine_all_elements();
    Space* ref_space = space->dup(ref_mesh);
    int order_increase = 1;
    ref_space->copy_orders(space, order_increase);

    // Solve the reference problem.
    // The NULL pointer means that we do not want the resulting coefficient vector.
    solve_linear(ref_space, &wf, matrix_solver, ref_sln, NULL, is_complex);

    // Project the reference solution on the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    // NULL means that we do not want to know the resulting coefficient vector.
    project_global(space, H2D_HCURL_NORM, ref_sln, sln, NULL, is_complex); 

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution.
    // NOTE: It is possible to visualize the real and imaginary parts,
    // see RealFilter and ImagFilter in filter.h.
    v_view.show(sln); 
    o_view.show(space);

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors.
    info("Calculating error (est).");
    Adapt hp(space, H2D_HCURL_NORM);
    hp.set_solutions(sln, ref_sln);
    double err_est_rel = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100.;
 
    // Report results.
    info("ndof: %d, ref_ndof: %d, err_est_rel_total: %g%%", 
         get_num_dofs(space), get_num_dofs(ref_space), err_est_rel);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(get_num_dofs(space), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else {
      info("Adapting the coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (get_num_dofs(space) >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);

  info("Total running time: %g s", cpu_time.accumulated());

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

