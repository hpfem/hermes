#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This is a long version of test example "11-adapt-long": function solve_linear_adapt() is not used.
// This test makes sure that example 11-adapt-system-long works correctly.

const int P_INIT_U = 2;                  // Initial polynomial degree for u.
const int P_INIT_V = 2;                  // Initial polynomial degree for v.
const int INIT_REF_BDY = 3;              // Number of initial boundary refinements
const bool MULTI = true;                 // MULTI = true  ... use multi-mesh,
                                         // MULTI = false ... use single-mesh.
                                         // Note: In the single mesh option, the meshes are
                                         // forced to be geometrically the same but the
                                         // polynomial degrees can still vary.
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                  // Adaptive strategy:
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
const int MESH_REGULARITY = -1;  // Maximum allowed level of hanging nodes:
                                 // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                 // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                 // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                 // Note that regular meshes are not supported, this is due to
                                 // their notoriously bad performance.
const double CONV_EXP = 1;       // Default value is 1.0. This parameter influences the selection of
                                 // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const int MAX_ORDER = 10;        // Maximum allowed element degree
const double ERR_STOP = 0.5;     // Stopping criterion for adaptivity (rel. error tolerance between the
                                 // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;     // Adaptivity process stops when the number of degrees of freedom grows over
                                 // this limit. This is mainly to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
const double D_u = 1;
const double D_v = 1;
const double SIGMA = 1;
const double LAMBDA = 1;
const double KAPPA = 1;
const double K = 100;

// Boundary condition types.
BCType bc_types(int marker) { return BC_ESSENTIAL; }

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y) { return 0;}

// Exact solution.
#include "exact_solution.cpp"

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh u_mesh, v_mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &u_mesh);
  if (MULTI == false) u_mesh.refine_towards_boundary(1, INIT_REF_BDY);

  // Create initial mesh (master mesh).
  v_mesh.copy(&u_mesh);

  // Initial mesh refinements in the v_mesh towards the boundary.
  if (MULTI == true) v_mesh.refine_towards_boundary(1, INIT_REF_BDY);

  // Create H1 spaces with default shapeset for both displacement components.
  H1Space* u_space = new H1Space(&u_mesh, bc_types, essential_bc_values, P_INIT_U);
  H1Space* v_space = new H1Space(MULTI ? &v_mesh : &u_mesh, bc_types, essential_bc_values, P_INIT_V);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1));
  wf.add_matrix_form(1, 0, callback(bilinear_form_1_0));
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1));
  wf.add_vector_form(0, linear_form_0, linear_form_0_ord);
  wf.add_vector_form(1, linear_form_1, linear_form_1_ord);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, get_num_dofs(Tuple<Space *>(u_space, v_space)), mat, rhs, solver);

  // Adaptivity loop.
  Solution *u_sln = new Solution();
  Solution *v_sln = new Solution();
  Solution *u_ref_sln = new Solution();
  Solution *v_ref_sln = new Solution();
  ExactSolution u_exact(&u_mesh, uexact);
  ExactSolution v_exact(&v_mesh, vexact);
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
    info("Solving on reference mesh.");

    // Construct globally refined reference mesh
    // and setup reference space.
    Mesh *u_ref_mesh = new Mesh();
    u_ref_mesh->copy(u_space->get_mesh());
    u_ref_mesh->refine_all_elements();
    Space* u_ref_space = u_space->dup(u_ref_mesh);
    int order_increase = 1;
    u_ref_space->copy_orders(u_space, order_increase);
    Mesh *v_ref_mesh = new Mesh();
    v_ref_mesh->copy(v_space->get_mesh());
    v_ref_mesh->refine_all_elements();
    Space* v_ref_space = v_space->dup(v_ref_mesh);
    v_ref_space->copy_orders(v_space, order_increase);

    // Solve the reference problem.
    solve_linear(Tuple<Space *>(u_ref_space, v_ref_space), &wf, matrix_solver, 
                 Tuple<Solution *>(u_ref_sln, v_ref_sln));

    // Project the reference solution on the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    // NULL means that we do not want to know the resulting coefficient vector.
    project_global(Tuple<Space *>(u_space, v_space), 
                   Tuple<int>(H2D_H1_NORM, H2D_H1_NORM), 
                   Tuple<MeshFunction *>(u_ref_sln, v_ref_sln), 
                   Tuple<Solution *>(u_sln, v_sln), NULL); 

    // Calculate element errors.
    info("Calculating error (est).");
    Adapt hp(Tuple<Space *>(u_space, v_space), 
             Tuple<int>(H2D_H1_NORM, H2D_H1_NORM));
    hp.set_solutions(Tuple<Solution *>(u_sln, v_sln), 
                     Tuple<Solution *>(u_ref_sln, v_ref_sln));
    hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL);
 
    // Calculate error estimate for each solution component.
    double u_err_est_abs = calc_abs_error(u_sln, u_ref_sln, H2D_H1_NORM);
    double u_norm_est = calc_norm(u_ref_sln, H2D_H1_NORM);
    double v_err_est_abs = calc_abs_error(v_sln, v_ref_sln, H2D_H1_NORM);
    double v_norm_est = calc_norm(v_ref_sln, H2D_H1_NORM);
    double err_est_abs_total = sqrt(u_err_est_abs*u_err_est_abs + v_err_est_abs*v_err_est_abs);
    double norm_est_total = sqrt(u_norm_est*u_norm_est + v_norm_est*v_norm_est);
    double err_est_rel_total = err_est_abs_total / norm_est_total * 100.;

    // Calculate exact error for each solution component.   
    double err_exact_abs_total = 0;
    double norm_exact_total = 0;
    double u_err_exact_abs = calc_abs_error(u_sln, &u_exact, H2D_H1_NORM);
    double u_norm_exact = calc_norm(&u_exact, H2D_H1_NORM);
    err_exact_abs_total += u_err_exact_abs * u_err_exact_abs;
    norm_exact_total += u_norm_exact * u_norm_exact;
    double v_err_exact_abs = calc_abs_error(v_sln, &v_exact, H2D_H1_NORM);
    double v_norm_exact = calc_norm(&v_exact, H2D_H1_NORM);
    err_exact_abs_total += v_err_exact_abs * v_err_exact_abs;
    norm_exact_total += v_norm_exact * v_norm_exact;
    err_exact_abs_total = sqrt(err_exact_abs_total);
    norm_exact_total = sqrt(norm_exact_total);
    double err_exact_rel_total = err_exact_abs_total / norm_exact_total * 100.;

    // Report results.
    info("ndof[0]: %d, ref_ndof[0]: %d, err_est_rel[0]: %g%%", 
         u_space->get_num_dofs(), u_ref_space->get_num_dofs(),
         u_err_est_abs/u_norm_est*100);
    info("err_exact_rel[0]: %g%%", u_err_exact_abs/u_norm_exact*100);
    info("ndof[1]: %d, ref_ndof[1]: %d, err_est_rel[1]: %g%%", 
         v_space->get_num_dofs(), v_ref_space->get_num_dofs(),
         v_err_est_abs/v_norm_est*100);
    info("err_exact_rel[1]: %g%%", v_err_exact_abs/v_norm_exact*100);
    info("ndof: %d, ref_ndof: %d, err_est_rel_total: %g%%", 
         get_num_dofs(Tuple<Space *>(u_space, v_space)), 
         get_num_dofs(Tuple<Space *>(u_ref_space, v_ref_space)), err_est_rel_total);

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < ERR_STOP) done = true;
    else {
      info("Adapting the coarse mesh.");
      done = hp.adapt(Tuple<RefinementSelectors::Selector *>(&selector, &selector), 
                      THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (get_num_dofs(Tuple<Space *>(u_space, v_space)) >= NDOF_STOP) done = true;
    }

    as++;
  }
  while (done == false);

  int ndof = get_num_dofs(Tuple<Space *>(u_space, v_space));

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  printf("ndof allowed = %d\n", 1200);
  printf("ndof actual = %d\n", ndof);
  if (ndof < 1200) {      // ndofs was 1158 at the time this test was created
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

