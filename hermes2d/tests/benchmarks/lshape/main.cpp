#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

/** \addtogroup t_bench_lshape Benchmarks/LShape
 *  \{
 *  \brief This test makes sure that the benchmark "lshape" works correctly.
 *
 *  \section s_params Parameters
 *  - P_INIT=1
 *  - THERSHOLD=0.3
 *  - STRATEGY=0
 *  - CAND_LIST=HP_ANISO
 *  - MESH_REGULARITY=-1
 *  - ERR_STOP=0.1
 *  - CONV_EXP=1.0
 *  - NDOF_STOP=40000
 *  - ERROR_WEIGHTS=(H: 1; P: 1; ANISO: 1)
 *
 *  \section s_res Results
 *  - DOFs: 595
 */

const int P_INIT = 4;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial mesh refinements.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO_H;        // Predefined list of element refinement candidates. Possible values are
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
                                                  // error behavior err \approx const1*exp(-const2*pow(NDOF, CONV_EXP)).
const double ERR_STOP = 0.01;                     // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Exact solution.
#include "exact_solution.cpp"

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_ESSENTIAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return fn(x, y);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("lshape.mesh", &mesh);

  // Perform initial mesh refinement.
  mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space* space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form), H2D_SYM);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;  
  init_matrix_solver(matrix_solver, get_num_dofs(space), mat, rhs, solver);

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
    solve_linear(ref_space, &wf, matrix_solver, ref_sln);

    // Project the reference solution on the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    // NULL means that we do not want to know the resulting coefficient vector.
    project_global(space, H2D_H1_NORM, ref_sln, sln, NULL); 

    // Calculate element errors.
    info("Calculating error (est).");
    Adapt hp(space, H2D_H1_NORM);
    hp.set_solutions(sln, ref_sln);
    hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL);
 
    // Calculate error estimate for each solution component.
    double err_est_abs = calc_abs_error(sln, ref_sln, H2D_H1_NORM);
    double norm_est = calc_norm(ref_sln, H2D_H1_NORM);
    double err_est_rel = err_est_abs / norm_est * 100.;

    // Report results.
    info("ndof: %d, ref_ndof: %d, err_est_rel_total: %g%%", 
         get_num_dofs(space), get_num_dofs(ref_space), err_est_rel);

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

  int ndof = get_num_dofs(space);

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  int n_dof_allowed = 1400;
  printf("n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);
  if (ndof <= n_dof_allowed) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}

