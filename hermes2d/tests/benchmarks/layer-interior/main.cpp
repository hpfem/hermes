#include "hermes2d.h"

using namespace RefinementSelectors;

/** \addtogroup t_bench_layer_int Benchmarks/Layer-Interior
 *  \{
 *  \brief This test makes sure that the benchmark "layer-interior" works correctly.
 *
 *  \section s_params Parameters
 *   - INIT_REF_NUM=1
 *   - P_INIT=2
 *   - THRESHOLD=0.3
 *   - STRATEGY=0
 *   - CAND_LIST=H2D_HP_ANISO_H
 *   - MESH_REGULARITY=-1
 *   - CONV_EXP=1.0
 *   - ERR_STOP=1.0
 *   - NDOF_STOP=60000
 *   - matrix_solver = SOLVER_UMFPACK
 *   - SLOPE = 60
 *
 *  \section s_res Results
 *   - DOFs: 1433
 *   - Adaptivity steps: 13 
 */

const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
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
const double CONV_EXP = 0.5;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
double SLOPE = 60;                                // Slope of the layer.

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
  mloader.load("square_quad.mesh", &mesh);     // quadrilaterals
  // mloader.load("square_tri.mesh", &mesh);   // triangles

  // Perform initial mesh refinements.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form), H2D_SYM);
  wf.add_vector_form(callback(linear_form));

  // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  ExactSolution exact(&mesh, fndd);

  // Number of physical fields in the problem.
  int num_comps = Tuple<Space*>(&space).size();

  // Number of degreeso of freedom 
  int ndof = get_num_dofs(Tuple<Space*>(&space));

  // Number of exact solutions given.
  if (Tuple<ExactSolution *>(&exact).size() != 0 && Tuple<ExactSolution *>(&exact).size() != num_comps)
    error("Number of exact solutions does not match number of equations.");
  bool is_exact_solution;
  if (Tuple<ExactSolution *>(&exact).size() == num_comps) is_exact_solution = true;
  else is_exact_solution = false;

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu, graph_dof_exact, graph_cpu_exact;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = construct_refined_space(&space);

    // Assemble the reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    FeProblem* fep = new FeProblem(&wf, ref_space, is_linear);
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_solver(matrix_solver, matrix, rhs);
    fep->assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();

    // Solve the linear system of the reference problem. If successful, obtain the solution.
    info("Solving the matrix problem.");
    if(solver->solve()) vector_to_solution(solver->get_solution(), ref_space, &ref_sln);
    else error ("Matrix solver failed.\n");

    // Time measurement.
    cpu_time.tick();

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on the coarse mesh.");
    project_global(&space, H2D_H1_NORM, &ref_sln, &sln, matrix_solver);

    // Calculate element errors and total error estimate.
    info("Calculating error.");
    Adapt* adaptivity = new Adapt(&space, H2D_H1_NORM);
    adaptivity->set_solutions(&sln, &ref_sln);
    double err_est = adaptivity->calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;

    // Calculate exact error for each solution component.   
    double err_exact_abs[H2D_MAX_COMPONENTS];
    double norm_exact[H2D_MAX_COMPONENTS];
    double err_exact_abs_total = 0;
    double norm_exact_total = 0;
    double err_exact_rel_total;
    if (is_exact_solution == true) {
      err_exact_abs[0] = calc_abs_error(&sln, &exact, H2D_H1_NORM);
      norm_exact[0] = calc_norm(&exact, H2D_H1_NORM);
      err_exact_abs_total += err_exact_abs[0] * err_exact_abs[0];
      norm_exact_total += norm_exact[0] * norm_exact[0];

      err_exact_abs_total = sqrt(err_exact_abs_total);
      norm_exact_total = sqrt(norm_exact_total);
      err_exact_rel_total = err_exact_abs_total / norm_exact_total * 100.;
    }

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est: %g%%",
      get_num_dofs(&space), get_num_dofs(ref_space), err_est);
    if (is_exact_solution == true) info("err_exact_rel_total: %g%%", err_exact_rel_total);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(get_num_dofs(&space), err_est);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est);
    graph_cpu.save("conv_cpu_est.dat");
    if (is_exact_solution == true) {
      graph_dof_exact.add_values(get_num_dofs(&space), err_exact_rel_total);
      graph_dof_exact.save("conv_dof_exact.dat");
      graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel_total);
      graph_cpu_exact.save("conv_cpu_exact.dat");
    }

    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) done = true;
    else
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    if (get_num_dofs(&space) >= NDOF_STOP) done = true;

    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false)
      delete ref_space->mesh;
    delete ref_space;
    delete fep;

  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  ndof = get_num_dofs(&space);

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  int n_dof_allowed = 1450;
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

