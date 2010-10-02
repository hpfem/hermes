#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This test makes sure that example 11-adapt-system works correctly.

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
  H1Space u_space(&u_mesh, bc_types, essential_bc_values, P_INIT_U);
  H1Space v_space(MULTI ? &v_mesh : &u_mesh, bc_types, essential_bc_values, P_INIT_V);

  // Initialize the weak formulation.
  WeakForm wf(2);
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0));
  wf.add_matrix_form(0, 1, callback(bilinear_form_0_1));
  wf.add_matrix_form(1, 0, callback(bilinear_form_1_0));
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1));
  wf.add_vector_form(0, linear_form_0, linear_form_0_ord);
  wf.add_vector_form(1, linear_form_1, linear_form_1_ord);

  // Initialize coarse and reference mesh solutions.
  Solution u_sln, v_sln, u_ref_sln, v_ref_sln;

  // Initialize exact solutions.
  ExactSolution u_exact(&u_mesh, uexact);
  ExactSolution v_exact(&v_mesh, vexact);

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  ScalarView s_view_0("Solution[0]", new WinGeom(0, 0, 440, 350));
  s_view_0.show_mesh(false);
  OrderView  o_view_0("Mesh[0]", new WinGeom(450, 0, 420, 350));
  ScalarView s_view_1("Solution[1]", new WinGeom(880, 0, 440, 350));
  s_view_1.show_mesh(false);
  OrderView  o_view_1("Mesh[1]", new WinGeom(1330, 0, 420, 350));

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est, 
              graph_dof_exact, graph_cpu_exact;

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Tuple<Space *>* ref_spaces = construct_refined_spaces(Tuple<Space *>(&u_space, &v_space));

    // Assemble the reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    FeProblem* fep = new FeProblem(&wf, *ref_spaces, is_linear);
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
    fep->assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();
    
    // Solve the linear system of the reference problem. If successful, obtain the solutions.
    if(solver->solve()) Solution::vector_to_solutions(solver->get_solution(), *ref_spaces, 
                                            Tuple<Solution *>(&u_ref_sln, &v_ref_sln));
    else error ("Matrix solver failed.\n");
  
    // Time measurement.
    cpu_time.tick();

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    project_global(Tuple<Space *>(&u_space, &v_space), Tuple<Solution *>(&u_ref_sln, &v_ref_sln), 
                   Tuple<Solution *>(&u_sln, &v_sln), matrix_solver); 
   
    // View the coarse mesh solution and polynomial orders.
    s_view_0.show(&u_sln); 
    o_view_0.show(&u_space);
    s_view_1.show(&v_sln); 
    o_view_1.show(&v_space);

    // Calculate element errors.
    info("Calculating error estimate and exact error."); 
    Adapt* adaptivity = new Adapt(Tuple<Space *>(&u_space, &v_space), Tuple<ProjNormType>(HERMES_H1_NORM, HERMES_H1_NORM));
    adaptivity->set_solutions(Tuple<Solution *>(&u_sln, &v_sln), Tuple<Solution *>(&u_ref_sln, &v_ref_sln));
    
    // Calculate error estimate for each solution component and the total error estimate.
    Tuple<double> err_est_rel;
    double err_est_rel_total = adaptivity->calc_err_est(err_est_rel, 
                               HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_ABS) * 100;

    // Calculate exact error for each solution component and the total exact error.
    Tuple<double> err_exact_rel;
    double err_exact_rel_total = adaptivity->calc_err_exact(err_exact_rel, 
                                 HERMES_TOTAL_ERROR_REL, Tuple<Solution *>(&u_exact, &v_exact)) * 100;

    // Time measurement.
    cpu_time.tick();

    // Report results.
    info("ndof_coarse[0]: %d, ndof_fine[0]: %d",
         u_space.get_num_dofs(), (*ref_spaces)[0]->get_num_dofs());
    info("err_est_rel[0]: %g%%, err_exact_rel[0]: %g%%", err_est_rel[0]*100, err_exact_rel[0]*100);
    info("ndof_coarse[1]: %d, ndof_fine[1]: %d",
         v_space.get_num_dofs(), (*ref_spaces)[1]->get_num_dofs());
    info("err_est_rel[1]: %g%%, err_exact_rel[1]: %g%%", err_est_rel[1]*100, err_exact_rel[1]*100);
    info("ndof_coarse_total: %d, ndof_fine_total: %d",
         get_num_dofs(Tuple<Space *>(&u_space, &v_space)), get_num_dofs(*ref_spaces));
    info("err_est_rel_total: %g%%, err_est_exact_total: %g%%", err_est_rel_total, err_exact_rel_total);

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(get_num_dofs(Tuple<Space *>(&u_space, &v_space)), err_est_rel_total);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel_total);
    graph_cpu_est.save("conv_cpu_est.dat");
    graph_dof_exact.add_values(get_num_dofs(Tuple<Space *>(&u_space, &v_space)), err_exact_rel_total);
    graph_dof_exact.save("conv_dof_exact.dat");
    graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel_total);
    graph_cpu_exact.save("conv_cpu_exact.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel_total < ERR_STOP) 
      done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(Tuple<RefinementSelectors::Selector *>(&selector, &selector), 
                               THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (get_num_dofs(Tuple<Space *>(&u_space, &v_space)) >= NDOF_STOP) done = true;

    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false)
      for(int i = 0; i < ref_spaces->size(); i++)
        delete (*ref_spaces)[i]->mesh;
    delete ref_spaces;
    delete fep;
    
    // Increase counter.
    as++;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  int ndof = get_num_dofs(Tuple<Space *>(&u_space, &v_space));

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

