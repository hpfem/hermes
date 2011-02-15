#define HERMES_REPORT_ALL
#include "hermes2d.h"

// This example shows how to use the automatic h-adaptivity based on a Kelly-type error estimator for an elliptic
// problem with all types of standard boundary conditions.
//
// PDE: Laplace equation -Laplace u = f, where f = CONST_F.
//
// BC: u = T1 ... fixed temperature on Gamma_left (Dirichlet)
//     du/dn = 0 ... insulated wall on Gamma_inner (Neumann)
//     du/dn = C ... constant heat flux on Gamma_inner (Neumann)
//     du/dn = H*(u - T0) ... heat flux on Gamma_bottom (Newton)
//
// Note that the last BC can be written in the form  du/dn - H*u = -H*T0.
//
// The following parameters can be changed:

const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int CORNER_REF_LEVEL = 0;                  // Number of mesh refinements towards the re-entrant corner.
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
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double ERR_STOP = 0.1;                      // Stopping criterion for adaptivity (rel. error tolerance between the

const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Element material markers and associated diffusion coefficients.
const int MATERIAL_1 = 1; const double EPS_1 = 1.0;
const int MATERIAL_2 = 2; const double EPS_2 = 10.0;

// Boundary markers and boundary condition data.
const int BDY_BOTTOM = 1, BDY_OUTER = 2, BDY_LEFT = 3, BDY_INNER = 4;
const double T1 = 30.0;                   // Prescribed temperature on Gamma_left.
const double T0 = 20.0;                   // Outer temperature on Gamma_bottom.
const double H  = 0.05;                   // Heat flux on Gamma_bottom.
const double CONST_GAMMA_OUTER = 1.00;    // Heat flux on Gamma_inner.

// Right-hand side.
const double CONST_F = -1.0;              

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i=0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_vertex(3, CORNER_REF_LEVEL);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(BDY_LEFT);
  bc_types.add_bc_neumann(Hermes::vector<int>(BDY_OUTER, BDY_INNER));
  bc_types.add_bc_newton(BDY_BOTTOM);

  // Enter Dirichlet boudnary values.
  BCValues bc_values;
  bc_values.add_const(BDY_LEFT, T1);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form_1), HERMES_SYM, MATERIAL_1);
  wf.add_matrix_form(callback(bilinear_form_2), HERMES_SYM, MATERIAL_2);
  wf.add_matrix_form_surf(callback(bilinear_form_surf_bottom), BDY_BOTTOM);
  wf.add_vector_form_surf(callback(linear_form_surf_bottom), BDY_BOTTOM);
  wf.add_vector_form_surf(callback(linear_form_surf_outer), BDY_OUTER);

  // Initialize coarse and reference mesh solution.
  Solution sln;
  
  // Initialize views.
  ScalarView sview("Solution", new WinGeom(0, 0, 410, 600));
  sview.fix_scale_width(50);
  sview.show_mesh(false);
  OrderView  oview("Polynomial orders", new WinGeom(420, 0, 400, 600));

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();
  
  // Initialize matrix solver.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);
       
    // Assemble reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, &space, is_linear);
    dp->assemble(matrix, rhs);
    
    // Time measurement.
    cpu_time.tick();
    
    // Solve the linear system of the reference problem. 
    // If successful, obtain the solution.
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
    else error ("Matrix solver failed.\n");
    
    // Time measurement.
    cpu_time.tick();
    
    // View the coarse mesh solution and polynomial orders.
    sview.show(&sln);
    oview.show(&space);
    
    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);
    
    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    KellyTypeAdapt* adaptivity = new KellyTypeAdapt(&space);
    adaptivity->add_error_estimator_surf(callback(kelly_interface_estimator));
    adaptivity->add_error_estimator_surf(callback(kelly_newton_boundary_estimator), BDY_BOTTOM);
    adaptivity->add_error_estimator_surf(callback(kelly_neumann_boundary_estimator), BDY_OUTER);
    adaptivity->add_error_estimator_surf(callback(kelly_zero_neumann_boundary_estimator), BDY_INNER);
    
    double err_est_rel = adaptivity->calc_err_est(&sln, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;
                                                  
    // Report results.
    info("ndof: %d, err_est_rel: %g%%", Space::get_num_dofs(&space), err_est_rel);
    
    // Time measurement.
    cpu_time.tick();
    
    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");
    
    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;
    
    // Increase the counter of performed adaptivity steps.
    if (done == false)  as++;
    
    // Clean up.
    delete adaptivity;
    delete dp;                                    
  }
  while (done == false);
  
  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;
  
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  // Wait for all views to be closed.
  View::wait();
  
  return 0;
}
