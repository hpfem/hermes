#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example shows how to run adaptive hp-FEM, h-FEM and p-FEM with
// basic control parameters. The underlying problem is a planar model
// of an electrostatic micromotor (MEMS). You may want to experiment with 
// various types of adaptivity via the options H2D_P_ISO, H2D_P_ANISO, 
// H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H, H2D_HP_ANISO_P, 
// and H2D_HP_ANISO. See the User Documentation for more details. 
//   Uniform initial polynomial degree of mesh elements can be set using
// the variable P_INIT. Before using adaptivity, you have to define a refinement 
// selector as shown below. The function adapt() takes the selector
// as a parameter, along with THRESHOLD, STRATEGY, and MESH_REGULARITY. 
//   Additional control parameters are available, these will be demonstrated
// in the following tutorial examples. In this example, two types of convergence  
// graphs are created -- error estimate wrt. the number of degrees of freedom 
// (DOF), and error estimate wrt. CPU time. Later we will show how to output 
// the error wrt. exact solution when exact solution is available. 
//   This example also demonstrates how to define different material parameters
// in various parts of the computational domain, and how to measure time. 
//
// PDE: -div[eps_r(x,y) grad phi] = 0
//      eps_r = EPS_1 in Omega_1 (surrounding air)
//      eps_r = EPS_2 in Omega_2 (moving part of the motor)
//
// BC: phi = 0 V on Gamma_1 (left edge and also the rest of the outer boundary
//     phi = VOLTAGE on Gamma_2 (boundary of stator)
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization. 
const bool VTK_OUTPUT = true;                     // Set to "true" to enable VTK output.
const int P_INIT = 2;                             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.2;                     // This is a quantitative parameter of the adapt(...) function and
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
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
                                                  // H2D_HP_ANISO_P, H2D_HP_ANISO. See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
                                                  // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const int OUTER_BDY = 1, STATOR_BDY = 2;

// Problem parameters.
const int MATERIAL_1 = 1;
const int MATERIAL_2 = 2;
const double EPS_1 = 1.0;       // Relative electric permittivity in Omega_1.
const double EPS_2 = 10.0;      // Relative electric permittivity in Omega_2.
const double VOLTAGE = 50.0;    // Voltage on the stator.

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("motor.mesh", &mesh);

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_dirichlet(Hermes::vector<int>(OUTER_BDY, STATOR_BDY));

  // Enter Dirichlet boundary values.
  BCValues bc_values;
  bc_values.add_const(STATOR_BDY, VOLTAGE);
  bc_values.add_const(OUTER_BDY, 0.0);

  // Create an H1 space with default shapeset.
  H1Space space(&mesh, &bc_types, &bc_values, P_INIT);

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(biform1), HERMES_SYM, MATERIAL_1);
  wf.add_matrix_form(callback(biform2), HERMES_SYM, MATERIAL_2);

  // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln;
  
  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

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

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = construct_refined_space(&space);

    // Initialize matrix solver.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Assemble reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space, is_linear);
    dp->assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();
    
    // Solve the linear system of the reference problem. 
    // If successful, obtain the solution.
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), ref_space, &ref_sln);
    else error ("Matrix solver failed.\n");

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

    // Time measurement.
    cpu_time.tick();
   
    // VTK output.
    if (VTK_OUTPUT) {
      // Output solution in VTK format.
      Linearizer lin;
      char* title = new char[100];
      sprintf(title, "sln-%d.vtk", as);
      lin.save_solution_vtk(&sln, title, "Potential", false);
      info("Solution in VTK format saved to file %s.", title);

      // Output mesh and element orders in VTK format.
      Orderizer ord;
      sprintf(title, "ord-%d.vtk", as);
      ord.save_orders_vtk(&space, title);
      info("Element orders in VTK format saved to file %s.", title);
    }

    // View the coarse mesh solution and polynomial orders.
    if (HERMES_VISUALIZATION) {
      sview.show(&sln);
      oview.show(&space);
    }

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt* adaptivity = new Adapt(&space);
    bool solutions_for_adapt = true;
    // In the following function, the Boolean parameter "solutions_for_adapt" determines whether 
    // the calculated errors are intended for use with adaptivity (this may not be the case, for example,
    // when error wrt. an exact solution is calculated). The default value is solutions_for_adapt = true, 
    // The last parameter "error_flags" determine whether the total and element errors are treated as 
    // absolute or relative. Its default value is error_flags = HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL.
    // In subsequent examples and benchmarks, these two parameters will be often used with 
    // their default values, and thus they will not be present in the code explicitly.
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, 
                         HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
      Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel);

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
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      
      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;

    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false) delete ref_space->get_mesh();
    delete ref_space;
    delete dp;
    
  }
  while (done == false);
  
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  sview.set_title("Fine mesh solution");
  sview.show_mesh(false);
  sview.show(&ref_sln);
  
  // Wait for all views to be closed.
  View::wait();

  return 0;
}

