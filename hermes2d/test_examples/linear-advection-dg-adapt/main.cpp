#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  This example solves a linear advection equation using Dicontinuous Galerkin (DG) method.
//  It is intended to show how evalutation of surface matrix forms that take basis functions defined
//  on different elements work. It is the same example as linear-advection-dg, but with automatic adaptivity.
//
//  PDE: \nabla \cdot (\Beta u) = 0, where \Beta = (-x_2, x_1) / |x| represents a circular counterclockwise flow field.
//
//  Domain: Square (0, 1) x (0, 1).
//
//  BC:		Dirichlet, u = 1 where \Beta(x) \cdot n(x) < 0, that is on [0,0.5] x {0}, and g = 0 anywhere else.
//				
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF = 2;                           
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 0;                             
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.2;                     
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
const int STRATEGY = 1;                           
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
// H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;          
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;                   
// Stopping criterion for adaptivity.
const double ERR_STOP = 3.0;                      
// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0. 
const double CONV_EXP = 1.0;                      
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;                      
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  
// Name of the iterative method employed by AztecOO (ignored by the other solvers). 
// Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* iterative_method = "bicgstab";        
// Name of the preconditioner employed by AztecOO (ignored by the other solvers). 
// Possibilities: none, jacobi, neumann, least-squares, or a
// preconditioner from IFPACK (see solver/aztecoo.h).
const char* preconditioner = "jacobi";            

int main(int argc, char* args[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF; i++) mesh.refine_all_elements();
  
  // Create an L2 space.
  L2Space<double> space(&mesh, P_INIT);
  
  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
  
  // Disable weighting of refinement candidates.
  selector.set_error_weights(1, 1, 1);

  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof_est, graph_cpu_est;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Display the mesh.
  OrderView oview("Coarse mesh", new WinGeom(0, 0, 440, 350));
  oview.show(&space);
  BaseView<double> bview("Distribution of polynomial orders", new WinGeom(450, 0, 440, 350));
  bview.show(&space);

  Solution<double> sln;
  Solution<double> ref_sln;

  // Initialize the weak formulation.
  CustomWeakForm wf("Bdy_bottom_left", &mesh);

  ScalarView view1("Solution", new WinGeom(900, 0, 450, 350));
  view1.fix_scale_width(60);

  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh
    // and setup reference space.
    Space<double>* ref_space = Space<double>::construct_refined_space(&space);

    info("Solving on reference mesh.");

    // Initialize the FE problem.
    DiscreteProblem<double>* dp = new DiscreteProblem<double>(&wf, ref_space);
    
    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix<double>* matrix = create_matrix<double>(matrix_solver);
    Vector<double>* rhs = create_vector<double>(matrix_solver);
    LinearMatrixSolver<double>* solver = create_linear_solver<double>(matrix_solver, matrix, rhs);

    // Assemble the linear problem.
    info("Assembling (ndof: %d).", Space<double>::get_num_dofs(ref_space));
    dp->assemble(matrix, rhs);

    // Solve the linear system. If successful, obtain the solution.
    info("Solving.");
    if(solver->solve()) Solution<double>::vector_to_solution(solver->get_sln_vector(), ref_space, &ref_sln);
    else error ("Matrix solver failed.\n");
    
    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection<double>::project_global(&space, &ref_sln, &sln, matrix_solver, HERMES_L2_NORM);  

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution.
    view1.show(&sln);
    oview.show(&space);
    bview.show(&space);
    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt<double>* adaptivity = new Adapt<double>(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;
 
    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
         Space<double>::get_num_dofs(&space), Space<double>::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof_est.add_values(Space<double>::get_num_dofs(&space), err_est_rel);
    graph_dof_est.save("conv_dof_est.dat");
    graph_cpu_est.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu_est.save("conv_cpu_est.dat");

    // If err_est_rel too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting the coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if (Space<double>::get_num_dofs(&space) >= NDOF_STOP) 
      {
        done = true;
        break;
      }
    }

    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false)
      delete ref_space->get_mesh();
    delete ref_space;
    delete dp;

    as++;
  }
  while (done == false);

  info("Total running time: %g s", cpu_time.accumulated());
  
  // Wait for keyboard or mouse input.
  View::wait();
  return 0;
}
