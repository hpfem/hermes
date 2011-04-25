#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

//  This example shows how to save visualization data if you are working 
//  on a distant computer and cannot use ScalarView, OrderView, or 
//  related classes directly. The two basic ways are:
//    * to use the methods Solution::save(), which saves a complete 
//      Solution including Mesh and element orders. Then you can fetch the 
//      file and use Solution::load() to restore the Solution completely 
//      in the memory of your local machine. 
//    * Linearizer::save_data() which only saves linearized data for direct 
//      OpenGL processing. After fetching the file, you can use the methods
//      ScalarView::Linearizer::load_data() and ScalarView::show_linearizer_data()
//      on your local machine.
//  The underlying model for computation is the tutorial example P03-timedep/01-cathedral-ie. 

int OUTPUT_FREQUENCY = 20;                        // Number of time steps between saving data.

const int P_INIT = 2;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 3;                   // Number of initial uniform mesh refinements towards the boundary.
const double time_step = 300.0;                   // Time step in seconds.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double TEMP_INIT = 10;     // Temperature of the ground (also initial temperature).
const double ALPHA = 10;         // Heat flux coefficient for Newton's boundary condition.
const double LAMBDA = 1e2;       // Thermal conductivity of the material.
const double HEATCAP = 1e2;      // Heat capacity.
const double RHO = 3000;         // Material density.
const double T_FINAL = 18000;    // Length of time interval in seconds.

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Hermes2D hermes2d;

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("cathedral.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary("Boundary air", INIT_REF_NUM_BDY);
  mesh.refine_towards_boundary("Boundary ground", INIT_REF_NUM_BDY);

  // Previous time level solution (initialized by the external temperature).
  Solution tsln(&mesh, TEMP_INIT);

  // Initialize the weak formulation.
  double current_time = 0;
  CustomWeakFormHeatRK1 wf("Boundary air", ALPHA, LAMBDA, HEATCAP, RHO, time_step, 
                           &current_time, TEMP_INIT, T_FINAL, &tsln);

  // Initialize boundary conditions.
  DefaultEssentialBCConst essential_bc("Boundary ground", TEMP_INIT);
  EssentialBCs bcs(&essential_bc);

  // Initialize an H1 space with default shepeset.
  H1Space space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Initialize the FE problem.
  DiscreteProblem dp(&wf, &space);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);

  // Initial coefficient vector for the Newton's method.  
  scalar* coeff_vec = new scalar[ndof];
  memset(coeff_vec, 0, ndof*sizeof(scalar));

  // Initialize views.
  ScalarView Tview("Temperature", new WinGeom(0, 0, 450, 600));
  char title[100];
  sprintf(title, "Time %3.5f", current_time);
  Tview.set_min_max_range(0,20);
  Tview.set_title(title);
  Tview.fix_scale_width(3);

  // Time stepping:
  int nsteps = (int)(T_FINAL/time_step + 0.5);
  bool jacobian_changed = true;
  for(int ts = 1; ts <= nsteps; ts++)
  {
    info("---- Time step %d, time %3.5f", ts, current_time);

    // Perform Newton's iteration.
    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs, 
        jacobian_changed)) error("Newton's iteration failed.");
    jacobian_changed = false;

    // Translate the resulting coefficient vector into the Solution sln.
    Solution::vector_to_solution(coeff_vec, &space, &tsln);

    if (ts % OUTPUT_FREQUENCY == 0) {
      Linearizer lin;
      int item = H2D_FN_VAL_0;
      double eps = HERMES_EPS_NORMAL;
      double max_abs = -1.0;
      MeshFunction* xdisp = NULL; 
      MeshFunction* ydisp = NULL;
      double dmult = 1.0;
      lin.process_solution(&tsln, item, eps, max_abs, xdisp, ydisp, dmult);
      char* filename = new char[100];
      sprintf(filename, "tsln_%d.lin", ts);

      // Save Linearizer data.
      lin.save_data(filename);
      info("Linearizer data saved to file %s.", filename);

      // Save complete Solution.
      sprintf(filename, "tsln_%d.dat", ts);
      bool compress = false;   // Gzip compression not used as it only works on Linux.
      tsln.save(filename, compress);
      info("Complete Solution saved to file %s.", filename);
    }

    // Update the time variable.
    current_time += time_step;
  }

  info("Let's assume that the remote computation has finished and you fetched the *.lin files.");
  info("Visualizing Linearizer data from file tsln_40.lin.");

  // First use ScalarView to read and show the Linearizer data.
  WinGeom* win_geom_1 = new WinGeom(0, 0, 450, 600);
  ScalarView sview_1("Saved Linearizer data", win_geom_1);
  sview_1.lin.load_data("tsln_40.lin");
  sview_1.set_min_max_range(0,20);
  sview_1.fix_scale_width(3);
  sview_1.show_linearizer_data();

  info("Visualizing Solution from file tsln_60.dat.");

  Solution sln_from_file;
  sln_from_file.load("tsln_60.dat");
  WinGeom* win_geom_2 = new WinGeom(460, 0, 450, 600);
  ScalarView sview_2("Saved Solution data", win_geom_2);
  sview_2.set_min_max_range(0,20);
  sview_2.fix_scale_width(3);
  sview_2.show(&sln_from_file);

  info("Visualizing Mesh and Orders extracted from the Solution.");
 
  int p_init = 1;

  H1Space space_from_file(sln_from_file.get_mesh(), p_init);
  space_from_file.set_element_orders(sln_from_file.get_element_orders());
  OrderView oview("Saved Solution -> Space", new WinGeom(920, 0, 450, 600));
  oview.show(&space_from_file);

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  // Wait for the view to be closed.
  View::wait();
  return 0;
}
