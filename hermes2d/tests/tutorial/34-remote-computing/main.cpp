#define H2D_REPORT_INFO
#define H2D_REPORT_FILE "application.log"
#include "hermes2d.h"

// This test makes sure that example 34-remote-computing works correctly.

int OUTPUT_FREQUENCY = 20;                        // Number of time steps between saving data.

const int P_INIT = 4;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 1;                   // Number of initial uniform mesh refinements towards the boundary.
const double TAU = 300.0;                         // Time step in seconds.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_UMFPACK, SOLVER_PETSC,
                                                  // SOLVER_MUMPS, and more are coming.

// Problem parameters.
const double T_INIT = 10;        // Temperature of the ground (also initial temperature).
const double ALPHA = 10;         // Heat flux coefficient for Newton's boundary condition.
const double LAMBDA = 1e5;       // Thermal conductivity of the material.
const double HEATCAP = 1e6;      // Heat capacity.
const double RHO = 3000;         // Material density.
const double FINAL_TIME = 18000; // Length of time interval (24 hours) in seconds.

// Global time variable.
double TIME = 0;

// Time-dependent exterior temperature.
double temp_ext(double t) {
  return T_INIT + 10. * sin(2*M_PI*t/FINAL_TIME);
}

// Boundary markers.
int bdy_ground = 1;
int bdy_air = 2;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (marker == bdy_ground) return BC_ESSENTIAL;
  else return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return T_INIT;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("cathedral.mesh", &mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(bdy_air, INIT_REF_NUM_BDY);

  // Initialize an H1 space with default shepeset.
  H1Space space(&mesh, bc_types, essential_bc_values, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d.", ndof);

  // Set initial condition.
  Solution tsln;
  tsln.set_const(&mesh, T_INIT);

  // Initialize weak formulation.
  WeakForm wf;
  wf.add_matrix_form(bilinear_form<double, double>, bilinear_form<Ord, Ord>);
  wf.add_matrix_form_surf(bilinear_form_surf<double, double>, bilinear_form_surf<Ord, Ord>, bdy_air);
  wf.add_vector_form(linear_form<double, double>, linear_form<Ord, Ord>, HERMES_ANY, &tsln);
  wf.add_vector_form_surf(linear_form_surf<double, double>, linear_form_surf<Ord, Ord>, bdy_air);

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, &space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Initialize views.
  ScalarView Tview("Temperature", new WinGeom(0, 0, 450, 600));
  char title[100];
  sprintf(title, "Time %3.5f, exterior temperature %3.5f", TIME, temp_ext(TIME));
  Tview.set_min_max_range(0,20);
  Tview.set_title(title);
  Tview.fix_scale_width(3);

  // Time stepping:
  int nsteps = (int)(FINAL_TIME/TAU + 0.5);
  bool rhs_only = false;
  for(int ts = 1; ts <= nsteps; ts++)
  {
    info("---- Time step %d, time %3.5f, ext_temp %g", ts, TIME, temp_ext(TIME));

    // First time assemble both the stiffness matrix and right-hand side vector,
    // then just the right-hand side vector.
    if (rhs_only == false) info("Assembling the stiffness matrix and right-hand side vector.");
    else info("Assembling the right-hand side vector (only).");
    dp.assemble(matrix, rhs, rhs_only);
    rhs_only = true;

    // Solve the linear system and if successful, obtain the solution.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution::vector_to_solution(solver->get_solution(), &space, &tsln);
    else 
      error ("Matrix solver failed.\n");

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
    TIME += TAU;
  }

  info("Let's assume that the remote computation has finished and you fetched the *.lin files.");
  info("Visualizing Linearizer data from file tsln_40.lin.");

  // First use ScalarView to read and show the Linearizer data.
  WinGeom* win_geom_1 = new WinGeom(0, 0, 450, 600);
  ScalarView sview_1("Saved Linearizer data", win_geom_1);
  sview_1.lin.load_data("tsln_40.lin");
  sview_1.set_min_max_range(0,20);
  sview_1.fix_scale_width(3);
  //sview_1.show_linearizer_data();

  info("Visualizing Solution from file tsln_60.dat.");

  Solution sln_from_file;
  sln_from_file.load("tsln_60.dat");
  WinGeom* win_geom_2 = new WinGeom(460, 0, 450, 600);
  ScalarView sview_2("Saved Solution data", win_geom_2);
  sview_2.set_min_max_range(0,20);
  sview_2.fix_scale_width(3);

  info("Visualizing Mesh and Orders extracted from the Solution.");
 
  int p_init = 1;
  // The NULLs are for bc_types() and essential_bc_values().
  H1Space space_from_file(sln_from_file.get_mesh(), NULL, NULL, p_init);
  space_from_file.set_element_orders(sln_from_file.get_element_orders());
  WinGeom* win_geom_3 = new WinGeom(920, 0, 450, 600);
  OrderView oview("Saved Solution -> Space", win_geom_3);

  // Wait for the view to be closed.
  // View::wait();

  // It will "Exception: SegFault" if we do not use View::wait() or View::close(). 
  sview_1.close(); 
  sview_2.close(); 
  oview.close(); 

  bool success = true;

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1
  if (success == true) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }

  delete win_geom_1;
  delete win_geom_2;
  delete win_geom_3;
}
