#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "disc.h"
#include "basic.h"

int main(int argc, char* argv[])
{

  /*** RECEIVE DATA ***/

  // Check the number of command line arguments.
  if(argc != 2) error("Configuration file missing.");

  // Open configuration file.
  FILE* f = fopen(argv[1], "r");
  if(f == NULL) error("Cannot open file %s.", argv[1]);

  // Initialize the ModuleBasic class.
  ModuleBasic B;

  /*** READ MESH (TEMPORARILY AS A STRING )***/

  // Set mesh filename.
  B.set_mesh_str("\na = 1.0  # size of the mesh\nb = sqrt(2)/2\n\nvertices =\n{\n  { 0, -a },    # vertex 0\n  { a, -a },    # vertex 1\n  { -a, 0 },    # vertex 2\n  { 0, 0 },     # vertex 3\n  { a, 0 },     # vertex 4\n  { -a, a },    # vertex 5\n  { 0, a },     # vertex 6\n  { a*b, a*b }  # vertex 7\n}\n\nelements =\n{\n  { 0, 1, 4, 3, 0 },  # quad 0\n  { 3, 4, 7, 0 },     # tri 1\n  { 3, 7, 6, 0 },     # tri 2\n  { 2, 3, 6, 5, 0 }   # quad 3\n}\n\nboundaries =\n{\n  { 0, 1, 1 },\n  { 1, 4, 2 },\n  { 3, 0, 4 },\n  { 4, 7, 2 },\n  { 7, 6, 2 },\n  { 2, 3, 4 },\n  { 6, 5, 2 },\n  { 5, 2, 3 }\n}\n\ncurves =\n{\n  { 4, 7, 45 },  # +45 degree circular arcs\n  { 7, 6, 45 }\n}\n");

  /*** READ DATA - GENERAL ***/

  // Read number of initial uniform mesh refinements.
  int init_ref_num;
  if(!Get(f, &init_ref_num)) error("Could not read number of initial mesh refinements.");
  info("init_ref_num: %d", init_ref_num);
  B.set_initial_mesh_refinement(init_ref_num);

  // Read initial polynomial degree of elements.
  int init_p;
  if(!Get(f, &init_p)) error("Could not read initial polynomial degree.");
  info("init_p: %d", init_p);
  B.set_initial_poly_degree(init_p);

  // Read matrix solver.
  char* str = new char[255];
  if(!Get(f, str)) error("Could not read matrix solver name.");
  info("matrix solver: %s", str);
  B.set_matrix_solver(std::string(str));
  delete str;

  /*** READ DATA - MODEL-SPECIFIC ***/

  // Read list of material markers.
  int n_mat_markers;
  if(!Get(f, &n_mat_markers)) error("Could not read number of material markers.");
  info("n_mat_markers: %d", n_mat_markers);
  if(n_mat_markers <= 0) error("At least one material marker must be given.");
  std::vector<int> mat_markers;
  for (int i = 0; i < n_mat_markers; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a material marker.");
    info("mat_marker[%d]: %d", i, tmp);
    mat_markers.push_back(tmp);
  }
  B.set_material_markers(mat_markers);

  // Read list of c1 constants.
  std::vector<double> c1_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c1 constant.");
    info("c1_array[%d]: %g", i, tmp);
    c1_array.push_back(tmp);
  }
  B.set_c1_array(c1_array);

  // Read list of c2 constants.
  std::vector<double> c2_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c2 constant.");
    info("c2_array[%d]: %g", i, tmp);
    c2_array.push_back(tmp);
  }
  B.set_c2_array(c1_array);

  // Read list of c3 constants.
  std::vector<double> c3_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c3 constant.");
    info("c3_array[%d]: %g", i, tmp);
    c3_array.push_back(tmp);
  }
  B.set_c3_array(c1_array);

  // Read list of c4 constants.
  std::vector<double> c4_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c4 constant.");
    info("c4_array[%d]: %g", i, tmp);
    c4_array.push_back(tmp);
  }
  B.set_c4_array(c1_array);

  // Read list of c5 constants.
  std::vector<double> c5_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c5 constant.");
    info("c5_array[%d]: %g", i, tmp);
    c5_array.push_back(tmp);
  }
  B.set_c5_array(c1_array);

  // Read Dirichlet boundary markers.
  int n_bc_dirichlet;
  if(!Get(f, &n_bc_dirichlet)) error("Could not read number of Dirichlet boundary markers.");
  info("n_bc_dirichlet: %d", n_bc_dirichlet);
  std::vector<int> bdy_markers_dirichlet;
  for (int i = 0; i < n_bc_dirichlet; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a VALUE boundary marker.");
    info("bdy_markers_dirichlet[%d]: %d", i, tmp);
    bdy_markers_dirichlet.push_back(tmp);
  }
  if (n_bc_dirichlet > 0) B.set_dirichlet_markers(bdy_markers_dirichlet);

  // Read Dirichlet boundary values.
  std::vector<double> bdy_values_dirichlet;
  for (int i = 0; i < n_bc_dirichlet; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a Dirichlet boundary value.");
    info("bdy_values_dirichlet[%d]: %g", i, tmp);
    bdy_values_dirichlet.push_back(tmp);
  }
  if (n_bc_dirichlet > 0) B.set_dirichlet_values(bdy_values_dirichlet);

  // Read Neumann boundary markers.
  int n_bc_neumann;
  if(!Get(f, &n_bc_neumann)) error("Could not read number of Neumann boundary markers.");
  info("n_bc_neumann: %d", n_bc_neumann);
  std::vector<int> bdy_markers_neumann;
  for (int i = 0; i < n_bc_neumann; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a Neumann boundary marker.");
    info("bdy_markers_neumann[%d]: %d", i, tmp);
    bdy_markers_neumann.push_back(tmp);
  }
  if (n_bc_neumann > 0) B.set_neumann_markers(bdy_markers_neumann);

  // Read list of Neumann boundary values.
  std::vector<double> bdy_values_neumann;
  for (int i = 0; i < n_bc_neumann; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a Neumann boundary value.");
    info("bdy_values_neumann[%d]: %g", i, tmp);
    bdy_values_neumann.push_back(tmp);
  }
  if (n_bc_neumann > 0) B.set_neumann_values(bdy_values_neumann);

  // Read Newton boundary markers.
  int n_bc_newton;
  if(!Get(f, &n_bc_newton)) error("Could not read number of Newton boundary markers.");
  info("n_bc_newton: %d", n_bc_newton);
  std::vector<int> bdy_markers_newton;
  for (int i = 0; i < n_bc_newton; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a Newton boundary marker.");
    info("bdy_markers_newton[%d]: %d", i, tmp);
    bdy_markers_newton.push_back(tmp);
  }
  if (n_bc_newton > 0) B.set_newton_markers(bdy_markers_newton);

  // Read list of Newton boundary value pairs.
  std::vector<double_pair> bdy_values_newton;
  for (int i = 0; i < n_bc_newton; i++) {
    double tmp1, tmp2;
    if(!Get(f, &tmp1)) error("Could not read a Newton boundary value (first in pair).");
    if(!Get(f, &tmp2)) error("Could not read a Newton boundary value (second in pair).");
    info("bdy_values_newton[%d]: %g %g", i, tmp1, tmp2);
    bdy_values_newton.push_back(std::make_pair(tmp1, tmp2));
  }
  if (n_bc_newton > 0) B.set_newton_values(bdy_values_newton);

  /*** SOLVE THE PROBLEM ***/

  // Variables for time measurement.
  double assembly_time, solver_time;

  // Solve.
  bool success = B.calculate(assembly_time, solver_time);
  info("=====================");
  info("Computation finished.");
  info("assembly time: %g s", assembly_time);
  info("solver time: %g s", solver_time);
  if (!success) error("Computation failed.");

  /*** SHOW THE SOLUTION, GRADIENT, and SPACE ***/

  // Solution.
  Solution* sln = B.get_solution();
  ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
  view.show_mesh(false);
  view.show(sln);

  // Gradient magnitude.
  ScalarView gradview("Gradient", new WinGeom(445, 0, 440, 350));
  MagFilter grad(Hermes::Tuple<MeshFunction *>(sln, sln), Hermes::Tuple<int>(H2D_FN_DX, H2D_FN_DY));
  gradview.show_mesh(false);
  gradview.show(&grad);

  // Space.
  Space* space = B.get_space();
  OrderView oview("Mesh", new WinGeom(890, 0, 440, 350));
  oview.show(space);

  // Wait for the views to be closed.
  View::wait();

  return 1;
}
