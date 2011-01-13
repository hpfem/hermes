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

#include "read_model_data.cpp"

  /*** SOLVE THE PROBLEM ***/

  // Variables for time measurement.
  double assembly_time, solver_time;

  // Solve.
  bool success = B.calculate();
  info("=====================");
  info("Computation finished.");
  info("assembly time: %g s", B.get_assembly_time());
  info("solver time: %g s", B.get_solver_time());
  if (!success) error("Computation failed.");

  /*** SHOW THE SOLUTION, GRADIENT, and SPACE ***/

  // Solution.
  Solution sln;
  B.get_solution(&sln);
  ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
  view.show_mesh(false);
  view.show(&sln);

  // Gradient magnitude.
  ScalarView gradview("Gradient", new WinGeom(445, 0, 440, 350));
  MagFilter grad(Hermes::vector<MeshFunction *>(&sln, &sln), Hermes::vector<int>(H2D_FN_DX, H2D_FN_DY));
  gradview.show_mesh(false);
  gradview.show(&grad);

  // Space.
  BCTypes bctypes;
  BCValues bcvalues;
  Mesh m;
  H1Space space(B.get_mesh(), &bctypes, &bcvalues, 1); // FIXME: this is a hack since constructor 
                                                       // to Space needs some mesh. The mesh is not used.
  B.get_space(&space);
  OrderView oview("Mesh", new WinGeom(890, 0, 440, 350));
  oview.show(&space);

  // Wait for the views to be closed.
  View::wait();

  return 1;
}
