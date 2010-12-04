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

  // Read number of initial uniform mesh refinements.
  int init_ref_num;
  if(!Get(f, &init_ref_num)) error("Could not read number of initial mesh refinements.");

  // Read initial polynomial degree of elements.
  int init_p;
  if(!Get(f, &init_p)) error("Could not read number of initial polynomial degree.");

  // Read number of material markers.
  int n_mat_markers;
  if(!Get(f, &n_mat_markers)) error("Could not read number of material markers.");
  if(n_mat_markers <= 0) error("At least one material marker must be given.");

  // Read list of material markers.
  std::vector<int> mat_markers;
  for (int i = 0; i < n_mat_markers; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a material marker.");
    mat_markers.push_back(tmp);
  }

  // Read list of c1 constants.
  std::vector<double> c1_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c1 constant.");
    c1_array.push_back(tmp);
  }

  // Read list of c2 constants.
  std::vector<double> c2_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c2 constant.");
    c2_array.push_back(tmp);
  }

  // Read list of c3 constants.
  std::vector<double> c3_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c3 constant.");
    c3_array.push_back(tmp);
  }

  // Read list of c4 constants.
  std::vector<double> c4_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c4 constant.");
    c4_array.push_back(tmp);
  }

  // Read list of c5 constants.
  std::vector<double> c5_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a c5 constant.");
    c5_array.push_back(tmp);
  }

  // Read number of Dirichlet boundary markers.
  int n_bc_dirichlet;
  if(!Get(f, &n_bc_dirichlet)) error("Could not read number of Dirichlet boundary markers.");

  // Read list of Dirichlet boundary markers.
  std::vector<int> bdy_markers_dirichlet;
  for (int i = 0; i < n_bc_dirichlet; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a VALUE boundary marker.");
    bdy_markers_dirichlet.push_back(tmp);
  }

  // Read list of Dirichlet boundary values.
  std::vector<double> bdy_values_dirichlet;
  for (int i = 0; i < n_bc_dirichlet; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a Dirichlet boundary value.");
    bdy_values_dirichlet.push_back(tmp);
  }

  // Read number of Neumann boundary markers.
  int n_bc_neumann;
  if(!Get(f, &n_bc_neumann)) error("Could not read number of Neumann boundary markers.");

  // Read list of Neumann boundary markers.
  std::vector<int> bdy_markers_neumann;
  for (int i = 0; i < n_bc_neumann; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a Neumann boundary marker.");
    bdy_markers_neumann.push_back(tmp);
  }

  // Read list of Neumann boundary values.
  std::vector<double> bdy_values_neumann;
  for (int i = 0; i < n_bc_neumann; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a Neumann boundary value.");
    bdy_values_neumann.push_back(tmp);
  }

  // Read number of Newton boundary markers.
  int n_bc_newton;
  if(!Get(f, &n_bc_newton)) error("Could not read number of Newton boundary markers.");

  // Read list of Newton boundary markers.
  std::vector<int> bdy_markers_newton;
  for (int i = 0; i < n_bc_newton; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a Newton boundary marker.");
    bdy_markers_newton.push_back(tmp);
  }

  // Read list of Newton boundary value pairs.
  std::vector<double_pair> bdy_values_newton;
  for (int i = 0; i < n_bc_newton; i++) {
    double tmp1, tmp2;
    if(!Get(f, &tmp1)) error("Could not read a Newton boundary value (first in pair).");
    if(!Get(f, &tmp2)) error("Could not read a Newton boundary value (second in pair).");
    bdy_values_newton.push_back(std::make_pair(tmp1, tmp2));
  }

  /*** FEED THE DATA INTO THE ELECTROSTATICS MODULE ***/

  // Initialize the Electrostatics class.
  Basic B;

  // Set mesh filename.
  B.set_mesh_str("\na = 1.0  # size of the mesh\nb = sqrt(2)/2\n\nvertices =\n{\n  { 0, -a },    # vertex 0\n  { a, -a },    # vertex 1\n  { -a, 0 },    # vertex 2\n  { 0, 0 },     # vertex 3\n  { a, 0 },     # vertex 4\n  { -a, a },    # vertex 5\n  { 0, a },     # vertex 6\n  { a*b, a*b }  # vertex 7\n}\n\nelements =\n{\n  { 0, 1, 4, 3, 0 },  # quad 0\n  { 3, 4, 7, 0 },     # tri 1\n  { 3, 7, 6, 0 },     # tri 2\n  { 2, 3, 6, 5, 0 }   # quad 3\n}\n\nboundaries =\n{\n  { 0, 1, 1 },\n  { 1, 4, 2 },\n  { 3, 0, 4 },\n  { 4, 7, 2 },\n  { 7, 6, 2 },\n  { 2, 3, 4 },\n  { 6, 5, 2 },\n  { 5, 2, 3 }\n}\n\ncurves =\n{\n  { 4, 7, 45 },  # +45 degree circular arcs\n  { 7, 6, 45 }\n}\n");
  
  // Set number of initial uniform mesh refinements.
  B.set_initial_mesh_refinement(init_ref_num);

  // Set initial polynomial degree of elements.
  B.set_initial_poly_degree(init_p);

  // Set material markers.
  B.set_material_markers(mat_markers);

  // Set c1 array.
  B.set_c1_array(c1_array);

  // Set c2 array.
  B.set_c2_array(c2_array);

  // Set c3 array.
  B.set_c3_array(c3_array);

  // Set c4 array.
  B.set_c4_array(c4_array);

  // Set c5 array.
  B.set_c5_array(c5_array);

  if (n_bc_dirichlet > 0) {
    // Set Dirichlet boundary markers.
    B.set_dirichlet_markers(bdy_markers_dirichlet);

    // Set Dirichlet boundary values.
    B.set_dirichlet_values(bdy_values_dirichlet);
  }

  if (n_bc_neumann > 0) {
    // Set Neumann boundary markers.
    B.set_neumann_markers(bdy_markers_neumann);

    // Set Neumann boundary values.
    B.set_neumann_values(bdy_values_neumann);
  }

  if (n_bc_newton > 0) {
    // Set Newton boundary markers.
    B.set_newton_markers(bdy_markers_newton);

    // Set Newton boundary values.
    B.set_newton_values(bdy_values_newton);
  }

  /*** SOLVE THE PROBLEM ***/

  Solution phi;
  bool success = B.calculate(&phi);
  if (!success) error("Computation failed.");

  /*** SHOW THE SOLUTION ***/

  ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
  view.show(&phi);

  // Compute and show gradient magnitude.
  // (Note that the gradient at the re-entrant
  // corner needs to be truncated for visualization purposes.)
  ScalarView gradview("Gradient", new WinGeom(450, 0, 440, 350));
  MagFilter grad(Hermes::Tuple<MeshFunction *>(&phi, &phi), Hermes::Tuple<int>(H2D_FN_DX, H2D_FN_DY));
  gradview.show(&grad);

  // Wait for the views to be closed.
  View::wait();

  return 1;
}
