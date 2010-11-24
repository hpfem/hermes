#include "hermes2d.h"
#include "disc.h"
#include "electrostatics.h"

int main(int argc, char* argv[])
{

  /*** RECEIVE DATA ***/

  // Check the number of command line arguments.
  if(argc != 2) error("Configuration file missing.");
  // Open configuration file.
  FILE* f = fopen(argv[1], "r");
  if(f == NULL) error("Cannot open file %s.", argv[1]);
  // Number of initial uniform mesh refinements.
  int init_ref_num;
  if(!Get(f, &init_ref_num)) error("Could not read number of initial mesh refinements.");
  // Initial polynomial degree of elements.
  int init_p;
  if(!Get(f, &init_p)) error("Could not read number of initial polynomial degree.");
  // Number of material markers.
  int n_mat_markers;
  if(!Get(f, &n_mat_markers)) error("Could not read number of material markers.");
  if(n_mat_markers <= 0) error("At least one material merker must be given.");
  // Array of material markers.
  std::vector<int> mat_markers;
  for (int i = 0; i < n_mat_markers; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a material marker.");
    mat_markers.push_back(tmp);
  }
  // Array of material constants (permittivities).
  std::vector<double> permittivity_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a permittivity.");
    permittivity_array.push_back(tmp);
  }
  // Array of material constants (charge densities).
  std::vector<double> charge_density_array;
  for (int i = 0; i < n_mat_markers; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a charge density.");
    charge_density_array.push_back(tmp);
  }
  // Number of VALUE boundary conditions.
  int n_bc_value;
  if(!Get(f, &n_bc_value)) error("Could not read number of VALUE boundary markers.");
  if(n_bc_value <= 0) error("At least one VALUE boundary marker must be given.");
  // VALUE boundary markers.
  std::vector<int> bc_markers_value;
  for (int i = 0; i < n_bc_value; i++) {
    int tmp;
    if(!Get(f, &tmp)) error("Could not read a VALUE boundary marker.");
    bc_markers_value.push_back(tmp);
  }
  // Boundary values (electric potential).
  std::vector<double> bc_values;
  for (int i = 0; i < n_bc_value; i++) {
    double tmp;
    if(!Get(f, &tmp)) error("Could not read a boundary value.");
    bc_values.push_back(tmp);
  }
  // Number of DERIVATIVE boundary conditions.
  int n_bc_derivative;
  if(!Get(f, &n_bc_derivative)) error("Could not read number of DERIVATIVE boundary markers.");
  // DERIVATIVE boundary markers.
  std::vector<int> bc_markers_derivative;
  std::vector<double> bc_derivatives;
  // Boundary derivatives.
  if(n_bc_derivative >= 1) {
    for (int i = 0; i < n_bc_derivative; i++) {
      int tmp;
      if(!Get(f, &tmp)) error("Could not read a DERIVATIVE boundary marker.");
      bc_markers_derivative.push_back(tmp);
    }
    for (int i = 0; i < n_bc_derivative; i++) {
      double tmp;
      if(!Get(f, &tmp)) error("Could not read a boundary derivative.");
      bc_derivatives.push_back(tmp);
    }
  }


  /*** FEED THE DATA INTO THE ELECTROSTATICS MODULE ***/

  // Initialize the Electrostatics class.
  Electrostatics E;

  // Set mesh filename.
  E.set_mesh_str("\na = 1.0  # size of the mesh\nb = sqrt(2)/2\n\nvertices =\n{\n  { 0, -a },    # vertex 0\n  { a, -a },    # vertex 1\n  { -a, 0 },    # vertex 2\n  { 0, 0 },     # vertex 3\n  { a, 0 },     # vertex 4\n  { -a, a },    # vertex 5\n  { 0, a },     # vertex 6\n  { a*b, a*b }  # vertex 7\n}\n\nelements =\n{\n  { 0, 1, 4, 3, 0 },  # quad 0\n  { 3, 4, 7, 0 },     # tri 1\n  { 3, 7, 6, 0 },     # tri 2\n  { 2, 3, 6, 5, 0 }   # quad 3\n}\n\nboundaries =\n{\n  { 0, 1, 1 },\n  { 1, 4, 2 },\n  { 3, 0, 4 },\n  { 4, 7, 2 },\n  { 7, 6, 2 },\n  { 2, 3, 4 },\n  { 6, 5, 2 },\n  { 5, 2, 3 }\n}\n\ncurves =\n{\n  { 4, 7, 45 },  # +45 degree circular arcs\n  { 7, 6, 45 }\n}\n");
  
  // Set number of initial uniform mesh refinements.
  E.set_initial_mesh_refinement(init_ref_num);

  // Set initial polynomial degree of elements.
  E.set_initial_poly_degree(init_p);

  // Set material markers (also checks compatibility with the mesh file).
  E.set_material_markers(mat_markers);

  // Set permittivity array.
  E.set_permittivity_array(permittivity_array);

  // Set charge density array.
  E.set_charge_density_array(charge_density_array);

  // Set VALUE boundary markers (also check with the mesh file).
  E.set_boundary_markers_value(bc_markers_value);

  // Set boundary values.
  E.set_boundary_values(bc_values);

  // Set DERIVATIVE boundary markers (also check with the mesh file).
  E.set_boundary_markers_derivative(bc_markers_derivative);

  // Set boundary derivatives.
  E.set_boundary_derivatives(bc_derivatives);


  /*** SOLVE THE PROBLEM ***/

  Solution phi;
  bool success = E.calculate(&phi);
  if (!success) error("Computation failed.");

  /*** SHOW THE SOLUTION ***/

  ScalarView view("Solution", new WinGeom(0, 0, 440, 350));
  view.show(&phi);

  // Compute and show gradient magnitude.
  // (Note that the gradient at the re-entrant
  // corner needs to be truncated for visualization purposes.)
  ScalarView gradview("Gradient", new WinGeom(450, 0, 400, 350));
  MagFilter grad(Tuple<MeshFunction *>(&phi, &phi), Tuple<int>(H2D_FN_DX, H2D_FN_DY));
  gradview.show(&grad);

  // Wait for the views to be closed.
  View::wait();

  return 1;
}
