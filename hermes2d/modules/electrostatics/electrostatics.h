#include "hermes2d.h"

// This is a simple electrostatics module based on the Hermes library. The purpose
// of the module is to provide higher-level API where the user is not exposed to 
// wesk forms and other technical details. 
//
// PDE: Poisson equation -div(epsilon grad phi) = rho 
//      epsilon... permittivity (given, element-wise constant)
//      phi...     electric potential (unknown)
//      rho...     charge density (element-wise constant) 
//
// BC: piecewise-constant on Dirichlet boundaries
//     piecewise constant normal derivative on Neumann boundaries

// This class allows the user to define the mesh, boundary conditions,
// material constants and right-hand side via methods starting with "set".
// After that, the problem is solved by calling the method calculate().
class Electrostatics {
public:
  Electrostatics();
  ~Electrostatics();
  // Set mesh as a string (see example at the end of this file).
  void set_mesh_str(const std::string &mesh);
  // Set number of initial unifrm mesh refinements.
  void set_initial_mesh_refinement(int init_ref_num);
  // Set initial polynomial degree in all elements.
  void set_initial_poly_degree(int p);
  // Set material markers. These markers define subdomains.
  void set_material_markers(const std::vector<int> &mat_markers);
  // Set constant permittivities for all material subdomains.
  void set_permittivity_array(const std::vector<double> &p_array);
  // Set constant charge densities for all material subdomains.
  void set_charge_density_array(const std::vector<double> &cd_array);
  // Set an array of boundary markers for boundary conditions with 
  // prescribed electric potential (Dirichlet conditions). 
  void set_boundary_markers_value(const std::vector<int> &bdy_markers_val);
  // Set an array of constant values for Dirichlet boundary markers.
  void set_boundary_values(const std::vector<double> &bc_val);
  // Set an array of boundary markers for boundary conditions with 
  // prescribed surface charge density (normal derivative of electric potential).
  // These are often called Neumann boundary conditions. 
  void set_boundary_markers_derivative(const std::vector<int> &bdy_markers_der);
  // Set an array of normal derivatives Neumann boundary markers.
  void set_boundary_derivatives(const std::vector<double> &bc_der);
  // SOlve the problem and return the solution.
  bool calculate(Solution* phi);

  BCTypes bc_types;                         // Associates BC markers with BC boundary types

private:
  std::string mesh_str;
  int init_ref_num;
  int init_p;
  std::vector<int> mat_markers;             // Array of material markers (>= 0).
                                            // Example: [1, 3, 0] 
  std::vector<int> mat_permut;              // For a material marker, this array gives 
                                            // its index in the list of material constants.
                                            // Example (for the above array): [-1, 0, -1, 1, 2]
  std::vector<double> permittivity_array;   // Array of corresponding permittivities.
  std::vector<double> charge_density_array; // Array of corresponding charge densities.
  std::vector<double> bc_val;               // List of Dirichlet boundary values (in the same order)..
  std::vector<double> bc_der;               // List of Neumann boundary values (in the same order).
  std::vector<int> bc_permut;               // For any boundary marker (either Dirichlet or Neumann, 
                                            // this array gives its index in the corresponding list of
                                            // (either Dirichlet or Neumann) boundary conditions.
                                            // Example (for the above Dirichlet and Neumann arrays): 
                                            // [-1, 2, 0, 0, 1, 1]


  Mesh* mesh;
  H1Space* space;
};

/* Mesh string example. This is the mesh from the Hermes tutorial example 01-mesh,
see http://hpfem.org/hermes/doc/src/hermes2d/tutorial-1/mesh.html.

"\na = 1.0  # size of the mesh\nb = sqrt(2)/2\n\nvertices =\n{\n  { 0, -a },    # vertex 0\n  { a, -a },    # vertex 1\n  { -a, 0 },    # vertex 2\n  { 0, 0 },     # vertex 3\n  { a, 0 },     # vertex 4\n  { -a, a },    # vertex 5\n  { 0, a },     # vertex 6\n  { a*b, a*b }  # vertex 7\n}\n\nelements =\n{\n  { 0, 1, 4, 3, 0 },  # quad 0\n  { 3, 4, 7, 0 },     # tri 1\n  { 3, 7, 6, 0 },     # tri 2\n  { 2, 3, 6, 5, 0 }   # quad 3\n}\n\nboundaries =\n{\n  { 0, 1, 1 },\n  { 1, 4, 2 },\n  { 3, 0, 4 },\n  { 4, 7, 2 },\n  { 7, 6, 2 },\n  { 2, 3, 4 },\n  { 6, 5, 2 },\n  { 5, 2, 3 }\n}\n\ncurves =\n{\n  { 4, 7, 45 },  # +45 degree circular arcs\n  { 7, 6, 45 }\n}\n"

 */
