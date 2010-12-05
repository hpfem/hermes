#include "hermes2d.h"
#include "numerical_flux.h"

// This is a module for Euler equations. 
// The purpose of modules like this is to provide a higher-level API where 
// the user is not exposed to wesk forms, spaces, solver initialization, and other 
// technical details. 
//
// PDE: \partial w / \partial t + \partial f_x / \partial x + \partial f_y / \partial y.
//      w = (w0, w1, w2, w3) = (rho, rho * v_x, rho * v_y, E).
//      f_x (w0, w1, w2, w3) = (w1, w1^2/w0 + p, w1*w2/w0, w1/w0*(w3 + p)).
//      f_y (w0, w1, w2, w3) = (w2, w2^2/w0 + p, w1*w2/w0, w2/w0*(w3 + p)).
//      p (w0, w1, w2, w3) = (\kappa - 1) * (w3 - (w1^2 + w2^2)/ (2 * w0)).
//
// Possible BC: "Solid wall" boundary condition on the bottom and top part of the domain.
//              Prescribed state (w0, w1 w2, w3) on the inlet part of the domain.
//              Prescribed state (w0, w1, w2, w3) on the outlet.
//              
// This class allows the user to define the mesh, boundary condition values
// via methods starting with "set".
// After that, the problem is solved by calling the method calculate().

class HERMES_API CompressibleEuler {
public:
  CompressibleEuler(bool use_vector_valued_forms = true, bool calculate_time_derivative_approximation = true);

  ~CompressibleEuler();

  // Set mesh as a string (see example at the end of this file).
  void set_mesh_str(const std::string &mesh);

  // Set mesh filename.
  void set_mesh_filename(const std::string &mesh);

  // Set number of initial unifrm mesh refinements.
  void set_initial_mesh_refinement(int init_ref_num);

  // Set initial polynomial degree in all elements.
  void set_initial_poly_degree(Ord2 p);
  void set_initial_poly_degree(int p);

  // Set the value used in the CFL condition.
  void set_CFL_value(double CFL);

  // Set the initial time step (checked vs. the CFL condition).
  void set_time_step(double tau);

  // Set the length of the time interval.
  void set_time_interval_length(double T);

  // Set the state.
  void set_state(double p_ext, double rho_ext, double v1_ext, double v2_ext);

  // Set the specific heat ratio.
  void set_kappa(double kappa);

  // Set the boundary markers.
  void set_inlet_outlet_marker(int marker);
  void set_solid_wall_marker(int marker);

  // Solve the problem and return the solution.
  bool calculate(Hermes::Tuple<Solution*> sln);

  // This class associates BC markers with BC boundary types.
  // Only default in the current implementation.
  BCTypes* bc_types;

  // Use of filters for visualization.
  // If called with set == true, filters are used and velocity field and pressure
  // are displayed, otherwise the solution components are.
  // By default filtres are not used.
  void use_filters_for_visual(bool set = false);

  // Saving of screenshots.
  // saving_frequency is 0 by default (nothing is saved).
  void save_every_th_solution_visualization(unsigned int frequency); 

  // Set the usage of vector valued forms.
  // This is enabled by default (in the constructor).
  void set_use_vector_valued_forms(bool set = true);

  // Set the calculation of the time derivative approximation.
  // This is enabled by default (in the constructor).
  void set_calculate_time_derivative_approximation(bool set = true);

private:
  // At the beginning of calculate() checks that everything is set.
  bool all_set();

  // Calculation enhancements.
  bool use_vector_valued_forms;
  bool calculate_time_derivative_approximation;

  // Use of filters. Set to true by the method set_filters_for_visual().
  bool filters_for_visual;

  // If the user wishes to save screenshots of the solution.
  int saving_frequency;

  std::string mesh_str;
  std::string mesh_filename;
  int init_ref_num;
  Ord2 init_p;

  // Problem parameters.
  double CFL;
  double tau, T;
  double p_ext, rho_ext, v1_ext, v2_ext;

  NumericalFlux* num_flux;

  // Boundary value markers.
  int inlet_outlet_marker;
  int solid_wall_marker;

  // Finite element mesh.
  // TODO: change this when multimesh assembling is done.
  Mesh* mesh;

  // Finite element space;
  Hermes::Tuple<Space*>* spaces;

};

/* Mesh string example. This is the mesh from the Hermes tutorial example 01-mesh,
see http://hpfem.org/hermes/doc/src/hermes2d/tutorial-1/mesh.html.

"\na = 1.0  # size of the mesh\nb = sqrt(2)/2\n\nvertices =\n{\n  { 0, -a },    # vertex 0\n  { a, -a },    # vertex 1\n  { -a, 0 },    # vertex 2\n  { 0, 0 },     # vertex 3\n  { a, 0 },     # vertex 4\n  { -a, a },    # vertex 5\n  { 0, a },     # vertex 6\n  { a*b, a*b }  # vertex 7\n}\n\nelements =\n{\n  { 0, 1, 4, 3, 0 },  # quad 0\n  { 3, 4, 7, 0 },     # tri 1\n  { 3, 7, 6, 0 },     # tri 2\n  { 2, 3, 6, 5, 0 }   # quad 3\n}\n\nboundaries =\n{\n  { 0, 1, 1 },\n  { 1, 4, 2 },\n  { 3, 0, 4 },\n  { 4, 7, 2 },\n  { 7, 6, 2 },\n  { 2, 3, 4 },\n  { 6, 5, 2 },\n  { 5, 2, 3 }\n}\n\ncurves =\n{\n  { 4, 7, 45 },  # +45 degree circular arcs\n  { 7, 6, 45 }\n}\n"

 */
