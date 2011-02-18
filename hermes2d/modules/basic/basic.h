#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"

// The module is a library.
#ifndef EXPORT_HERMES_MODULE
#define EXPORT_HERMES_MODULE
#endif

#include "hermes2d.h"

// This is a simple generic module for a linear second-order PDE based on the Hermes 
// library. The purpose of modules like this is to provide a higher-level API where 
// the user is not exposed to wesk forms, spaces, solver initialization, and other 
// technical details. 
//
// PDE: -div(c1 \nabla u) + (c2, c3) \cdot \nabla u + c4 u = c5 
//      c1 ... equation parameter, element-wise constant
//      c2 ... equation parameter, element-wise constant
//      c3 ... equation parameter, element-wise constant
//      c4 ... equation parameter, element-wise constant
//      c2 ... equation parameter, element-wise constant
//
// Possible BC: Dirichlet with piecewise-constant values (u = const)
//              Neumann with piecewise-constant normal derivatives (du/dn = const)
//              Newton with piecewise-constant parameters (const_1 u + du/dn = const_2)
//
// This class allows the user to define the mesh, boundary conditions,
// material constants and right-hand side via methods starting with "set".
// After that, the problem is solved by calling the method calculate().

typedef std::pair<double, double> double_pair;

class HERMES_MODULE_API ModuleBasic {
public:
  ModuleBasic();

  ~ModuleBasic();

  // Set mesh as a string (see example at the end of this file).
  void set_mesh_str(const std::string &mesh);

  // Set number of initial unifrm mesh refinements.
  void set_initial_mesh_refinement(int init_ref_num);

  // Set initial polynomial degree in all elements.
  void set_initial_poly_degree(int p);

  // Set material markers. These markers define subdomains.
  void set_material_markers(const std::vector<int> &mat_markers);

  // Set parameter c1 for all material subdomains.
  void set_c1_array(const std::vector<double> &c1_array);

  // Set parameter c2 for all material subdomains.
  void set_c2_array(const std::vector<double> &c2_array);

  // Set parameter c3 for all material subdomains.
  void set_c3_array(const std::vector<double> &c3_array);

  // Set parameter c4 for all material subdomains.
  void set_c4_array(const std::vector<double> &c4_array);

  // Set parameter c5 for all material subdomains.
  void set_c5_array(const std::vector<double> &c5_array);

  // Set Dirichlet boundary markers.
  void set_dirichlet_markers(const std::vector<int> &bdy_markers_dirichlet);

  // Set Dirichlet boundary values.
  void set_dirichlet_values(const std::vector<int> &bdy_markers_dirichlet,
                            const std::vector<double> &bdy_values_dirichlet);

  // Set Neumann boundary markers.
  void set_neumann_markers(const std::vector<int> &bdy_markers_neumann);

  // Set Neumann boundary values.
  void set_neumann_values(const std::vector<double> &bdy_values_neumann);

  // Set Newton boundary markers.
  void set_newton_markers(const std::vector<int> &bdy_markers_newton);

  // Set Newton boundary value pairs.
  void set_newton_values(const std::vector<double_pair> &bdy_values_newton);

  // Sanity check of material markers and material constants.
  void materials_sanity_check();

  // Setting mesh.
  void set_mesh(Mesh* m);

  // Get mesh.
  Mesh* get_mesh();

  // Solve the problem and return the solution.
  bool calculate();

  // This class associates BC markers with BC boundary types.
  BCTypes bc_types;

  // This class provides values for Dirichlet (essential) BC markers.
  BCValues bc_values;

  // Get mesh string.
  const char* get_mesh_string();

  // Clear mesh string.
  void clear_mesh_string();

  // Get solution.
  void get_solution(Solution *s);

  // Get space.
  void get_space(H1Space* s);

  // Set matrix solver.
  void set_matrix_solver(std::string solver_name);

  // Get matrix solver.
  MatrixSolverType get_matrix_solver();

  // Perform basic sanity checks, create mesh, perform 
  // uniform refinements, create space, register weak forms.
  void create_mesh_space_forms();

  // Get assembly time.
  double get_assembly_time();

  // Get solver time.
  double get_solver_time();


protected:
  std::string mesh_str;
  int init_ref_num;
  int init_p;
  std::vector<int> mat_markers;                // Array of material markers (>= 0).
                                               // Example: [1, 3, 0] 
  std::vector<int> mat_permut;                 // For a material marker, this array gives 
                                               // its index in the list of material constants.
                                               // Example (for the above array): [-1, 0, -1, 1, 2]
  std::vector<double> c1_array;                // Equation parameter c1 (array for all material subdomains).
  std::vector<double> c2_array;                // Equation parameter c2 (array for all material subdomains).
  std::vector<double> c3_array;                // Equation parameter c3 (array for all material subdomains).
  std::vector<double> c4_array;                // Equation parameter c4 (array for all material subdomains).
  std::vector<double> c5_array;                // Equation parameter c5 (array for all material subdomains).
  std::vector<int> bdy_markers_dirichlet;      // List of Dirichlet boundary markers.
  std::vector<double> bdy_values_dirichlet;    // List of Dirichlet boundary values.
  std::vector<int> bdy_markers_neumann;        // List of Neumann boundary markers.
  std::vector<double> bdy_values_neumann;      // List of Neumann boundary values.
  std::vector<int> bdy_markers_newton;         // List of Newton boundary markers.
  std::vector<double_pair> bdy_values_newton;  // List of Newton boundary value pairs.
  std::vector<int> bc_permut;                  // Boundary conditions are entered in the same order as 
                                               // boundary markers. The markers need to be integers greater
                                               // than zero, but otherwise they are the user's choice. 
                                               // Therefore we introduce a bc_permut array that for any 
                                               // boundary marker gives its index in the list of boundary 
                                               // conditions. 
  MatrixSolverType matrix_solver;              // Possibilities: SOLVER_AMESOS, SOLVER_ZATECOO, SOLVER_MUMPS, 
                                               // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

  // Finite element mesh.
  Mesh* mesh;

  // Finite element space;
  H1Space* space;

  // Weak form.
  WeakForm* wf;

  // Solution.
  Solution* sln;

  // Assembly and solver times.
  double assembly_time, solver_time;
};

/* Mesh string example. This is the mesh from the Hermes tutorial example 01-mesh,
see http://hpfem.org/hermes/doc/src/hermes2d/tutorial-1/mesh.html.

"\na = 1.0  # size of the mesh\nb = sqrt(2)/2\n\nvertices =\n{\n  { 0, -a },    # vertex 0\n  { a, -a },    # vertex 1\n  { -a, 0 },    # vertex 2\n  { 0, 0 },     # vertex 3\n  { a, 0 },     # vertex 4\n  { -a, a },    # vertex 5\n  { 0, a },     # vertex 6\n  { a*b, a*b }  # vertex 7\n}\n\nelements =\n{\n  { 0, 1, 4, 3, 0 },  # quad 0\n  { 3, 4, 7, 0 },     # tri 1\n  { 3, 7, 6, 0 },     # tri 2\n  { 2, 3, 6, 5, 0 }   # quad 3\n}\n\nboundaries =\n{\n  { 0, 1, 1 },\n  { 1, 4, 2 },\n  { 3, 0, 4 },\n  { 4, 7, 2 },\n  { 7, 6, 2 },\n  { 2, 3, 4 },\n  { 6, 5, 2 },\n  { 5, 2, 3 }\n}\n\ncurves =\n{\n  { 4, 7, 45 },  # +45 degree circular arcs\n  { 7, 6, 45 }\n}\n"

 */
