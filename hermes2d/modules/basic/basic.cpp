#include "basic.h"

// FIXME: Because of callbacks that depend on global variables
// we need to define a few global arrays below, but this is 
// a temporary solution. Like this, one cannot have two instances 
// of the basic module at the same time -- their global variables 
// would interfere with each other.
Hermes::Tuple<int> _global_mat_markers;
std::vector<double> _global_c1_array;
std::vector<double> _global_c2_array;
std::vector<double> _global_c3_array;
std::vector<double> _global_c4_array;
std::vector<double> _global_c5_array;
std::vector<double> _global_bdy_values_dirichlet;
std::vector<double> _global_bdy_values_neumann;
std::vector<double_pair> _global_bdy_values_newton;
BCTypes* _global_bc_types = NULL;
void *_global_data;

// Matrix solver.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_ZATECOO, SOLVER_MUMPS, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Weak form (volumetric, left-hand side).
template<typename Real, typename Scalar>
Scalar bilinear_form_vol(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                         Geom<Real> *e, ExtData<Scalar> *ext)
{
  int elem_marker = e->elem_marker;
  double c1, c2, c3, c4;
  if (elem_marker < 0) {
    // This is for Order calculation only:
    c1 = c2 = c3 = c4 = -5555.0;
  } else {
    int index = _global_mat_markers.find_index(elem_marker);
    c1 = _global_c1_array[index];
    c2 = _global_c2_array[index];
    c3 = _global_c3_array[index];
    c4 = _global_c4_array[index];
  }
  return  
    c1 * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v)
    + c2 * int_dudx_v<Real, Scalar>(n, wt, u, v)
    + c3 * int_dudy_v<Real, Scalar>(n, wt, u, v)
    + c4 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Weak form (volumetric, right-hand side).
template<typename Real, typename Scalar>
Scalar linear_form_vol(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext)
{
  int elem_marker = e->elem_marker;
  double c5;
  if (elem_marker < 0) {
    // This is for Order calculation only:
    c5 = -5555.0;
  } else {
    int index = _global_mat_markers.find_index(elem_marker);
    c5 = _global_c5_array[index];
  }
  return c5 * int_v<Real, Scalar>(n, wt, v);
}

// Weak form (surface, left-hand side).
template<typename Real, typename Scalar>
Scalar bilinear_form_surf_newton(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                          Geom<Real> *e, ExtData<Scalar> *ext)
{
  int edge_marker = e->edge_marker;
  int elem_marker = e->elem_marker;
  double const_newton_1;
  double c1;
  if (edge_marker < 0) {
    // This is for Order calculation only:
    const_newton_1 = -5555.0;
    c1 = -5555.0;
  } else {
    if (_global_bdy_values_newton.size() > 0) {
      double_pair newton_pair = _global_bdy_values_newton[_global_bc_types->find_index_newton(edge_marker)];
      const_newton_1 = newton_pair.first;
    }
    else error("Internal in module Basic: bilinear_form_surf_newton() should not have been called.");
    c1 = _global_c1_array[_global_mat_markers.find_index(elem_marker)];
  }
  return c1 * const_newton_1 * int_u_v<Real, Scalar>(n, wt, u, v);
}

// Weak form (surface, neumann, right-hand side).
template<typename Real, typename Scalar>
Scalar linear_form_surf_neumann(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                        Geom<Real> *e, ExtData<Scalar> *ext)
{
  int edge_marker = e->edge_marker;
  int elem_marker = e->elem_marker;
  double const_neumann;
  double c1;
  double result = 0;
  if (edge_marker < 0) {
    // This is for Order calculation only:
    const_neumann = -5555.0;
    c1 = -5555.0;
  } else {
    if (_global_bdy_values_neumann.size() > 0) {
      int index = _global_bc_types->find_index_neumann(edge_marker);
      const_neumann = _global_bdy_values_neumann[index];
    }
    else error("Internal in module Basic: linear_form_surf_neumann() should not have been called.");
    c1 = _global_c1_array[_global_mat_markers.find_index(elem_marker)];  
  }
  return c1 * const_neumann * int_v<Real, Scalar>(n, wt, v);
}

// Weak form (surface, newton, right-hand side).
template<typename Real, typename Scalar>
Scalar linear_form_surf_newton(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                        Geom<Real> *e, ExtData<Scalar> *ext)
{
  int edge_marker = e->edge_marker;
  int elem_marker = e->elem_marker;
  double const_newton_2;
  double c1;
  double result = 0;
  if (edge_marker < 0) {
    // This is for Order calculation only:
    const_newton_2 = -5555.0;
    c1 = -5555.0;
  } else {
    if (_global_bdy_values_newton.size() > 0) {
      int index = _global_bc_types->find_index_newton(edge_marker);
      double_pair newton_pair = _global_bdy_values_newton[index];
      const_newton_2 = newton_pair.second;
    }
    else error("Internal in module Basic: linear_form_surf_newton() should not have been called.");
    c1 = _global_c1_array[_global_mat_markers.find_index(elem_marker)];  
  }
  return c1 * const_newton_2 * int_v<Real, Scalar>(n, wt, v);
}

// Look up an integer number in an array.
bool find_index(const std::vector<int> &array, int x, int &i_out)
{
  for (int i=0; i < array.size(); i++)
    if (array[i] == x) {
      i_out = i;
      return true;
    }
  return false;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  Basic *self = (Basic *) _global_data;
  int idx = self->bc_types.find_index_dirichlet(ess_bdy_marker);
  return _global_bdy_values_dirichlet[idx];
}

// Constructor.
Basic::Basic()
{
  init_ref_num = -1;
  init_p = -1;
  mesh = new Mesh();
  space = NULL;

  // FIXME: these global arrays need to be removed.
  _global_data = this;
  _global_bc_types = &(this->bc_types);
}

// Destructor.
Basic::~Basic()
{
  delete this->mesh;
  delete this->space;
}

// Set mesh via a string. See basic.h for an example of such a string.
void Basic::set_mesh_str(const std::string &mesh)
{
    this->mesh_str = mesh;
}

// Set the number of initial uniform mesh refinements.
void Basic::set_initial_mesh_refinement(int init_ref_num) 
{
  this->init_ref_num = init_ref_num;
}

// Set initial poly degree in elements.
void Basic::set_initial_poly_degree(int p) 
{
  this->init_p = p;
}

// Set material markers, and check compatibility with mesh file.
void Basic::set_material_markers(const std::vector<int> &m_markers)
{
  this->mat_markers = m_markers;
  _global_mat_markers = m_markers;
}

// Set c1 array.
void Basic::set_c1_array(const std::vector<double> &c1_array)
{
  int n = c2_array.size();
  for (int i = 0; i < n; i++) 
  if (c2_array[i] <= 1e-10) error("The c1 array needs to be positive.");
  this->c1_array = c1_array;
  _global_c1_array = c1_array;
}

// Set c2 array.
void Basic::set_c2_array(const std::vector<double> &c2_array)
{
  this->c2_array = c2_array;
  _global_c2_array = c2_array;
}

// Set c3 array.
void Basic::set_c3_array(const std::vector<double> &c3_array)
{
  this->c3_array = c3_array;
  _global_c3_array = c3_array;
}

// Set c4 array.
void Basic::set_c4_array(const std::vector<double> &c4_array)
{
  this->c4_array = c4_array;
  _global_c4_array = c4_array;
}

// Set c5 array.
void Basic::set_c5_array(const std::vector<double> &c5_array)
{
  this->c5_array = c5_array;
  _global_c5_array = c5_array;
}

// Set Dirichlet boundary markers.
void Basic::set_dirichlet_markers(const std::vector<int> &bdy_markers_dirichlet)
{
  this->bdy_markers_dirichlet = bdy_markers_dirichlet;
  Hermes::Tuple<int> t;
  t = bdy_markers_dirichlet;
  this->bc_types.add_bc_dirichlet(t);
}

// Set Dirichlet boundary values.
void Basic::set_dirichlet_values(const std::vector<double> &bdy_values_dirichlet)
{
  this->bdy_values_dirichlet = bdy_values_dirichlet;
  _global_bdy_values_dirichlet = bdy_values_dirichlet;
}

// Set Neumann boundary markers.
void Basic::set_neumann_markers(const std::vector<int> &bdy_markers_neumann)
{
  this->bdy_markers_neumann = bdy_markers_neumann;
  Hermes::Tuple<int> t;
  t = bdy_markers_neumann;
  this->bc_types.add_bc_neumann(t);
}

// Set Neumann boundary values.
void Basic::set_neumann_values(const std::vector<double> &bdy_values_neumann)
{
  this->bdy_values_neumann = bdy_values_neumann;
  _global_bdy_values_neumann = bdy_values_neumann;
}

// Set Newton boundary markers.
void Basic::set_newton_markers(const std::vector<int> &bdy_markers_newton)
{
  this->bdy_markers_newton = bdy_markers_newton;
  Hermes::Tuple<int> t;
  t = bdy_markers_newton;
  this->bc_types.add_bc_newton(t);
}

// Set Newton boundary values.
void Basic::set_newton_values(const std::vector<double_pair> &bdy_values_newton)
{
  this->bdy_values_newton = bdy_values_newton;
  _global_bdy_values_newton = bdy_values_newton;
}

// Solve the problem.
bool Basic::calculate(Solution* phi) 
{
  /* SANITY CHECKS */

  this->bc_types.check_consistency();

  // Sanity check of material markers and material constants.
  int n_mat_markers = this->mat_markers.size();
  if (n_mat_markers != this->c1_array.size()) error("Wrong length of c1 array.");
  if (n_mat_markers != this->c2_array.size()) error("Wrong length of c2 array.");
  if (n_mat_markers != this->c3_array.size()) error("Wrong length of c3 array.");
  if (n_mat_markers != this->c4_array.size()) error("Wrong length of c4 array.");
  if (n_mat_markers != this->c5_array.size()) error("Wrong length of c5 array.");

  // Making sure that material markers are nonnegative (>= 0).
  for (int i = 0; i < n_mat_markers; i++) {
    if(this->mat_markers[i] < 0) error("Material markers must be nonnegative.");
  }

  /* BEGIN THE COMPUTATION */

  // Load the mesh.
  H2DReader mloader;
  mloader.load_str(this->mesh_str.c_str(), this->mesh);

  // Debug.
  /*
  MeshView m("", 0, 0, 400, 400);
  m.show(this->mesh);
  View::wait();
  */

  // Perform initial uniform mesh refinements.
  for (int i = 0; i < this->init_ref_num; i++) this->mesh->refine_all_elements();

  // Create an H1 space with default shapeset.
  this->space = new H1Space(this->mesh, &(this->bc_types), essential_bc_values, this->init_p);
  int ndof = Space::get_num_dofs(this->space);
  info("ndof = %d", ndof);

  // Debug.
  /*
  BaseView b("", new WinGeom(0, 0, 400, 400));
  b.show(this->space);
  View::wait();
  */

  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form_vol));
  wf.add_vector_form(callback(linear_form_vol));
  for (int i=0; i < this->bdy_values_neumann.size(); i++) {
    wf.add_vector_form_surf(callback(linear_form_surf_neumann), this->bdy_markers_neumann[i]);
  }
  for (int i=0; i < this->bdy_values_newton.size(); i++) {
    wf.add_matrix_form_surf(callback(bilinear_form_surf_newton), this->bdy_markers_newton[i]);
    wf.add_vector_form_surf(callback(linear_form_surf_newton), this->bdy_markers_newton[i]);
  }

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, this->space, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling the stiffness matrix and right-hand side vector.");
  dp.assemble(matrix, rhs);

  // Solve the linear system and if successful, obtain the solution.
  info("Solving the matrix problem.");
  if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), this->space, phi);
  else error ("Matrix solver failed.\n");

  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  return true;
}

