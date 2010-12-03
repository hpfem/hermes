#include "electrostatics.h"

// FIXME: Because of callbacks that depend on global variables
// we need to define this here, but it is a temporary solution.
// Like this, one cannot have two instances of the Electrostatics
// module -- their global variables would interfere with each other.
std::vector<int> _global_mat_markers;
std::vector<double> _global_permittivity_array;
std::vector<double> _global_charge_density_array;
std::vector<double> _global_bc_val;
std::vector<double> _global_bc_der;
std::vector<int> _global_mat_permut;
std::vector<int> _global_bc_permut;
BCTypes* _global_bc_types=NULL;

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.

// Weak forms.
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                     Geom<Real> *e, ExtData<Scalar> *ext)
{
  int mat_marker = e->marker;
  double permittivity = _global_permittivity_array[_global_mat_permut[mat_marker]];
  return permittivity * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                   Geom<Real> *e, ExtData<Scalar> *ext)
{
  int mat_marker = e->marker;
  double charge_density = _global_charge_density_array[_global_mat_permut[mat_marker]];
  return charge_density * int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                        Geom<Real> *e, ExtData<Scalar> *ext)
{
  int edge_marker = e->marker;
  double surf_charge_density = _global_bc_der[_global_bc_types->find_index_natural(edge_marker)];
  return surf_charge_density * int_v<Real, Scalar>(n, wt, v);
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

/*
   // DEPRECATED
// Boundary condition types.
// Note: "essential" means that solution value is prescribed.
BCType bc_types(int marker)
{
    if (marker == 0)
        return BC_NONE;
    int idx;
    if (find_index(_global_bdy_markers_val, marker, idx))
        return BC_ESSENTIAL;
    if (find_index(_global_bdy_markers_der, marker, idx))
        return BC_NATURAL;
    error("Wrong boundary marker");
}
*/

void *_global_data;

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
    Electrostatics *self = (Electrostatics *) _global_data;
    int idx = self->bc_types.find_index_essential(ess_bdy_marker);
    return _global_bc_val[idx];
}

// Constructor.
Electrostatics::Electrostatics()
{
  init_ref_num = -1;
  init_p = -1;
  mesh = new Mesh();
  space = NULL;

  // FIXME: fix these:
  _global_data = this;
  _global_bc_types = &(this->bc_types);
}

// Destructor.
Electrostatics::~Electrostatics()
{
  delete this->mesh;
  delete this->space;
}

// Set mesh via a string. See electrostatics.h for an example of such string.
void Electrostatics::set_mesh_str(const std::string &mesh)
{
    this->mesh_str = mesh;
}

// Set the number of initial uniform mesh refinements.
void Electrostatics::set_initial_mesh_refinement(int init_ref_num) 
{
  this->init_ref_num = init_ref_num;
}

// Set initial poly degree in elements.
void Electrostatics::set_initial_poly_degree(int p) 
{
  this->init_p = p;
}

// Set material markers, and check compatibility with mesh file.
void Electrostatics::set_material_markers(const std::vector<int> &m_markers)
{
    this->mat_markers = m_markers;
    _global_mat_markers = m_markers;
}

// Set permittivity array.
void Electrostatics::set_permittivity_array(const std::vector<double> &p_array)
{
    this->permittivity_array = p_array;
    _global_permittivity_array = p_array;
}

// Set charge density array.
void Electrostatics::set_charge_density_array(const std::vector<double> &cd_array)
{
    this->charge_density_array = cd_array;
    _global_charge_density_array = cd_array;
}

// Set VALUE boundary markers (also check with the mesh file).
void Electrostatics::set_boundary_markers_value(const std::vector<int>
            &bdy_markers_val)
{
    Tuple<int>  t;
    (std::vector<int>)t = bdy_markers_val;
    printf("assigning....");
    t.print();
    printf("but:");
    printf("%d ", bdy_markers_val.size());
    printf("%d ", bdy_markers_val[0]);
    this->bc_types.add_bc_essential(bdy_markers_val);
}

// Set boundary values.
void Electrostatics::set_boundary_values(const std::vector<double> &bc_val)
{
    this->bc_val = bc_val;
    _global_bc_val = bc_val;
}

// Set DERIVATIVE boundary markers (also check with the mesh file).
void Electrostatics::set_boundary_markers_derivative(const std::vector<int> &bdy_markers_der)
{
    this->bc_types.add_bc_natural(bdy_markers_der);
}

// Set boundary derivatives.
void Electrostatics::set_boundary_derivatives(const std::vector<double> &bc_der)
{
    this->bc_der = bc_der;
    _global_bc_der = bc_der;
}


// Solve the problem.
bool Electrostatics::calculate(Solution* phi) 
{
  /* SANITY CHECKS */

  this->bc_types.check_consistency();

  // Sanity check of material markers and material constants.
  int n_mat_markers = this->mat_markers.size();
  if (n_mat_markers != this->permittivity_array.size()) error("Wrong number of permittivities.");
  if (n_mat_markers != this->charge_density_array.size()) error("Wrong number of charge densities.");
  // Making sure that they are nonnegative (>= 0).
  for (int i=0; i < n_mat_markers; i++) {
    if(this->mat_markers[i] < 0) error("Material markers must be nonnegative.");
  }

  /* CREATE THE PERMUTATION ARRAY mat_permut[] */

  // Get maximum material marker.
  int max_mat_marker = -1;
  for (int i=0; i < n_mat_markers; i++) {
    if (this->mat_markers[i] > max_mat_marker) max_mat_marker = this->mat_markers[i];
  }
  // Create the permutation array and initiate it with minus ones.
  for (int i=0; i < max_mat_marker+1; i++) this->mat_permut.push_back(-1);
  // Fill it.
  for (int i=0; i < n_mat_markers; i++) this->mat_permut[this->mat_markers[i]] = i;

  /* CREATE THE PERMUTATION ARRAY bc_permut[] */

  // Define global permutation arrays.
  _global_mat_permut = this->mat_permut;

  /*
  // debug: print permutation arrays
  printf("mat_permut = ");
  for (int i=0; i < max_mat_marker+1; i++) printf("%d ", this->mat_permut[i]);
  printf("\nbc_permut = ");
  for (int i=0; i < max_bdy_marker+1; i++) printf("%d ", this->bc_permut[i]);
  */

  /* BEGIN THE COMPUTATION */

  // Load the mesh.
  H2DReader mloader;
  mloader.load_str(this->mesh_str.c_str(), this->mesh);
  /*
  MeshView m("", 0, 0, 400, 400);
  m.show(this->mesh);
  View::wait();
  */

  // Perform initial uniform mesh refinements.
  for (int i = 0; i < this->init_ref_num; i++) this->mesh->refine_all_elements();

  // Create an H1 space with default shapeset.
  this->space = new H1Space(this->mesh, &(this->bc_types), essential_bc_values,
          this->init_p);
  int ndof = Space::get_num_dofs(this->space);
  info("ndof = %d", ndof);
  /*
  BaseView b("", new WinGeom(0, 0, 400, 400));
  b.show(this->space);
  View::wait();
  */


  // Initialize the weak formulation.
  WeakForm wf;
  wf.add_matrix_form(callback(bilinear_form));
  wf.add_vector_form(callback(linear_form));
  wf.add_vector_form_surf(callback(linear_form_surf));

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

