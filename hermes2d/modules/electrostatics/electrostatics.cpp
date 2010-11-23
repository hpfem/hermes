// Constructor.
Electrostatics::Electrostatics()
{
  mesh_filename = NULL;
  init_ref_num = -1;
  init_p = -1;
  n_mat_markers = -1;
  mat_markers = NULL;
  permittivity_array = NULL;
  charge_density_array = NULL;
  n_bc_value = -1;
  bc_markers_value = NULL;
  bc_values = NULL;
  n_bc_derivative = -1;
  bc_markers_derivative = NULL;
  bc_derivative = NULL;
  mesh = NULL;
  space = NULL;
}

// Destructor.
Electrostatics::~Electrostatics()
{
  delete this->mesh;
  delete this->space;
  delete [] this->mesh_filename;
}

// Set mesh file name.
bool Electrostatics::set_mesh_filename(char* filename)
{
  strcpy(filename, this->mesh_filename);
}

// Set the number of initial uniform mesh refinements.
void set_initial_mesh_refinement(int init_ref_num) 
{
  this->init_ref_num = init_ref_num;
}

// Set initial poly degree in elements.
void Electrostatics::set_initial_poly_degree(int p) 
{
  this->init_p = p;
}

// Set material markers, and check compatibility with mesh file.
void Electrostatics::set_material_markers(int* mat_markers)
{

}

// Set permittivity array.
void Electrostatics set_permittivity_array(double* p_array)
{

}

// Set charge density array.
void Electrostatics set_charge_density_array(double* cd_array)
{

}

// Set VALUE boundary markers (also check with the mesh file).
void Electrostatics set_boundary_markers_value(int* bdy_markers_val)
{

}

// Set boundary values.
void Electrostatics set_boundary_values(double* bc_val)
{

}

// Set DERIVATIVE boundary markers (also check with the mesh file).
void Electrostatics::set_boundary_markers_derivative(int* bc_markers_der)
{

}

// Set boundary derivatives.
void Electrostatics set_boundary_derivatives(double* bc_der)
{

}

// Solve the problem.
void Electrostatics::calculate(Solution* phi) 
{
  // Load the mesh.
  H2DReader mloader;
  mloader.load(this->mesh_filename, this->mesh);

  // Perform initial uniform mesh refinements.
  for (int i = 0; i < this->init_ref_num; i++) this->mesh.refine_all_elements();

  // Create an H1 space with default shapeset.
  H1Space space(this->mesh, this->bc_types, this->essential_bc_values, this->init_p);
  int ndof = Space::get_num_dofs(this->space);
  info("ndof = %d", ndof);

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







// Weak forms.
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, 
                     Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                   Geom<Real> *e, ExtData<Scalar> *ext)
{
  return CONST_F*int_v<Real, Scalar>(n, wt, v);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, 
                        Geom<Real> *e, ExtData<Scalar> *ext)
{
  return CONST_GAMMA[e->marker - 1] * int_v<Real, Scalar>(n, wt, v);
}
