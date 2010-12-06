#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "compressible-euler.h"

// Matrix solver.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_ZATECOO, SOLVER_MUMPS, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

#include "forms.cpp"
#include "filters.cpp"

// Global variables. This is necessary with the current state of callbacks.
// It should be gotten rid of as this implies that only one instance of the
// thermomechanics module can run at a time.
double dummy_kappa = 1.4;
NumericalFlux _num_flux(dummy_kappa);
double _tau, _rho_ext, _v1_ext, _v2_ext, _p_ext;
int _solid_wall_marker, _inlet_outlet_marker;

// Constructor.
CompressibleEuler::CompressibleEuler(bool use_vector_valued_forms, bool calculate_time_derivative_approximation) :
  use_vector_valued_forms(use_vector_valued_forms), calculate_time_derivative_approximation(calculate_time_derivative_approximation),
    init_p(-1,-1)
{
  init_ref_num = -1;
  mesh = new Mesh;
  spaces = NULL;
  bc_types = new BCTypes();
  CFL = tau =  T = p_ext = rho_ext = v1_ext = v2_ext = -1;
  num_flux = NULL;
  filters_for_visual = false;
  saving_frequency = 0;
}

// Destructor.
CompressibleEuler::~CompressibleEuler()
{
  if(mesh != NULL)
    delete this->mesh;
  if(spaces != NULL) {
    for(unsigned int i = 0; i < spaces->size(); i++)
      delete (*this->spaces)[i];
    delete spaces;
  }
}

// Set mesh via a string. See CompressibleEuler.h for an example of such a string.
void CompressibleEuler::set_mesh_str(const std::string &mesh)
{
  this->mesh_str = mesh;
}

void CompressibleEuler::set_mesh_filename(const std::string &mesh)
{
  this->mesh_filename = mesh;
}

// Set the number of initial uniform mesh refinements.
void CompressibleEuler::set_initial_mesh_refinement(int init_ref_num) 
{
  this->init_ref_num = init_ref_num;
}

// Set initial poly degree in elements.
void CompressibleEuler::set_initial_poly_degree(Ord2 p) 
{
  this->init_p = p;
}
void CompressibleEuler::set_initial_poly_degree(int p) 
{
  this->init_p = Ord2(p, p);
}

void CompressibleEuler::set_CFL_value(double CFL)
{
  this->CFL = CFL;
}

void CompressibleEuler::set_time_step(double tau_)
{
  this->tau = tau_;
  _tau = tau_;
}

void CompressibleEuler::set_time_interval_length(double T)
{
  this->T = T;
}

void CompressibleEuler::set_state(double p_ext_, double rho_ext_, double v1_ext_, double v2_ext_)
{
  this->p_ext = p_ext_;
  this->rho_ext = rho_ext_;
  this->v1_ext = v1_ext_;
  this->v2_ext = v2_ext_;

  _p_ext = p_ext_;
  _rho_ext = rho_ext_;
  _v1_ext = v1_ext_;
  _v2_ext = v2_ext_;
}

void CompressibleEuler::set_kappa(double kappa)
{
  this->num_flux = new NumericalFlux(kappa); 
  _num_flux = *(this->num_flux);
}

void CompressibleEuler::set_inlet_outlet_marker(int marker)
{
  this->inlet_outlet_marker = marker;
  _inlet_outlet_marker = marker;
}

void CompressibleEuler::set_solid_wall_marker(int marker)
{
  this->solid_wall_marker = marker;
  _solid_wall_marker = marker;
}

void CompressibleEuler::save_every_th_solution_visualization(unsigned int frequency) 
{ this->saving_frequency = frequency; }

void CompressibleEuler::set_use_vector_valued_forms(bool set) 
{ this->use_vector_valued_forms = set; }

void CompressibleEuler::use_filters_for_visual(bool set) 
{ this->filters_for_visual = set; }

void CompressibleEuler::set_calculate_time_derivative_approximation(bool set) 
{ this->calculate_time_derivative_approximation = set; }

bool CompressibleEuler::all_set()
{
  if(init_ref_num == -1)
    return false;
  if(init_p.order_h == -1 || init_p.order_v == -1)
    return false;
  if(CFL == -1 || tau == -1 || T == -1 || p_ext == -1 || 
     rho_ext == -1 || v1_ext == -1 || v2_ext == -1)
    return false;
  if(num_flux == NULL)
    return false;
  return true;
}

// Solve the problem.
bool CompressibleEuler::calculate(Hermes::Tuple<Solution*> sln) 
{
  // Initial check.
  if(!all_set())
    return false;

  // Load the mesh.
  H2DReader mloader;
  mloader.load(this->mesh_filename.c_str(), this->mesh);

  // Perform initial uniform mesh refinements.
  for (int i = 0; i < this->init_ref_num; i++) this->mesh->refine_all_elements();

  // Create L2 spaces with default shapesets.
  this->spaces = new Hermes::Tuple<Space*>();
  this->spaces->push_back(new L2Space(this->mesh, this->bc_types, NULL, this->init_p));
  this->spaces->push_back(new L2Space(this->mesh, this->bc_types, NULL, this->init_p));
  this->spaces->push_back(new L2Space(this->mesh, this->bc_types, NULL, this->init_p));
  this->spaces->push_back(new L2Space(this->mesh, this->bc_types, NULL, this->init_p));
  
  int ndof = Space::get_num_dofs(*(this->spaces));
  info("ndof = %d", ndof);

  // Initialize solutions, set initial conditions.
  Solution sln_rho, sln_rho_v_x, sln_rho_v_y, sln_e, prev_rho, prev_rho_v_x, prev_rho_v_y, prev_e;
  sln_rho.set_exact(mesh, ic_density);
  sln_rho_v_x.set_exact(mesh, ic_density_vel_x);
  sln_rho_v_y.set_exact(mesh, ic_density_vel_y);
  sln_e.set_exact(mesh, ic_energy);
  prev_rho.set_exact(mesh, ic_density);
  prev_rho_v_x.set_exact(mesh, ic_density_vel_x);
  prev_rho_v_y.set_exact(mesh, ic_density_vel_y);
  prev_e.set_exact(mesh, ic_energy);

  // Initialize weak formulation.
  WeakForm wf(4);

  // Bilinear forms coming from time discretization by explicit Euler's method.
  wf.add_matrix_form(0, 0, callback(bilinear_form_0_0_time));
  wf.add_matrix_form(1, 1, callback(bilinear_form_1_1_time));
  wf.add_matrix_form(2, 2, callback(bilinear_form_2_2_time));
  wf.add_matrix_form(3, 3, callback(bilinear_form_3_3_time));

  // Volumetric linear forms.
  // Linear forms coming from the linearization by taking the Eulerian fluxes' Jacobian matrices 
  // from the previous time step.
  // First flux.
  // Unnecessary for FVM.
  if(init_p.order_h > 0 || init_p.order_v > 0) {
    wf.add_vector_form(0, callback(linear_form_0_1), HERMES_ANY, Hermes::Tuple<MeshFunction*>(&prev_rho_v_x));
    wf.add_vector_form(1, callback(linear_form_1_0_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_1_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_2_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_3_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(2, callback(linear_form_2_0_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_1_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_2_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_3_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_0_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_1_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_2_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_3_first_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    // Second flux.
    
    wf.add_vector_form(0, callback(linear_form_0_2), HERMES_ANY, Hermes::Tuple<MeshFunction*>(&prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_0_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_1_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_2_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(1, callback(linear_form_1_3_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(2, callback(linear_form_2_0_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_1_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_2_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y));
    wf.add_vector_form(2, callback(linear_form_2_3_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_0_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_1_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_2_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, callback(linear_form_3_3_second_flux), HERMES_ANY, 
                       Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  }

  // Volumetric linear forms coming from the time discretization.
  if(use_vector_valued_forms) {
    wf.add_vector_form(0, linear_form_vector, linear_form_order, HERMES_ANY, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(1, linear_form_vector, linear_form_order, HERMES_ANY, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(2, linear_form_vector, linear_form_order, HERMES_ANY, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form(3, linear_form_vector, linear_form_order, HERMES_ANY, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  }
  else {
    wf.add_vector_form(0, linear_form, linear_form_order, HERMES_ANY, &prev_rho);
    wf.add_vector_form(1, linear_form, linear_form_order, HERMES_ANY, &prev_rho_v_x);
    wf.add_vector_form(2, linear_form, linear_form_order, HERMES_ANY, &prev_rho_v_y);
    wf.add_vector_form(3, linear_form, linear_form_order, HERMES_ANY, &prev_e);
  }

  // Surface linear forms - inner edges coming from the DG formulation.
  if(use_vector_valued_forms) {
    wf.add_vector_form_surf(0, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(1, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(2, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(3, linear_form_interface_vector, linear_form_order, H2D_DG_INNER_EDGE, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  }
  else {
    wf.add_vector_form_surf(0, linear_form_interface_0, linear_form_order, H2D_DG_INNER_EDGE, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(1, linear_form_interface_1, linear_form_order, H2D_DG_INNER_EDGE, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(2, linear_form_interface_2, linear_form_order, H2D_DG_INNER_EDGE, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(3, linear_form_interface_3, linear_form_order, H2D_DG_INNER_EDGE, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  }

  // Surface linear forms - inlet / outlet edges.
  if(use_vector_valued_forms) {
    wf.add_vector_form_surf(0, bdy_flux_inlet_outlet_comp_vector, linear_form_order, inlet_outlet_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(1, bdy_flux_inlet_outlet_comp_vector, linear_form_order, inlet_outlet_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(2, bdy_flux_inlet_outlet_comp_vector, linear_form_order, inlet_outlet_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(3, bdy_flux_inlet_outlet_comp_vector, linear_form_order, inlet_outlet_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  }
  else {
    wf.add_vector_form_surf(0, bdy_flux_inlet_outlet_comp_0, linear_form_order, inlet_outlet_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(1, bdy_flux_inlet_outlet_comp_1, linear_form_order, inlet_outlet_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(2, bdy_flux_inlet_outlet_comp_2, linear_form_order, inlet_outlet_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(3, bdy_flux_inlet_outlet_comp_3, linear_form_order, inlet_outlet_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  }
  
  // Surface linear forms - Solid wall edges.
  if(use_vector_valued_forms) {
    wf.add_vector_form_surf(0, bdy_flux_solid_wall_comp_vector, linear_form_order, solid_wall_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(1, bdy_flux_solid_wall_comp_vector, linear_form_order, solid_wall_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(2, bdy_flux_solid_wall_comp_vector, linear_form_order, solid_wall_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(3, bdy_flux_solid_wall_comp_vector, linear_form_order, solid_wall_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  }
  else {
    wf.add_vector_form_surf(0, bdy_flux_solid_wall_comp_0, linear_form_order, solid_wall_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(1, bdy_flux_solid_wall_comp_1, linear_form_order, solid_wall_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(2, bdy_flux_solid_wall_comp_2, linear_form_order, solid_wall_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
    wf.add_vector_form_surf(3, bdy_flux_solid_wall_comp_3, linear_form_order, solid_wall_marker, 
                            Hermes::Tuple<MeshFunction*>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e));
  }

  // Initialize the FE problem.
  bool is_linear = true;
  DiscreteProblem dp(&wf, *spaces, is_linear);

  // If the FE problem is in fact a FV problem.
  if(init_p.order_h == 0 && init_p.order_v == 0)
    dp.set_fvm();

  // If the vector valued forms should be used.
  if(use_vector_valued_forms)
    dp.use_vector_valued_forms();

  // Filters for visualization of pressure and the two components of velocity.
  SimpleFilter pressure(calc_pressure_func, Hermes::Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter u(calc_u_func, Hermes::Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  SimpleFilter w(calc_w_func, Hermes::Tuple<MeshFunction*>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
  VectorView vview("Velocity", new WinGeom(0, 0, 600, 300));
  ScalarView sview("Pressure", new WinGeom(700, 0, 600, 300));

  // Visualization of the solution components.
  ScalarView s1("w1", new WinGeom(0, 0, 620, 300));
  s1.fix_scale_width(80);
  ScalarView s2("w2", new WinGeom(625, 0, 600, 300));
  s2.fix_scale_width(50);
  ScalarView s3("w3", new WinGeom(0, 350, 620, 300));
  s3.fix_scale_width(80);
  ScalarView s4("w4", new WinGeom(625, 350, 600, 300));
  s4.fix_scale_width(50);

  // Iteration number.
  int iteration = 0;
  
  // Set up the solver, matrix, and rhs according to the solver selection.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

  // Output of the approximate time derivative.
  std::ofstream time_der_out("time_der");

  for(double t = 0.0; t < T; t += tau)
  {
    info("---- Time step %d, time %3.5f.", iteration, t);

    iteration++;

    bool rhs_only = (iteration == 1 ? false : true);
    // Assemble stiffness matrix and rhs or just rhs.
    if (rhs_only == false) info("Assembling the stiffness matrix and right-hand side vector.");
    else info("Assembling the right-hand side vector (only).");
    dp.assemble(matrix, rhs, rhs_only);

    // Solve the matrix problem.
    info("Solving the matrix problem.");
    if(solver->solve())
      Solution::vector_to_solutions(solver->get_solution(), *spaces, Hermes::Tuple<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e));
    else
    error ("Matrix solver failed.\n");

    // Approximate the time derivative of the solution.
    if(calculate_time_derivative_approximation) {
      Adapt *adapt_for_time_der_calc = new Adapt(*spaces, Hermes::Tuple<ProjNormType>(HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM, HERMES_L2_NORM));
      bool solutions_for_adapt = false;
      double difference = adapt_for_time_der_calc->calc_err_est(Hermes::Tuple<Solution *>(&prev_rho, &prev_rho_v_x, &prev_rho_v_y, &prev_e), 
        Hermes::Tuple<Solution *>(&sln_rho, &sln_rho_v_x, &sln_rho_v_y, &sln_e), solutions_for_adapt, HERMES_TOTAL_ERROR_ABS | HERMES_ELEMENT_ERROR_ABS) / tau;
      delete adapt_for_time_der_calc;

      // Info about the approximate time derivative.
      if(iteration > 1)
      {
        info("Approximate the norm time derivative : %g.", difference);
        time_der_out << iteration << '\t' << difference << std::endl;
      }
    }

    // Determine the time step according to the CFL condition.
    // Only mean values on an element of each solution component are taken into account.
    double *solution_vector = solver->get_solution();
    double min_condition = 0;
    Element *e;
    for (int _id = 0, _max = mesh->get_max_element_id(); _id < _max; _id++) \
          if (((e) = mesh->get_element_fast(_id))->used) \
            if ((e)->active)
    {
      AsmList al;
      (*spaces)[0]->get_element_assembly_list(e, &al);
      double rho = solution_vector[al.dof[0]];
      (*spaces)[1]->get_element_assembly_list(e, &al);
      double v1 = solution_vector[al.dof[0]] / rho;
      (*spaces)[2]->get_element_assembly_list(e, &al);
      double v2 = solution_vector[al.dof[0]] / rho;
      (*spaces)[3]->get_element_assembly_list(e, &al);
      double energy = solution_vector[al.dof[0]];
      
      double condition = e->get_area() / (std::sqrt(v1*v1 + v2*v2) + calc_sound_speed(rho, rho*v1, rho*v2, energy));
      
      if(condition < min_condition || min_condition == 0.)
        min_condition = condition;
    }
    if(tau > min_condition)
      tau = min_condition;
    if(tau < min_condition * 0.9)
      tau = min_condition;

    // Copy the solutions into the previous time level ones.
    prev_rho.copy(&sln_rho);
    prev_rho_v_x.copy(&sln_rho_v_x);
    prev_rho_v_y.copy(&sln_rho_v_y);
    prev_e.copy(&sln_e);

    // Visualization.
    if(filters_for_visual) {
      pressure.reinit();
      u.reinit();
      w.reinit();
      sview.show(&pressure);
      vview.show(&u, &w);
      if(saving_frequency > 0 && iteration % saving_frequency == 0) {
        sview.save_numbered_screenshot("pressure_%d", iteration);
        vview.save_numbered_screenshot("velocity_%d", iteration);
      }
    }
    else {
      s1.show(&sln_rho);
      s2.show(&sln_rho_v_x);
      s3.show(&sln_rho_v_y);
      s4.show(&sln_e);
      if(saving_frequency > 0 && iteration % saving_frequency == 0) {
        s1.save_numbered_screenshot("rho_%d.bmp", iteration);
        s2.save_numbered_screenshot("rho_v_x%d.bmp", iteration);
        s3.save_numbered_screenshot("rho_v_y%d.bmp", iteration);
        s4.save_numbered_screenshot("e_%d.bmp", iteration);
      }
    }
    
    // If used, we need to clean the vector valued form caches.
    if(use_vector_valued_forms) DiscreteProblem::empty_form_caches();
  }
  time_der_out.close();
  return 0;
}
