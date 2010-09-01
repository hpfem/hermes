#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"

#include "hermes2d.h"

using namespace RefinementSelectors;

// TODO : header
// This test makes sure that the example "neutronics-4-group-adapt" works correctly.


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Number of energy discretization intervals (groups) that also defines the number of solution components, meshes, etc.,
const int N_GROUPS = 4;
#define for_each_group(g) for (int g = 0; g < N_GROUPS; g++)

const bool SOLVE_ON_COARSE_MESH = false;  // If true, coarse mesh FE problem is solved in every adaptivity step.
                                          // If false, projection of the fine mesh solution on the coarse mesh is used.
const int INIT_REF_NUM[N_GROUPS] = {      // Initial uniform mesh refinement for the individual solution components.
  1, 1, 1, 1                             
};
const int P_INIT[N_GROUPS] = {            // Initial polynomial orders for the individual solution components. 
  1, 1, 1, 1                              
};      
const double THRESHOLD = 0.3;             // This is a quantitative parameter of the adapt(...) function and
                                          // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                   // Adaptive strategy:
                                          // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                          //   error is processed. If more elements have similar errors, refine
                                          //   all to keep the mesh symmetric.
                                          // STRATEGY = 1 ... refine all elements whose error is larger
                                          //   than THRESHOLD times maximum element error.
                                          // STRATEGY = 2 ... refine all elements whose error is larger
                                          //   than THRESHOLD..
const CandList CAND_LIST = H2D_HP_ANISO;  // Predefined list of element refinement candidates. Possible values are
                                          // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                          // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                          // See User Documentation for details.
const int MESH_REGULARITY = -1;           // Maximum allowed level of hanging nodes:
                                          // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                          // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                          // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                          // Note that regular meshes are not supported, this is due to
                                          // their notoriously bad performance.
const double CONV_EXP = 1.0;              // Default value is 1.0. This parameter influences the selection of
                                          // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.1;              // Stopping criterion for adaptivity (rel. error tolerance between the
                                          // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;              // Adaptivity process stops when the number of degrees of freedom grows over
                                          // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 30;             // Adaptivity process stops when the number of adaptation steps grows over
                                          // this limit.

// Macros for simpler definition of tuples used in projections.
#define callback_pairs(a)      std::make_pair(callback(a)), std::make_pair(callback(a)), std::make_pair(callback(a)), std::make_pair(callback(a))

// Element markers.
const int marker_reflector = 1;
const int marker_core = 2;

// Boundary markers.
const int bc_vacuum = 1;
const int bc_sym = 2;

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

// Power iteration control.

double k_eff = 1.0;         // Initial eigenvalue approximation.
double TOL_PIT_CM = 5e-5;   // Tolerance for eigenvalue convergence when solving on coarse mesh.
double TOL_PIT_RM = 1e-6;   // Tolerance for eigenvalue convergence when solving on reference mesh.

// Physical data of the problem for the given number of energy groups (N_GROUPS).
#include "physical_parameters.cpp"
// Weak forms.
#include "forms.cpp"
// Norms in the axisymmetric coordinate system.
#include "norms.cpp"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// Fission source function.
inline void source_fn(int n, Tuple<scalar*> values, scalar* out)
{
  for (int i = 0; i < n; i++) {
    out[i] = 0.0;
    for_each_group(g)
      out[i] += nu[1][g] * Sf[1][g] * values.at(g)[i];
  }
}

/// Integral over the active core.
double integrate(MeshFunction* sln, int marker)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  
  double integral = 0.0;
  Element* e;
  Mesh* mesh = sln->get_mesh();
  
  for_all_active_elements(e, mesh)
  {
    if (e->marker == marker)
    {
      update_limit_table(e->get_mode());
      sln->set_active_element(e);
      RefMap* ru = sln->get_refmap();
      int o = sln->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o);
      sln->set_quad_order(o, H2D_FN_VAL);
      scalar *uval = sln->get_fn_values();
      double* x = ru->get_phys_x(o);
      double result = 0.0;
      h1_integrate_expression(x[i] * uval[i]);
      integral += result;
    }
  }
  
  return 2.0 * M_PI * integral;
}

/** \brief Power iteration. 
 *
 * Starts from an initial guess stored in the argument 'solutions' and updates it by the final result after the iteration
 * has converged, also updating the global eigenvalue 'k_eff'.
 *
 * \param[in]  spaces        Pointers to spaces on which the solutions are defined (one space for each energy group).
 * \param[in]  wf            Pointer to the weak form of the problem.
 * \param[in,out] slptr_solution   A set of Solution* pointers to solution components (neutron fluxes in each group). 
 *                                 Initial guess for the iteration on input, converged result on output.
 * \param[in,out] mfptr_solution   The same as above, only the type of the pointers is MeshFunction*.
 *                                 This is needed for the fission source filter, which accepts this type instead of Solution*.
 * \param[in]  tol           Relative difference between two successive eigenvalue approximations that stops the iteration.
 * \param[in]  matrix_solver Solver for the resulting matrix problem (one of the available types enumerated in common.h).
 * \return  number of iterations needed for convergence within the specified tolerance.
**/
int power_iteration(Tuple<Space *>& spaces, WeakForm *wf,
                    Tuple<Solution *>& slptr_solution, Tuple<MeshFunction *>& mfptr_solution,
                    double tol, MatrixSolverType matrix_solver = SOLVER_UMFPACK)
{
  // Sanity checks.
  if (slptr_solution.size() != N_GROUPS) 
    error("Wrong number of power iteration solutions for the given number of energy groups.");
  if (spaces.size() != N_GROUPS) 
    error("Spaces and solutions supplied to power_iteration do not match."); 
  if (slptr_solution.size() != mfptr_solution.size()) 
    error("Number of Solutions and corresponding MeshFunctions supplied to power_iteration does not match."); 
  
  // Initialize the linear problem.
  LinearProblem lp(wf, spaces);
  int ndof = get_num_dofs(spaces);
  
  // Select matrix solver.
  Matrix* mat; Vector* rhs; CommonSolver* solver;
  init_matrix_solver(matrix_solver, ndof, mat, rhs, solver);
  
  // The following variables will store pointers to solutions obtained at each iteration and will be needed for 
  // updating the eigenvalue. We will also need to use them in the fission source filter, so their MeshFunction* 
  // version is created as well.
  Tuple<Solution*> slptr_new_solution;
  Tuple<MeshFunction*> mfptr_new_solution;
  for_each_group(g) { 
    slptr_new_solution.push_back(new Solution);
    mfptr_new_solution.push_back(slptr_new_solution.back());
  }
  
  bool eigen_done = false; int it = 0;
  do {
    // Assemble the system matrix and rhs for the first time, then just update the rhs using the updated eigenpair.
    lp.assemble(mat, rhs, it == 0 ? false : true);
    
    // Solve the matrix problem to get a new approximation of the eigenvector.
    if (!solver->solve(mat, rhs)) error ("Matrix solver failed.\n");
    
    // Convert coefficients vector into a set of Solution pointers.
    for_each_group(g) slptr_new_solution[g]->set_fe_solution(spaces[g], rhs); 
    
    // Update fission sources.
    SimpleFilter new_source(source_fn, mfptr_new_solution);
    SimpleFilter old_source(source_fn, mfptr_solution);
    
    // Compute the eigenvalue for current iteration.
    double k_new = k_eff * (integrate(&new_source, marker_core) / integrate(&old_source, marker_core));
    
    info("      dominant eigenvalue (est): %g, rel error: %g", k_new, fabs((k_eff - k_new) / k_new));
    
    // Stopping criterion.
    if (fabs((k_eff - k_new) / k_new) < tol) eigen_done = true;
    
    // Update the final eigenvalue.
    k_eff = k_new;
    
    it++;
    
    // Store the new eigenvector approximation in the result.
    for_each_group(g) { slptr_solution[g]->copy(slptr_new_solution[g]); }
  }
  while (!eigen_done);
  
  // Free memory.
  for_each_group(g) delete slptr_new_solution[g];
  mat->free_data();
  rhs->free_data();
  //solver->free_data();  // FIXME: to be implemented. Default destructor is used now.
  delete solver;
  
  return it;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** General testing class. It can test
 *   > if the computed and expected results agree to within a specified tolerance (allowed_delta),
 *   > if the computed results do not overshoot the expected results by more than allowed_delta,
 *   > if the computed results do not undershoot the expected results by more than allowed_delta.
 * When, for the first time, the expected behavior is not achieved, 'passed' is set to false and nothing more is tested.
 * If used with a user-defined type, its class must overload all operators appearing in the intended test method.
**/
template<typename T>
struct TestSubject
{
  T allowed_delta;
  bool passed;
  
  TestSubject(T allowed_delta) : allowed_delta(allowed_delta), passed(true) {};
  
  void test_equality(const T& computed, const T& expected) { 
    if(passed) passed = (computed <= expected + allowed_delta && computed >= expected - allowed_delta);
  }
  void test_overshoot(const T& computed, const T& expected) { 
    if(passed) passed = (computed <= expected + allowed_delta);
  }
  void test_undershoot(const T& computed, const T& expected) { 
    if(passed) passed = (computed >= expected - allowed_delta);
  }
};

/// Structure that contains information about an extremum. It will be used
/// for testing the equality of computed and expected maxima, therefore it
/// overloads operators +,-,>=,<=..
struct Extremum
{
  double val, x, y;
  
  Extremum() : val(0.0), x(0.0), y(0.0) {};
  Extremum(double val, double x, double y) : val(val), x(x), y(y) {};
  
  inline Extremum operator+ (const Extremum &ex) const {
    return Extremum(val + ex.val, x + ex.x, y + ex.y);
  }
  inline Extremum operator- (const Extremum &ex) const {
    return Extremum(val - ex.val, x - ex.x, y - ex.y);
  }
  inline bool operator>= (const Extremum &ex) const {
    return (val >= ex.val && x >= ex.x && y >= ex.y);
  }
  inline bool operator<= (const Extremum &ex) const {
    return (val <= ex.val && x <= ex.x && y <= ex.y);
  }
  
  std::string str() { 
    char ret[50];
    sprintf(ret, "%lf, at (%lf,%lf)", val, x, y);
    return std::string(ret);
  }
};

/// Calculates maximum of a given function, including its coordinates.
Extremum get_peak(MeshFunction *sln)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  Element* e;
  Mesh* mesh = sln->get_mesh();
  
  scalar peak = 0.0;
  double pos_x = 0.0;
  double pos_y = 0.0;
  
  for_all_active_elements(e, mesh)
  {
    update_limit_table(e->get_mode());
    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();
    int o = sln->get_fn_order() + ru->get_inv_ref_order();
    limit_order(o);
    sln->set_quad_order(o, H2D_FN_VAL);
    scalar *uval = sln->get_fn_values();
    int np = quad->get_num_points(o);
    double* x = ru->get_phys_x(o);
    double* y = ru->get_phys_y(o);
    
    for (int i = 0; i < np; i++)
      if (uval[i] > peak) {
        peak = uval[i];
        pos_x = x[i];
        pos_y = y[i];
      }
  }
 
  return Extremum(peak, pos_x, pos_y);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();  

  // Use multimesh, i.e. create one mesh for each energy group.
  
  Tuple<Mesh *> meshes;
  for_each_group(g) meshes.push_back(new Mesh());
  
  // Load the mesh for the 1st group.
  H2DReader mloader;
  mloader.load("reactor.mesh", meshes[0]);

  for (int g = 1; g < N_GROUPS; g++) {
    // Obtain meshes for the 2nd to 4th group by cloning the mesh loaded for the 1st group.
    meshes[g]->copy(meshes[0]);
    // Initial uniform refinements.
    for (int i = 0; i < INIT_REF_NUM[g]; i++) meshes[g]->refine_all_elements();
  }
  for (int i = 0; i < INIT_REF_NUM[0]; i++) meshes[0]->refine_all_elements();
  
  // Create pointers to solutions on coarse and fine meshes and from the latest power iteration, respectively.
  Tuple<Solution*> slptr_coarse_slns, slptr_fine_slns, slptr_pow_iter_slns;
  // We will need to pass the power iteration solutions to methods like project_global,
  // which expect MeshFunction* pointers instead of just Solution*:
  Tuple<MeshFunction*> mfptr_pow_iter_slns;
  // Initialize all the new solution variables.
  for_each_group(g) 
  {
    slptr_coarse_slns.push_back(new Solution);
    slptr_fine_slns.push_back(new Solution);
    slptr_pow_iter_slns.push_back(new Solution);   
    slptr_pow_iter_slns[g]->set_const(meshes[g], 1.0);  // Starting point for the first power iteration.
    mfptr_pow_iter_slns.push_back(slptr_pow_iter_slns[g]);
  }
  // Define a macro for easier manipulation with Solution*/MeshFunction* pairs.
  #define mkptr(a) slptr_##a, mfptr_##a
  
  // Create the approximation spaces with the default shapeset.
  Tuple<Space *> spaces;
  for_each_group(g) spaces.push_back(new H1Space(meshes[g], bc_types, essential_bc_values, P_INIT[g]));

  // Initialize the weak formulation.
  WeakForm wf(N_GROUPS);
  wf.add_matrix_form(0, 0, callback(biform_0_0), H2D_SYM);
  wf.add_matrix_form(1, 1, callback(biform_1_1), H2D_SYM);
  wf.add_matrix_form(1, 0, callback(biform_1_0));
  wf.add_matrix_form(2, 2, callback(biform_2_2), H2D_SYM);
  wf.add_matrix_form(2, 1, callback(biform_2_1));
  wf.add_matrix_form(3, 3, callback(biform_3_3), H2D_SYM);
  wf.add_matrix_form(3, 2, callback(biform_3_2));
  wf.add_vector_form(0, callback(liform_0), marker_core, mfptr_pow_iter_slns);
  wf.add_vector_form(1, callback(liform_1), marker_core, mfptr_pow_iter_slns);
  wf.add_vector_form(2, callback(liform_2), marker_core, mfptr_pow_iter_slns);
  wf.add_vector_form(3, callback(liform_3), marker_core, mfptr_pow_iter_slns);
  wf.add_matrix_form_surf(0, 0, callback(biform_surf_0_0), bc_vacuum);
  wf.add_matrix_form_surf(1, 1, callback(biform_surf_1_1), bc_vacuum);
  wf.add_matrix_form_surf(2, 2, callback(biform_surf_2_2), bc_vacuum);
  wf.add_matrix_form_surf(3, 3, callback(biform_surf_3_3), bc_vacuum);
  
  // Initialize and solve coarse mesh problem.
  power_iteration(spaces, &wf, mkptr(pow_iter_slns), TOL_PIT_CM);
  // If SOLVE_ON_COARSE_MESH == true, we will store the results as the first coarse mesh solution;
  // otherwise, we will obtain this solution later by projecting the reference solution on the coarse mesh.
  if (SOLVE_ON_COARSE_MESH) 
    for_each_group(g) 
      slptr_coarse_slns[g]->copy(slptr_pow_iter_slns[g]);

  // Initialize the refinement selectors.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
  Tuple<RefinementSelectors::Selector*> selectors;
  for_each_group(g) selectors.push_back(&selector);
  
  // Adaptivity loop:
  int as = 1; bool done = false;
  double h1_err_est;
  do {
    
    info("---- Adaptivity step %d:", as);

    // Construct globally refined meshes and setup reference spaces on them.
    Tuple<Space *> ref_spaces;
    Tuple<Mesh *> ref_meshes;
    for_each_group(g) 
    { 
      ref_meshes.push_back(new Mesh());
      Mesh *ref_mesh = ref_meshes.back();      
      ref_mesh->copy(spaces[g]->get_mesh());
      ref_mesh->refine_all_elements();
      
      ref_spaces.push_back(spaces[g]->dup(ref_mesh));
      int order_increase = 1;
      ref_spaces[g]->copy_orders(spaces[g], order_increase);
    }

    // For the first time, project coarse mesh solutions on fine meshes to obtain 
    // a starting point for the fine mesh power iteration.
    if (as == 1) {
      info("Projecting initial coarse mesh solutions on fine meshes.");
      project_global(spaces, 
                      matrix_forms_tuple_t(callback_pairs(projection_biform)), 
                      vector_forms_tuple_t(callback_pairs(projection_liform)),
                      mfptr_pow_iter_slns, slptr_pow_iter_slns);
    }
    
    // Solve the fine mesh problem.
    power_iteration(ref_spaces, &wf, mkptr(pow_iter_slns), TOL_PIT_RM);
    
    // Store the results.
    for_each_group(g) slptr_fine_slns[g]->copy(slptr_pow_iter_slns[g]);
    
    // Either solve on coarse mesh or project the fine mesh solution on the coarse mesh.
    if (SOLVE_ON_COARSE_MESH) {
      if (as > 1) {
        power_iteration(spaces, &wf, mkptr(pow_iter_slns), TOL_PIT_CM);
        // Store the results.
        for_each_group(g) slptr_coarse_slns[g]->copy(slptr_pow_iter_slns[g]);
      }
    }
    else {
      info("Projecting fine mesh solutions on coarse meshes.");
      project_global(spaces, 
                      matrix_forms_tuple_t(callback_pairs(projection_biform)), 
                      vector_forms_tuple_t(callback_pairs(projection_liform)),
                      mfptr_pow_iter_slns, slptr_coarse_slns);
    }
    
    // Time measurement.
    cpu_time.tick();
    
    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    Adapt hp(spaces);
    hp.set_error_form(0, 0, callback(biform_0_0));
    hp.set_error_form(1, 1, callback(biform_1_1));
    hp.set_error_form(1, 0, callback(biform_1_0));
    hp.set_error_form(2, 2, callback(biform_2_2));
    hp.set_error_form(2, 1, callback(biform_2_1));
    hp.set_error_form(3, 3, callback(biform_3_3));
    hp.set_error_form(3, 2, callback(biform_3_2));
    
    // Calculate element errors and error estimate for adaptivity.
    info("Calculating error.");
    hp.set_solutions(slptr_coarse_slns, slptr_fine_slns);
    
    double energy_err_est = hp.calc_elem_errors(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;
    h1_err_est = error_total(error_fn_h1_axisym, norm_fn_h1_axisym, slptr_coarse_slns, slptr_fine_slns);
    double l2_err_est = error_total(error_fn_l2_axisym, norm_fn_l2_axisym, slptr_coarse_slns, slptr_fine_slns);
    
    // Time measurement.
    cpu_time.tick();
    double cta = cpu_time.accumulated();
    
    // Calculate H1 and L2 error estimates.
    std::vector<double> h1_errors, l2_errors;
    for_each_group(g) 
    {
      l2_errors.push_back( 100*l2_error_axisym(slptr_coarse_slns[g], slptr_fine_slns[g]) );
      h1_errors.push_back( 100*h1_error_axisym(slptr_coarse_slns[g], slptr_fine_slns[g]) );
    }
    
    // Millipercent eigenvalue error w.r.t. the reference value (see physical_parameters.cpp). 
    double keff_err = 1e5*fabs(k_eff - REF_K_EFF)/REF_K_EFF;
    
    cpu_time.tick(HERMES_SKIP);
    
    // If err_est too large, adapt the mesh.
    if (energy_err_est < ERR_STOP) break;
    else {
      info("Adapting the coarse mesh.");
      done = hp.adapt(selectors, THRESHOLD, STRATEGY, MESH_REGULARITY);
      if (get_num_dofs(spaces) >= NDOF_STOP) done = true;
    }
    
    // Free reference meshes and spaces.
    for_each_group(g) 
    {
      delete ref_spaces[g];
      delete ref_meshes[g];
    }
    
    as++;
    if (as >= MAX_ADAPT_NUM) done = true;
  }
  while(done == false);
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  // Display the results.
  
  info("Number of iterations: %d", as);
  info("Core eigenvalue: %lf.", k_eff);
  info("H1 error estimate: %lf.", h1_err_est);
  
  int ndofs[N_GROUPS];
  for_each_group(g) {
    ndofs[g] = spaces[g]->get_num_dofs();
    info("Number of DOFs for group %d: %d.", g+1, ndofs[g]);
  }
    
  Extremum maxima[N_GROUPS];
  for_each_group(g) { 
    maxima[g] = get_peak(slptr_coarse_slns[g]);
    info("Peak flux in group %d: %s.", g+1, maxima[g].str().c_str());
  }

  for_each_group(g) 
  {
    delete spaces[g]; delete meshes[g];
    delete slptr_coarse_slns[g], delete slptr_fine_slns[g]; delete slptr_pow_iter_slns[g];
  }
  
  // Test the results.
  
#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

  TestSubject<int> num_iter(2);
  num_iter.test_overshoot(as, 14);
  
  TestSubject<double> eigenvalue(1e-5);
  eigenvalue.test_equality(k_eff, 1.140910);
  
  TestSubject<double> error_estimate(1e-5);
  error_estimate.test_overshoot(h1_err_est, 0.002081);

  TestSubject<int> ndof(100);
  const int expected_ndofs[N_GROUPS] = {
    1204, 884, 792, 1104
  };
  for_each_group(g) ndof.test_overshoot(spaces[g]->get_num_dofs(), expected_ndofs[g]);

  TestSubject<Extremum> peak(Extremum(1e-3, 1e-3, 1e-3));
  const Extremum expected_maxima[N_GROUPS] = {
    Extremum(1.021108, 1.794896, 5.481176),
    Extremum(2.299487, 1.794896, 5.481176),
    Extremum(0.341151, 1.794896, 5.481176),
    Extremum(4.612591, 1.235793, 5.481176)
  };
  for_each_group(g) peak.test_equality(maxima[g], expected_maxima[g]);
  
  if (ndof.passed && peak.passed && eigenvalue.passed && num_iter.passed && error_estimate.passed) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
}
