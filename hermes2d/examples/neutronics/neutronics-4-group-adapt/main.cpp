#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"

#include "hermes2d.h"

using namespace RefinementSelectors;

// This example uses automatic adaptivity to solve a 4-group neutron diffusion equation in the reactor core.
// The eigenproblem is solved using power interations.
//
// The reactor neutronics in a general coordinate system is given by the following eigenproblem:
//
//  - \nabla \cdot D_g \nabla \phi_g + \Sigma_{Rg}\phi_g - \sum_{g' \neq g} \Sigma_s^{g'\to g} \phi_{g'} =
//  = \frac{\chi_g}{k_{eff}} \sum_{g'} \nu_{g'} \Sigma_{fg'}\phi_{g'}
//
// where 1/k_{eff} is eigenvalue and \phi_g, g = 1,...,4 are eigenvectors (neutron fluxes). The current problem
// is posed in a 3D cylindrical axisymmetric geometry, leading to a 2D problem with r-z as the independent spatial
// coordinates. Identifying r = x, z = y, the gradient in the weak form has the same components as in the
// x-y system, while all integrands are multiplied by 2\pi x (determinant of the transformation matrix).
//
// BC:
//
// homogeneous neumann on symmetry axis
// d D_g\phi_g / d n = - 0.5 \phi_g   elsewhere
//
// The eigenproblem is numerically solved using common technique known as the power method (power iterations):
//
//  1) Make an initial estimate of \phi_g and k_{eff}
//  2) For n = 1, 2,...
//         solve for \phi_g using previous k_prev
//         solve for new k_{eff}
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{new}
//               k_new =  k_prev -------------------------------------------------------------------------
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{prev}
//  3) Stop iterations when
//
//     |   k_new - k_prev  |
//     | ----------------- |  < epsilon
//     |       k_new       |
//
//  Author: Milan Hanus (University of West Bohemia, Pilsen, Czech Republic).



// Number of energy discretization intervals (groups) that also defines the number of solution components, meshes, etc.,
const int N_GROUPS = 4;
#define for_each_group(g) for (int g = 0; g < N_GROUPS; g++)

const bool SOLVE_ON_COARSE_MESH = false; // If true, coarse mesh FE problem is solved in every adaptivity step.
                                         // If false, projection of the fine mesh solution on the coarse mesh is used.
const int INIT_REF_NUM[N_GROUPS] = {     // Initial uniform mesh refinement for the individual solution components.
  1, 1, 1, 1                             
};
const int P_INIT[N_GROUPS] = {           // Initial polynomial orders for the individual solution components. 
  1, 1, 1, 1                              
};      
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1e-1;             // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 30;            // Adaptivity process stops when the number of adaptation steps grows over
                                         // this limit.

// Macros for simpler definition of tuples used in projections.
#define callback_pairs(a)      std::make_pair(callback(a)), std::make_pair(callback(a)), std::make_pair(callback(a)), std::make_pair(callback(a))

// Element markers.
const int marker_reflector = 1;
const int marker_core = 2;

// Boundary markers.
const int BDY_VACUUM = 1;
const int BDY_SYM = 2;

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

/// Fission source function.
inline void source_fn(int n, Hermes::vector<scalar*> values, scalar* out)
{
  for (int i = 0; i < n; i++) {
    out[i] = 0.0;
    for_each_group(g)
      out[i] += nu[1][g] * Sf[1][g] * values.at(g)[i];
  }
}

// Integral over the active core.
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

/// Calculate number of negative solution values.
int get_num_of_neg(MeshFunction *sln)
{
	Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  Element* e;
  Mesh* mesh = sln->get_mesh();

  int n = 0;

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

		for (int i = 0; i < np; i++)
			if (uval[i] < -1e-12)
				n++;
  }

  return n;
} 

/// \brief Power iteration. 
///
/// Starts from an initial guess stored in the argument 'solutions' and updates it by the final result after the iteration
/// has converged, also updating the global eigenvalue 'k_eff'.
///
/// \param[in]  spaces        Pointers to spaces on which the solutions are defined (one space for each energy group).
/// \param[in]  wf            Pointer to the weak form of the problem.
/// \param[in,out] slptr_solution   A set of Solution* pointers to solution components (neutron fluxes in each group). 
///                                 Initial guess for the iteration on input, converged result on output.
/// \param[in,out] mfptr_solution   The same as above, only the type of the pointers is MeshFunction*.
///                                 This is needed for the fission source filter, which accepts this type instead of Solution*.
/// \param[in]  tol           Relative difference between two successive eigenvalue approximations that stops the iteration.
/// \param[in,out] mat        Pointer to a matrix to which the system associated with the power iteration will be assembled.
/// \param[in,out] rhs        Pointer to a vector to which the right hand sides of the power iteration will be successively assembled.
/// \param[in]  solver        Solver for the resulting matrix problem (specified by \c mat and \c rhs).
///
/// \return  number of iterations needed for convergence within the specified tolerance.
///
int power_iteration(Hermes::vector<Space *>& spaces, WeakForm *wf,
                    Hermes::vector<Solution *>& slptr_solution, Hermes::vector<MeshFunction *>& mfptr_solution,
                    double tol, SparseMatrix *mat, Vector* rhs, Solver *solver)
{
  // Sanity checks.
  if (slptr_solution.size() != (unsigned) N_GROUPS) 
    error("Wrong number of power iteration solutions for the given number of energy groups.");
  if (spaces.size() != (unsigned) N_GROUPS) 
    error("Spaces and solutions supplied to power_iteration do not match."); 
  if (slptr_solution.size() != mfptr_solution.size()) 
    error("Number of Solutions and corresponding MeshFunctions supplied to power_iteration does not match."); 
  
  // Initialize the linear problem.
  bool is_linear = true;
  DiscreteProblem dp(wf, spaces, is_linear);
  int ndof = Space::get_num_dofs(spaces);
    
  // The following variables will store pointers to solutions obtained at each iteration and will be needed for 
  // updating the eigenvalue. We will also need to use them in the fission source filter, so their MeshFunction* 
  // version is created as well.
  Hermes::vector<Solution*> slptr_new_solution;
  Hermes::vector<MeshFunction*> mfptr_new_solution;
  for_each_group(g) { 
    slptr_new_solution.push_back(new Solution);
    mfptr_new_solution.push_back(slptr_new_solution.back());
  }
  
  // This power iteration will most probably run on a different mesh than the previous one and so will be different
  // the corresponding algebraic system. We will need to factorize it anew (but then, the L and U factors may be 
  // reused until the next adaptation changes the mesh again).
  // TODO: This could be solved more elegantly by defining a function Solver::reinit().
  solver->set_factorization_scheme(HERMES_FACTORIZE_FROM_SCRATCH);
  
  bool eigen_done = false; int it = 0;
  do {
    // Assemble the system matrix and rhs for the first time, then just update the rhs using the updated eigenpair.
    dp.assemble(mat, rhs, it == 0 ? false : true);
        
    // Solve the matrix problem to get a new approximation of the eigenvector.
    if (!solver->solve()) error ("Matrix solver failed.\n");
    
    // The matrix doesn't change within the power iteration loop, so the first computed LU factorization may be
    // completely reused in following iterations.
    solver->set_factorization_scheme(HERMES_REUSE_FACTORIZATION_COMPLETELY);
    
    // Convert coefficients vector into a set of Solution pointers.
    Solution::vector_to_solutions(solver->get_solution(), spaces, slptr_new_solution);
 
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
  
  return it;
}


// Macros for simpler reporting (four group case).
#define report_num_dofs(spaces) spaces[0]->Space::get_num_dofs(), spaces[1]->Space::get_num_dofs(),\
                                spaces[2]->Space::get_num_dofs(), spaces[3]->Space::get_num_dofs(), Space::get_num_dofs(spaces)
#define report_errors(errors) errors[0],errors[1],errors[2],errors[3]

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Use multimesh, i.e. create one mesh for each energy group.
  Hermes::vector<Mesh *> meshes;
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
  Hermes::vector<Solution*> slptr_coarse_slns, slptr_fine_slns, slptr_pow_iter_slns;
  // We will need to pass the power iteration solutions to methods like OGProjection::project_global,
  // which expect MeshFunction* pointers instead of just Solution*:
  Hermes::vector<MeshFunction*> mfptr_pow_iter_slns;
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

  // Enter boundary markers.
  BCTypes bc_types;
  bc_types.add_bc_neumann(BDY_SYM);
  bc_types.add_bc_newton(BDY_VACUUM);
  
  // Create the approximation spaces with the default shapeset.
  Hermes::vector<Space *> spaces;
  for_each_group(g) spaces.push_back(new H1Space(meshes[g], &bc_types, P_INIT[g]));

  // Initialize the weak formulation.
  WeakForm wf(N_GROUPS);
  wf.add_matrix_form(0, 0, callback(biform_0_0), HERMES_SYM);
  wf.add_matrix_form(1, 1, callback(biform_1_1), HERMES_SYM);
  wf.add_matrix_form(1, 0, callback(biform_1_0));
  wf.add_matrix_form(2, 2, callback(biform_2_2), HERMES_SYM);
  wf.add_matrix_form(2, 1, callback(biform_2_1));
  wf.add_matrix_form(3, 3, callback(biform_3_3), HERMES_SYM);
  wf.add_matrix_form(3, 2, callback(biform_3_2));
  wf.add_vector_form(0, callback(liform_0), marker_core, mfptr_pow_iter_slns);
  wf.add_vector_form(1, callback(liform_1), marker_core, mfptr_pow_iter_slns);
  wf.add_vector_form(2, callback(liform_2), marker_core, mfptr_pow_iter_slns);
  wf.add_vector_form(3, callback(liform_3), marker_core, mfptr_pow_iter_slns);
  wf.add_matrix_form_surf(0, 0, callback(biform_surf_0_0), BDY_VACUUM);
  wf.add_matrix_form_surf(1, 1, callback(biform_surf_1_1), BDY_VACUUM);
  wf.add_matrix_form_surf(2, 2, callback(biform_surf_2_2), BDY_VACUUM);
  wf.add_matrix_form_surf(3, 3, callback(biform_surf_3_3), BDY_VACUUM);
    
  // Initialize the discrete algebraic representation of the problem and its solver.
  //
  // Choose one of the available linear algebraic system solvers (possibilities:
  // SOLVER_UMFPACK, SOLVER_PETSC, SOLVER_MUMPS, SOLVER_SUPERLU, SOLVER_AMESOS, 
  // SOLVER_AZTECOO, depending on which solver libraries you have installed and enabled in Hermes).
#ifdef WITH_PETSC  
  MatrixSolverType matrix_solver = SOLVER_PETSC;
#else
  MatrixSolverType matrix_solver = SOLVER_UMFPACK;
#endif    
  // Create the matrix and right-hand side vector for the solver.
  SparseMatrix* mat = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  // Instantiate the solver itself.
  Solver* solver = create_linear_solver(matrix_solver, mat, rhs);
  

  // Solve the coarse mesh problem.
  info("Coarse mesh power iteration, %d + %d + %d + %d = %d ndof:", report_num_dofs(spaces));
  power_iteration(spaces, &wf, mkptr(pow_iter_slns), TOL_PIT_CM, mat, rhs, solver);
  // If SOLVE_ON_COARSE_MESH == true, we will store the results as the first coarse mesh solution;
  // otherwise, we will obtain this solution later by projecting the reference solution on the coarse mesh.
  if (SOLVE_ON_COARSE_MESH) 
    for_each_group(g) 
      slptr_coarse_slns[g]->copy(slptr_pow_iter_slns[g]);

  // Initialize views.
  /* for 1280x800 display */
  ScalarView view1("Neutron flux 1", new WinGeom(0, 0, 320, 400));
  ScalarView view2("Neutron flux 2", new WinGeom(330, 0, 320, 400));
  ScalarView view3("Neutron flux 3", new WinGeom(660, 0, 320, 400));
  ScalarView view4("Neutron flux 4", new WinGeom(990, 0, 320, 400));
  OrderView oview1("Mesh for group 1", new WinGeom(0, 450, 320, 500));
  OrderView oview2("Mesh for group 2", new WinGeom(330, 450, 320, 500));
  OrderView oview3("Mesh for group 3", new WinGeom(660, 450, 320, 500));
  OrderView oview4("Mesh for group 4", new WinGeom(990, 450, 320, 500));

  /* for adjacent 1280x800 and 1680x1050 displays
  ScalarView view1("Neutron flux 1", new WinGeom(0, 0, 640, 480));
  ScalarView view2("Neutron flux 2", new WinGeom(650, 0, 640, 480));
  ScalarView view3("Neutron flux 3", new WinGeom(1300, 0, 640, 480));
  ScalarView view4("Neutron flux 4", new WinGeom(1950, 0, 640, 480));
  OrderView oview1("Mesh for group 1", new WinGeom(1300, 500, 340, 500));
  OrderView oview2("Mesh for group 2", new WinGeom(1650, 500, 340, 500));
  OrderView oview3("Mesh for group 3", new WinGeom(2000, 500, 340, 500));
  OrderView oview4("Mesh for group 4", new WinGeom(2350, 500, 340, 500));
  */

  Hermes::vector<ScalarView *> sviews(&view1, &view2, &view3, &view4);
  Hermes::vector<OrderView *> oviews(&oview1, &oview2, &oview3, &oview4); 
  for_each_group(g) 
  { 
    sviews[g]->show_mesh(false);
    sviews[g]->set_3d_mode(true);
  }
  
  // DOF and CPU convergence graphs
  GnuplotGraph graph_dof("Error convergence", "NDOF", "log(error [%])");
  graph_dof.add_row("H1 error est.", "r", "-", "o");
  graph_dof.add_row("L2 error est.", "g", "-", "s");
  graph_dof.add_row("Keff error est.", "b", "-", "d");
  graph_dof.set_log_y();
  graph_dof.show_legend();
  graph_dof.show_grid();

  GnuplotGraph graph_dof_evol("Evolution of NDOF", "Adaptation step", "NDOF");
  graph_dof_evol.add_row("group 1", "r", "-", "o");
  graph_dof_evol.add_row("group 2", "g", "-", "x");
  graph_dof_evol.add_row("group 3", "b", "-", "+");
  graph_dof_evol.add_row("group 4", "m", "-", "*");
  graph_dof_evol.set_log_y();
  graph_dof_evol.set_legend_pos("bottom right");
  graph_dof_evol.show_grid();

  GnuplotGraph graph_cpu("Error convergence", "CPU time [s]", "log(error [%])");
  graph_cpu.add_row("H1 error est.", "r", "-", "o");
  graph_cpu.add_row("L2 error est.", "g", "-", "s");
  graph_cpu.add_row("Keff error est.", "b", "-", "d");
  graph_cpu.set_log_y();
  graph_cpu.show_legend();
  graph_cpu.show_grid();

  // Initialize the refinement selectors.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
  Hermes::vector<RefinementSelectors::Selector*> selectors;
  for_each_group(g) selectors.push_back(&selector);

  // Adaptivity loop:
  int as = 1; bool done = false;
  do {

    info("---- Adaptivity step %d:", as);

    // Construct globally refined meshes and setup reference spaces on them.
    Hermes::vector<Space *> ref_spaces;
    Hermes::vector<Mesh *> ref_meshes;
    for_each_group(g) 
    { 
      ref_meshes.push_back(new Mesh());
      Mesh *ref_mesh = ref_meshes.back();      
      ref_mesh->copy(spaces[g]->get_mesh());
      ref_mesh->refine_all_elements();
      
      int order_increase = 1;
      ref_spaces.push_back(spaces[g]->dup(ref_mesh, order_increase));
    }

#ifdef WITH_PETSC    
    // PETSc assembling is currently slow for larger matrices, so we switch to 
    // UMFPACK when matrices of order >8000 start to appear.
    if (Space::get_num_dofs(ref_spaces) > 8000 && matrix_solver != SOLVER_UMFPACK)
    {
      // Delete the old solver.
      delete mat;
      delete rhs;
      delete solver;
      
      // Create a new one.
      matrix_solver = SOLVER_UMFPACK;
      mat = create_matrix(matrix_solver);
      rhs = create_vector(matrix_solver);
      solver = create_linear_solver(matrix_solver, mat, rhs);
    }
#endif    

    // For the first time, project coarse mesh solutions on fine meshes to obtain 
    // a starting point for the fine mesh power iteration.
    scalar* coeff_vec = new scalar[Space::get_num_dofs(spaces)];
    if (as == 1) {
      info("Projecting initial coarse mesh solutions on fine meshes.");
      OGProjection::project_global(spaces, 
                     Hermes::vector< std::pair<WeakForm::matrix_form_val_t, WeakForm::matrix_form_ord_t> >(callback_pairs(projection_biform)), 
                     Hermes::vector< std::pair<WeakForm::vector_form_val_t, WeakForm::vector_form_ord_t> >(callback_pairs(projection_liform)),
                     mfptr_pow_iter_slns, coeff_vec, matrix_solver);
      Solution::vector_to_solutions(coeff_vec, spaces, slptr_pow_iter_slns);
    }
    delete coeff_vec;
    
    // Solve the fine mesh problem.
    info("Fine mesh power iteration, %d + %d + %d + %d = %d ndof:", report_num_dofs(ref_spaces));
    power_iteration(ref_spaces, &wf, mkptr(pow_iter_slns), TOL_PIT_RM, mat, rhs, solver);
    
    // Store the results.
    for_each_group(g) slptr_fine_slns[g]->copy(slptr_pow_iter_slns[g]);

    // Either solve on coarse mesh or project the fine mesh solution on the coarse mesh.
    if (SOLVE_ON_COARSE_MESH) {
      if (as > 1) {
        info("Coarse mesh power iteration, %d + %d + %d + %d = %d ndof:", report_num_dofs(spaces));
        power_iteration(spaces, &wf, mkptr(pow_iter_slns), TOL_PIT_CM, mat, rhs, solver);
        // Store the results.
        for_each_group(g) slptr_coarse_slns[g]->copy(slptr_pow_iter_slns[g]);
      }
    }
    else {
      scalar* coeff_vec = new scalar[Space::get_num_dofs(spaces)];
      info("Projecting fine mesh solutions on coarse meshes.");
      OGProjection::project_global(spaces,Hermes::vector< std::pair<WeakForm::matrix_form_val_t, WeakForm::matrix_form_ord_t> >(callback_pairs(projection_biform)), 
                     Hermes::vector< std::pair<WeakForm::vector_form_val_t, WeakForm::vector_form_ord_t> >(callback_pairs(projection_liform)),
                     mfptr_pow_iter_slns, coeff_vec, matrix_solver);
      Solution::vector_to_solutions(coeff_vec, spaces, slptr_coarse_slns);
      delete coeff_vec;
    }

    // Time measurement.
    cpu_time.tick();

    // View the coarse mesh solution and meshes.
    for_each_group(g) { 
      sviews[g]->show(slptr_coarse_slns[g]); 
      oviews[g]->show(spaces[g]);
    }

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Report the number of negative eigenfunction values.
    info("Num. of negative values: %d, %d, %d, %d",
         get_num_of_neg(slptr_coarse_slns[0]), get_num_of_neg(slptr_coarse_slns[1]),
         get_num_of_neg(slptr_coarse_slns[2]), get_num_of_neg(slptr_coarse_slns[3]));

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
    double energy_err_est = hp.calc_err_est(slptr_coarse_slns, slptr_fine_slns) * 100;
    double h1_err_est = error_total(error_fn_h1_axisym, norm_fn_h1_axisym, slptr_coarse_slns, slptr_fine_slns);
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
    
    // Report results.
    info("ndof_coarse: %d + %d + %d + %d = %d", report_num_dofs(spaces));

    // Millipercent eigenvalue error w.r.t. the reference value (see physical_parameters.cpp). 
  	double keff_err = 1e5*fabs(k_eff - REF_K_EFF)/REF_K_EFF;

  	info("per-group err_est_coarse (H1): %g%%, %g%%, %g%%, %g%%", report_errors(h1_errors));
  	info("per-group err_est_coarse (L2): %g%%, %g%%, %g%%, %g%%", report_errors(l2_errors));
  	info("total err_est_coarse (energy): %g%%", energy_err_est);
  	info("total err_est_coarse (H1): %g%%", h1_err_est);
  	info("total err_est_coarse (L2): %g%%", l2_err_est);
  	info("k_eff err: %g milli-percent", keff_err);

    // Add entry to DOF convergence graph.
    int ndof_coarse = Space::get_num_dofs(spaces);
    graph_dof.add_values(0, ndof_coarse, h1_err_est);
    graph_dof.add_values(1, ndof_coarse, l2_err_est);
    graph_dof.add_values(2, ndof_coarse, keff_err);

    // Add entry to CPU convergence graph.
    graph_cpu.add_values(0, cta, h1_err_est);
    graph_cpu.add_values(1, cta, l2_err_est);
    graph_cpu.add_values(2, cta, keff_err);

    for_each_group(g)
      graph_dof_evol.add_values(g, as, spaces[g]->Space::get_num_dofs());

    cpu_time.tick(HERMES_SKIP);

    // If err_est too large, adapt the mesh.
    if (energy_err_est < ERR_STOP) break;
    else {
      info("Adapting the coarse mesh.");
      done = hp.adapt(selectors, THRESHOLD, STRATEGY, MESH_REGULARITY);
      if (Space::get_num_dofs(spaces) >= NDOF_STOP) done = true;
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
  
  for_each_group(g) 
  {
    delete spaces[g]; delete meshes[g];
    delete slptr_coarse_slns[g], delete slptr_fine_slns[g]; delete slptr_pow_iter_slns[g];
  }
  
  delete mat;
  delete rhs;
  delete solver;

  graph_dof.save("conv_dof.gp");
  graph_cpu.save("conv_cpu.gp");
  graph_dof_evol.save("dof_evol.gp");

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
