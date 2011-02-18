#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

/** \addtogroup t_bench_stabilized-advection-diffusion Benchmarks/stabilized-advection-diffusion
 *  \{
 *  \brief This test makes sure that the benchmark "stabilized-advection-diffusion" works correctly.
 *
 *  \section s_params Parameters
 *   - P_INIT = 0
 *   - INIT_REF_NUM = 1
 *   - INIT_BDY_REF_NUM = 2
 *   - ORDER_INCREASE = 1
 *   - THRESHOLD=0.2
 *   - STRATEGY=0
 *   - CAND_LIST=H2D_HP_ANISO;
 *   - MESH_REGULARITY=-1
 *   - CONV_EXP=1.0
 *   - ERR_STOP=5.0
 *   - NDOF_STOP=60000
 *   - matrix_solver = SOLVER_UMFPACK
 *   - method = DG
 *
 *  \section s_res Results
 *   - DOFs: 36
 *   - Adaptivity steps: 9 
 */

enum GalerkinMethod
{
  CG,
  CG_STAB_SUPG,    // assumes H2D_SECOND_DERIVATIVES_ENABLED, linear elements
  CG_STAB_GLS,     // assumes H2D_SECOND_DERIVATIVES_ENABLED, linear or quadratic elements
  CG_STAB_SGS,     // assumes H2D_SECOND_DERIVATIVES_ENABLED, linear elements
  CG_STAB_SGS_ALT, // assumes H2D_SECOND_DERIVATIVES_ENABLED, linear elements
  DG
};

const std::string method_names[6] = 
{
  "unstabilized continuous Galerkin",
  "streamline upwind Petrov-Galerkin",
  "Galerkin least squares",
  "subgrid scale stabilized Galerkin",
  "alternative subgrid scale stabilized Galerkin",
  "discontinuous Galerkin"
};

const int P_INIT = 0;                             // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements. 
const int INIT_BDY_REF_NUM = 2;                   // Number of initial refinements towards boundary. If INIT_BDY_REF_NUM == 0, 
                                                  // the first solution will be performed on a mesh (INIT_REF_NUM + 1) times 
                                                  // globally refined.
const int ORDER_INCREASE = 1;                     // Order increase for the refined space. If no change of order is allowed 
                                                  // (not even for computing the reference solution), set to 0.
const double THRESHOLD = 0.2;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
                                                  // STRATEGY = -1... do not refine.
                                                  // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                                  //   error is processed. If more elements have similar errors, refine
                                                  //   all to keep the mesh symmetric.
                                                  // STRATEGY = 1 ... refine all elements whose error is larger
                                                  //   than THRESHOLD times maximum element error.
                                                  // STRATEGY = 2 ... refine all elements whose error is larger
                                                  //   than THRESHOLD.
                                                  // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates. Possible values are
                                                  // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                                  // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 5.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // exact and the coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.

GalerkinMethod method = DG;
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK. 

// Problem parameters.
const double EPSILON = 0.01;                      // Diffusivity.
const double B1 = 0., B2 = 1.;                    // Advection direction, div(B) = 0.

// Boundary markers.
const int NONZERO_DIRICHLET = 1;
const int BOUNDARY_LAYER = 2;

// Boundary condition types.
BCType bc_types(int marker)
{
  if (method != DG)
    return BC_ESSENTIAL;
  else
    return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
template<typename Real, typename Scalar>
Scalar essential_bc_values(int ess_bdy_marker, Real x, Real y)
{
    if (ess_bdy_marker == NONZERO_DIRICHLET) return sin(M_PI*x);
    else return 0;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Exact solution:

static double exact_sln(double x, double y, double& dx, double& dy)
{
  double D = sqrt(1.+4.*sqr(EPSILON*M_PI));
  double m1 = (1.-D) / (2.*EPSILON);
  double m2 = (1.+D) / (2.*EPSILON);
  
  double e12y = exp(m1+m2*y);
  double e21y = exp(m2+m1*y);
  double ie1me2 = 1. / (exp(m1) - exp(m2));
  
  dx = ie1me2 * (e12y - e21y) * M_PI * cos(M_PI*x);
  dy = ie1me2 * (e12y*m2 - e21y*m1) * sin(M_PI*x);
  
  return ie1me2 * (e12y - e21y) * sin(M_PI*x);
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square_quad.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();
  mesh.refine_towards_boundary(BOUNDARY_LAYER, INIT_BDY_REF_NUM);
  //mesh.refine_towards_boundary(NONZERO_DIRICHLET, INIT_BDY_REF_NUM/2);

  // Create a space and refinement selector appropriate for the selected discretization method.
  Space *space;
  ProjBasedSelector *selector;
  ProjNormType norm;
  if (method != DG)
  {
    space = new H1Space(&mesh, bc_types, essential_bc_values, P_INIT);
    selector = new L2ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
    norm = HERMES_L2_NORM;  // WARNING: In order to compare the errors with DG, L2 norm should be here.
  }
  else
  {
    space = new L2Space(&mesh, bc_types, NULL, Ord2(P_INIT));
    selector = new L2ProjBasedSelector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);
    norm = HERMES_L2_NORM;
    // Disable weighting of refinement candidates.
    selector->set_error_weights(1, 1, 1);
  }

  // Initialize the weak formulation.
  info("Discretization method: %s", method_names[method].c_str());
    
  if (method != CG && method != DG)
  {
    if (STRATEGY > -1 && CAND_LIST != H2D_H_ISO && CAND_LIST != H2D_H_ANISO)
      error("The %s method may be used only with h-refinement.", method_names[method].c_str());
  
    int eff_order = (STRATEGY == -1) ? P_INIT : P_INIT + ORDER_INCREASE;
    
    if (method != CG_STAB_GLS)
    {
      if (eff_order > 1)
        error("The %s method may be used only with at most 1st order elements.", method_names[method].c_str());
    }
    else
    {
      if (eff_order > 2)
        error("The %s method may be used only with at most 2nd order elements.", method_names[method].c_str());
    }
  }
  
  WeakForm wf;
  switch(method)
  {
    case CG:
      wf.add_matrix_form(callback(cg_biform));
      break;
    case CG_STAB_SUPG:
      wf.add_matrix_form(callback(cg_biform));
      wf.add_matrix_form(callback(stabilization_biform_supg));
      break; 
    case CG_STAB_GLS:
      wf.add_matrix_form(callback(cg_biform));
      wf.add_matrix_form(callback(stabilization_biform_gls));
      break;
    case CG_STAB_SGS:
      wf.add_matrix_form(callback(cg_biform));
      wf.add_matrix_form(callback(stabilization_biform_sgs));
      break;
    case CG_STAB_SGS_ALT:
      wf.add_matrix_form(callback(cg_biform));
      wf.add_matrix_form(callback(stabilization_biform_sgs_alt));
      break;
    case DG:
      wf.add_matrix_form(callback(dg_volumetric_biform_advection));
      wf.add_matrix_form(callback(dg_volumetric_biform_diffusion));
      wf.add_matrix_form_surf(callback(dg_interface_biform_advection), H2D_DG_INNER_EDGE);
      wf.add_matrix_form_surf(callback(dg_interface_biform_diffusion), H2D_DG_INNER_EDGE);
      wf.add_matrix_form_surf(callback(dg_boundary_biform_advection));
      wf.add_matrix_form_surf(callback(dg_boundary_biform_diffusion));
      wf.add_vector_form_surf(callback(dg_boundary_liform_advection));
      wf.add_vector_form_surf(callback(dg_boundary_liform_diffusion));
      break;
  }
  
  // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln;
  
  // Set exact solution.
  ExactSolution exact(&mesh, exact_sln);
  
  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu, graph_dof_exact, graph_cpu_exact;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Setup data structures for solving the discrete algebraic problem.
  SparseMatrix* matrix = create_matrix(matrix_solver);
  Vector* rhs = create_vector(matrix_solver);
  Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  Space* actual_sln_space;
  do
  {
    info("---- Adaptivity step %d:", as);

    if (STRATEGY == -1)
      actual_sln_space = space;
    else
      // Construct globally refined reference mesh and setup reference space.
      actual_sln_space = construct_refined_space(space, ORDER_INCREASE);

    // Assemble the reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, actual_sln_space, is_linear);
    dp->assemble(matrix, rhs);
    
    // Solve the linear system of the reference problem. 
    // If successful, obtain the solution.
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), actual_sln_space, &ref_sln);
    else error ("Matrix solver failed.\n");
    
    // Instantiate adaptivity and error calculation driver. Space is used only for adaptivity, it is ignored when 
    // STRATEGY == -1 and only the exact error is calculated by this object.
    Adapt* adaptivity = new Adapt(space, norm);
    
    // Calculate and report exact error.
    bool solutions_for_adapt = false;
    double err_exact_rel = adaptivity->calc_err_exact(&ref_sln, &exact, solutions_for_adapt, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;
    info("ndof_fine: %d, err_exact_rel: %g%%", Space::get_num_dofs(actual_sln_space), err_exact_rel);

    // Time measurement.
    cpu_time.tick();

    if (STRATEGY == -1) done = true;  // Do not adapt.
    else
    {  
      // Project the fine mesh solution onto the coarse mesh.
      info("Projecting reference solution on coarse mesh.");
      OGProjection::project_global(space, &ref_sln, &sln, matrix_solver, norm); 

      // Calculate element errors and total error estimate.
      info("Calculating error estimate."); 
      bool solutions_for_adapt = true;
      double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt, HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

      // Report results.
      info("ndof_coarse: %d, err_est_rel: %g%%", Space::get_num_dofs(space), err_est_rel);

      // Time measurement.
      cpu_time.tick();
      
      // Add entry to DOF and CPU convergence graphs.
      graph_dof.add_values(Space::get_num_dofs(space), err_est_rel);
      graph_dof.save("conv_dof_est.dat");
      graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
      graph_cpu.save("conv_cpu_est.dat");
      graph_dof_exact.add_values(Space::get_num_dofs(space), err_exact_rel);
      graph_dof_exact.save("conv_dof_exact.dat");
      graph_cpu_exact.add_values(cpu_time.accumulated(), err_exact_rel);
      graph_cpu_exact.save("conv_cpu_exact.dat");

      // Skip graphing time.
      cpu_time.tick(HERMES_SKIP);
      
      // If err_est too large, adapt the mesh.
      if (err_exact_rel < ERR_STOP) done = true;
      else 
      {
        info("Adapting coarse mesh.");
        done = adaptivity->adapt(selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
  
        // Increase the counter of performed adaptivity steps.
        if (done == false)  as++;
      }
      if (Space::get_num_dofs(space) >= NDOF_STOP) done = true;
      
      if(done == false) 
      {
        delete actual_sln_space->get_mesh();
        delete actual_sln_space;
      }
    }
    
    // Clean up.
    delete adaptivity;
    delete dp;
  }
  while (done == false);
  
  int ndof = Space::get_num_dofs(space);

  if (space != actual_sln_space) 
  {
    delete space;
    delete actual_sln_space->get_mesh();
  }
  delete actual_sln_space;
  delete solver;
  delete matrix;
  delete rhs;
  delete selector;
  
  verbose("Total running time: %g s", cpu_time.accumulated());

  int n_dof_allowed = 40;
  printf("n_dof_actual = %d\n", ndof);
  printf("n_dof_allowed = %d\n", n_dof_allowed);
  if (ndof <= n_dof_allowed) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
