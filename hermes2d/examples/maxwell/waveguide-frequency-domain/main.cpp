#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace RefinementSelectors;

// This example solves adaptively the electric field in a simplified microwave oven.
// The waves are generated using a harmonic surface current on the right-most edge.
// (Such small cavity is present in every microwave oven). There is a circular
// load located in the middle of the main cavity, defined through a different
// permittivity -- see function in_load(...). One can either use a mesh that is
// aligned to the load via curvilinear elements (ALIGN_MESH = true), or an unaligned
// mesh (ALIGN_MESH = false). Convergence graphs are saved both wrt. the dof number
// and cpu time.
//
// PDE: time-harmonic Maxwell's equations;
//      there is circular load in the middle of the large cavity, whose permittivity
//      is different from the rest of the domain.
//
// Domain: square cavity with another small square cavity attached from outside
//         on the right.
//
// Meshes: you can either use "oven_load_circle.mesh" containing curved elements
//         aligned with the circular load, or "oven_load_square.mesh" which is not
//         aligned.
//
// BC: perfect conductor on the boundary except for the right-most edge of the small
//     cavity, where a harmonic surface current is prescribed
//
// The following parameters can be changed:

const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int P_INIT = 2;                             // Initial polynomial degree. NOTE: The meaning is different from
                                                  // standard continuous elements in the space H1. Here, P_INIT refers
                                                  // to the maximum poly order of the tangential component, and polynomials
                                                  // of degree P_INIT + 1 are present in element interiors. P_INIT = 0
                                                  // is for Whitney elements.
const bool ALIGN_MESH = true;                     // if ALIGN_MESH == true, curvilinear elements aligned with the
                                                  // circular load are used, otherwise one uses a non-aligned mesh.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of the adapt(...) function and
                                                  // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 0;                           // Adaptive strategy:
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
                                                  // See the User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 2.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double e_0   = 8.8541878176 * 1e-12;
const double mu_0   = 1.256 * 1e-6;
const double e_r = 1.0;
const double mu_r = 1.0;
const double rho = 3820.0;
const double Cp = 7.531000;
const double freq = 1.0*2450000000.0;
const double omega = 2 * M_PI * freq;
const double c = 1 / sqrt(e_0 * mu_0);
const double kappa  = 2 * M_PI * freq * sqrt(e_0 * mu_0);
const double J = 0.0000033333;

//  Boundary markers.
const std::string BDY_PERFECT_CONDUCTOR = "2";
const std::string BDY_CURRENT = "1";

/* WEAK FORM FOR ERROR CALCULATION - TO BE USED IN ADAPTIVITY.
// error calculation
template<typename Real, typename Scalar>
Scalar hcurl_form_kappa(int n, double *wt, Func<Scalar> *u_ext[], Func<Scalar> *u, Func<Scalar> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_curl_e_curl_f<Scalar, Scalar>(n, wt, u, v) + sqr(kappa) * int_e_f<Scalar, Scalar>(n, wt, u, v);
}
*/

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  if (ALIGN_MESH) mloader.load("oven_load_circle.mesh", &mesh);
  else mloader.load("oven_load_square.mesh", &mesh);
  
  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Initialize the weak formulation.
  CustomWeakForm wf(e_0, mu_0, mu_r, kappa, omega, J, ALIGN_MESH);

  // Initialize boundary conditions
  DefaultEssentialBCConst bc_essential(BDY_PERFECT_CONDUCTOR, std::complex<double>(0.0, 0.0));
  EssentialBCs bcs(&bc_essential);

  // Create an Hcurl space with default shapeset.
  HcurlSpace space(&mesh, &bcs, P_INIT);
  int ndof = Space::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize coarse and reference mesh solution.
  Solution sln, ref_sln;

  // Initialize refinements selector.
  HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize views.
  VectorView eview("Electric field", new WinGeom(0, 0, 580, 400));
  OrderView  oview("Polynomial orders", new WinGeom(590, 0, 550, 400));
  
  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;
  
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Adaptivity loop:
  int as = 1; 
  bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space* ref_space = Space::construct_refined_space(&space);

    // Initialize matrix solver.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Assemble the reference problem.
    info("Solving on reference mesh.");
    bool is_linear = true;
    DiscreteProblem* dp = new DiscreteProblem(&wf, ref_space, is_linear);
    dp->assemble(matrix, rhs);

    // Time measurement.
    cpu_time.tick();
    
    // Solve the linear system of the reference problem. 
    // If successful, obtain the solution.
    if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), ref_space, &ref_sln);
    else error ("Matrix solver failed.\n");

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection::project_global(&space, &ref_sln, &sln, matrix_solver); 

    // Time measurement.
    cpu_time.tick();
   
    // Show real part of the solution.
    AbsFilter abs(&sln);
    eview.set_min_max_range(0, 4e3);
    eview.show(&abs);
    oview.show(&space);

    // Skip visualization time.
    cpu_time.tick(HERMES_SKIP);

    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt* adaptivity = new Adapt(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
      Space::get_num_dofs(&space), Space::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // Add entry to DOF and CPU convergence graphs.
    graph_dof.add_values(Space::get_num_dofs(&space), err_est_rel);
    graph_dof.save("conv_dof_est.dat");
    graph_cpu.add_values(cpu_time.accumulated(), err_est_rel);
    graph_cpu.save("conv_cpu_est.dat");

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      
      // Increase the counter of performed adaptivity steps.
      if (done == false)  as++;
    }
    if (Space::get_num_dofs(&space) >= NDOF_STOP) done = true;

    // Clean up.
    delete solver;
    delete matrix;
    delete rhs;
    delete adaptivity;
    if(done == false) delete ref_space->get_mesh();
    delete ref_space;
    delete dp;
    
  }
  while (done == false);
  
  verbose("Total running time: %g s", cpu_time.accumulated());

  // Show the reference solution - the final result.
  eview.set_title("Fine mesh solution");
  eview.show(&ref_sln);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

