#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

using namespace Hermes::Hermes2D::RefinementSelectors;

// This test makes sure that example 13-complex-adapt works correctly.

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;
// Initial polynomial degree of all mesh elements.
const int P_INIT = 1;
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies (see below).
const double THRESHOLD = 0.3;
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
// More adaptive strategies can be created in adapt_ortho_h1.cpp.
const int STRATEGY = 0;
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
// H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
// See User Documentation for details.
const CandList CAND_LIST = H2D_HP_ANISO;
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;
// Default value is 1.0. This parameter influences the selection of
// cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double CONV_EXP = 1.0;
// Stopping criterion for adaptivity (rel. error tolerance between the
// reference mesh and coarse mesh solution in percent).
const double ERR_STOP = 10.0;
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;
// Name of the iterative method employed by AztecOO (ignored
// by the other solvers). 
// Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* iterative_method = "bicgstab";
// Name of the preconditioner employed by AztecOO (ignored by
// the other solvers).
// Possibilities: none, jacobi, neumann, least-squares, or a
// preconditioner from IFPACK (see solver/aztecoo.h).
const char* preconditioner = "least-squares";
// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;

// Problem parameters.
const double MU_0 = 4.0*M_PI*1e-7;
const double MU_IRON = 1e3 * MU_0;
const double GAMMA_IRON = 6e6;
const double J_EXT = 1e6;
const double FREQ = 5e3;
const double OMEGA = 2 * M_PI * FREQ;

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<std::complex<double> > bc_essential("Dirichlet", std::complex<double>(0.0, 0.0));
  EssentialBCs<std::complex<double> > bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<std::complex<double> > space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakForm wf("Air", MU_0, "Iron", MU_IRON, GAMMA_IRON, 
    "Wire", MU_0, std::complex<double>(J_EXT, 0.0), OMEGA);

  // Initialize coarse and reference mesh solution.
  Solution<std::complex<double> > sln, ref_sln;

  // Initialize refinement selector.
  H1ProjBasedSelector<std::complex<double> > selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  Space<std::complex<double> >* ref_space = Space<std::complex<double> >::construct_refined_space(&space);
  
  DiscreteProblem<std::complex<double> > dp(&wf, ref_space);
  dp.set_adaptivity_cache();
    
  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::NewtonSolver<std::complex<double> > newton(&dp, matrix_solver_type);

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    ref_space = Space<std::complex<double> >::construct_refined_space(&space);
    dp.set_spaces(ref_space);
    int ndof_ref = ref_space->get_num_dofs();

    // Initialize reference problem.
    info("Solving on reference mesh.");

    // Time measurement.
    cpu_time.tick();

    // Initial coefficient vector for the Newton's method.  
    std::complex<double>* coeff_vec = new std::complex<double>[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(std::complex<double>));

    // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
    // For iterative solver.
    if (matrix_solver_type == SOLVER_AZTECOO)
    {
      newton.set_iterative_method(iterative_method);
      newton.set_preconditioner(preconditioner);
    }
    try{
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.printMsg();
      error("Newton's iteration failed.");
    }
    Hermes::Hermes2D::Solution<std::complex<double> >::vector_to_solution(newton.get_sln_vector(), ref_space, &ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection<std::complex<double> >::project_global(&space, &ref_sln, &sln, matrix_solver_type); 

    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt<std::complex<double> >* adaptivity = new Adapt<std::complex<double> >(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
      space.get_num_dofs(), ref_space->get_num_dofs(), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (space.get_num_dofs() >= NDOF_STOP) done = true;

    // Clean up.
    delete [] coeff_vec;
    delete adaptivity;

    // Increase counter.
    as++;
  }
  while (done == false);

  verbose("Total running time: %g s", cpu_time.accumulated());

  ndof = space.get_num_dofs();

  if (ndof == 86) // Tested value as of 12 Jul 2011.
  {
    printf("Success!\n");
    return TEST_SUCCESS;
  }
  else 
  {
    printf("Failure!\n");
    return TEST_FAILURE;
  }
}
