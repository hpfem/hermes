#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

using namespace Hermes::Hermes2D::RefinementSelectors;

// This test makes sure that example 13-complex-adapt works correctly.

const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
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
                                                  // See User Documentation for details.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1.0;                      // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "least-squares";     // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers).
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  //  preconditioner from IFPACK (see solver/aztecoo.h)
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Problem parameters.
const double MU_0 = 4.0*M_PI*1e-7;
const double MU_IRON = 1e3 * MU_0;
const double GAMMA_IRON = 6e6;
const double J_EXT = 1e6;
const double FREQ = 5e3;
const double OMEGA = 2 * M_PI * FREQ;

int main(int argc, char* argv[])
{
  // Instantiate a class with global functions.
  Global<std::complex<double> > hermes2d;

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
 Hermes::Hermes2D::DefaultEssentialBCConst<std::complex<double> > bc_essential("Dirichlet", std::complex<double>(0.0, 0.0));
  EssentialBCs<std::complex<double> > bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<std::complex<double> > space(&mesh, &bcs, P_INIT);
  int ndof = Space<std::complex<double> >::get_num_dofs(&space);
  info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakForm wf("Air", MU_0, "Iron", MU_IRON, GAMMA_IRON, 
                    "Wire", MU_0, std::complex<double>(J_EXT, 0.0), OMEGA);

  // Initialize coarse and reference mesh solution.
  Solution<std::complex<double> > sln, ref_sln;

  // Initialize refinement selector.
  H1ProjBasedSelector<std::complex<double> > selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Initialize the matrix solver.
  SparseMatrix<std::complex<double> >* matrix = create_matrix<std::complex<double> >(matrix_solver);
  Vector<std::complex<double> >* rhs = create_vector<std::complex<double> >(matrix_solver);
  Solver<std::complex<double> >* solver = create_linear_solver(matrix_solver, matrix, rhs);
  
  if (matrix_solver == SOLVER_AZTECOO) {
    dynamic_cast<AztecOOSolver<std::complex<double> >*>(solver)->set_solver(iterative_method);
    dynamic_cast<AztecOOSolver<std::complex<double> >*>(solver)->set_precond(preconditioner);
    // Using default iteration parameters (see solver/aztecoo.h).
  }
  
  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    info("---- Adaptivity step %d:", as);

    // Construct globally refined reference mesh and setup reference space.
    Space<std::complex<double> >* ref_space = Space<std::complex<double> >::construct_refined_space(&space);
    int ndof_ref = Space<std::complex<double> >::get_num_dofs(ref_space);

    // Initialize matrix solver.
    SparseMatrix<std::complex<double> >* matrix = create_matrix<std::complex<double> >(matrix_solver);
    Vector<std::complex<double> >* rhs = create_vector<std::complex<double> >(matrix_solver);
    Solver<std::complex<double> >* solver = create_linear_solver<std::complex<double> >(matrix_solver, matrix, rhs);

    // Initialize reference problem.
    info("Solving on reference mesh.");
    DiscreteProblem<std::complex<double> > dp(&wf, ref_space);

    // Time measurement.
    cpu_time.tick();

    // Initial coefficient vector for the Newton's method.  
    std::complex<double>* coeff_vec = new std::complex<double>[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(std::complex<double>));

    // Perform Newton's iteration.
    if (!hermes2d.solve_newton(coeff_vec, &dp, solver, matrix, rhs)) error("Newton's iteration failed.");

    // Translate the resulting coefficient vector into the Solution sln.
    Solution<std::complex<double> >::vector_to_solution(coeff_vec, ref_space, &ref_sln);
  
    // Project the fine mesh solution onto the coarse mesh.
    info("Projecting reference solution on coarse mesh.");
    OGProjection<std::complex<double> >::project_global(&space, &ref_sln, &sln, matrix_solver); 
   
    // Calculate element errors and total error estimate.
    info("Calculating error estimate."); 
    Adapt<std::complex<double> >* adaptivity = new Adapt<std::complex<double> >(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // Report results.
    info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
      Space<std::complex<double> >::get_num_dofs(&space), Space<std::complex<double> >::get_num_dofs(ref_space), err_est_rel);

    // Time measurement.
    cpu_time.tick();

    // If err_est too large, adapt the mesh.
    if (err_est_rel < ERR_STOP) done = true;
    else 
    {
      info("Adapting coarse mesh.");
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
    }
    if (Space<std::complex<double> >::get_num_dofs(&space) >= NDOF_STOP) done = true;

    // Clean up.
    delete [] coeff_vec;
    delete adaptivity;
    if (done == false)
      delete ref_space->get_mesh();
    delete ref_space;
    
    // Increase counter.
    as++;
  }
  while (done == false);
  
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  // Clean up.
  delete solver;
  delete matrix;
  delete rhs;

  ndof = Space<std::complex<double> >::get_num_dofs(&space);

  printf("ndof allowed = %d\n", 360);
  printf("ndof actual = %d\n", ndof);
  if (ndof < 650) {      // ndofs was 352 at the time this test was created
    printf("Success!\n");
    return TEST_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return TEST_FAILURE;
  }
}
