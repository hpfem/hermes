#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::RefinementSelectors;

//  This is a test of adaptivity.
//
//  PDE: -Laplace u - f = 0.
//
//  Known exact solution, see class CustomExactSolution in definitions.cpp.
//
//  Domain: square domain (0, pi) x (0, pi), mesh file square_quad.mesh.
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

int P_INIT = 1;                                   // Initial polynomial degree of all mesh elements.
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
                                                  // See User Documentation.
const int MESH_REGULARITY = -1;                   // Maximum allowed level of hanging nodes:
                                                  // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                                  // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                                  // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                                  // Note that regular meshes are not supported, this is due to
                                                  // their notoriously bad performance.
const double CONV_EXP = 1.0;                      // Default value is 1.0. This parameter influences the selection of
                                                  // cancidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 1e-4;                     // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;                      // Adaptivity process stops when the number of degrees of freedom grows
                                                  // over this limit. This is to prevent h-adaptivity to go on forever.
Hermes::MatrixSolverType matrix_solver = Hermes::SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("square_quad.mesh", &mesh);

  // Avoid zero ndof situation.
  if(P_INIT == 1) {
    if(is_hp(CAND_LIST)) P_INIT++;
    else mesh.refine_element_id(0, 0);
  }

  // Define exact solution.
  CustomExactSolution exact_sln(&mesh);

  // Initialize the weak formulation.
  Hermes1DFunction<double> lambda(1.0);
  CustomFunction f;
  DefaultWeakFormPoisson<double> wf(Hermes::HERMES_ANY, &lambda, &f);

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc_essential("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<double> space(&mesh, &bcs, P_INIT);

  // Initialize approximate solution.
  Solution<double> sln;

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Assemble the discrete problem.
  DiscreteProblem<double> dp(&wf, &space);

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    // Construct globally refined reference mesh and setup reference space.
    Space<double>* ref_space = Space<double>::construct_refined_space(&space);
    int ndof_ref = ref_space->get_num_dofs();

    dp.set_space(ref_space);

    // Initial coefficient vector for the Newton's method.
    double* coeff_vec = new double[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(double));

    NewtonSolver<double> newton(&dp);
    newton.set_verbose_output(false);

    Solution<double> ref_sln;
    try{
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.printMsg();
    }
    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, &ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<double> ogProjection;
    ogProjection.project_global(&space, &ref_sln, &sln);

    // Calculate element errors and total error estimate.
    Adapt<double> adaptivity(&space);
    double err_est_rel = adaptivity.calc_err_est(&sln, &ref_sln) * 100;

    // Calculate exact error.
    double err_exact_rel = Global<double>::calc_rel_error(&sln, &exact_sln, HERMES_H1_NORM) * 100;

    // Report results.

    // If err_est too large, adapt the mesh. The NDOF test must be here, so that the solution may be visualized
    // after ending due to this criterion.
    if(err_exact_rel < ERR_STOP || space.get_num_dofs() >= NDOF_STOP)
      done = true;
    else
      done = adaptivity.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

    // Increase the counter of adaptivity steps.
    if(done == false)
      as++;

    // Clean up.
    delete [] coeff_vec;

    if(done == false)
      delete ref_space->get_mesh();
    delete ref_space;
  }
  while (done == false);

  if(space.get_mesh()->get_num_active_elements() == 1)
  {
    return 0;
  }
  else {
    return -1;
  }
}