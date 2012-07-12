#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::RefinementSelectors;

// This test makes sure that example 14-hcurl-adapt works correctly.

// Initial polynomial degree. NOTE: The meaning is different from
// standard continuous elements in the space H1. Here, P_INIT refers
// to the maximum poly order of the tangential component, and polynomials
// of degree P_INIT + 1 are present in element interiors. P_INIT = 0
// is for Whitney elements.
const int P_INIT = 2;
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.

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
const int STRATEGY = 1;

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

// Problem parameters.
const double MU_R   = 1.0;
const double KAPPA  = 1.0;
const double LAMBDA = 1.0;

// Bessel functions, exact solution, and weak forms.
#include "../definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("../lshape3q.mesh", &mesh);    // quadrilaterals
  //mloader.load("../lshape3t.mesh", &mesh);  // triangles

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++)  mesh.refine_all_elements();

  // Initialize boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<std::complex<double> > bc_essential(Hermes::vector<std::string>("Corner_horizontal",
    "Corner_vertical"), 0);
  EssentialBCs<std::complex<double> > bcs(&bc_essential);

  // Create an Hcurl space with default shapeset.
  HcurlSpace<std::complex<double> > space(&mesh, &bcs, P_INIT);
  int ndof = space.get_num_dofs();

  // Initialize the weak formulation.
  CustomWeakForm wf(MU_R, KAPPA);

  // Initialize coarse and reference mesh solutions.
  Solution<std::complex<double> > sln, ref_sln;

  // Initialize exact solution.
  CustomExactSolution sln_exact(&mesh);

  // Initialize refinement selector.
  HcurlProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  DiscreteProblem<std::complex<double> > dp(&wf, &space);
  
  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    // Construct globally refined reference mesh and setup reference space.
    Space<std::complex<double> >* ref_space = Space<std::complex<double> >::construct_refined_space(&space);
    dp.set_space(ref_space);
    int ndof_ref = ref_space->get_num_dofs();

    // Initial coefficient vector for the Newton's method.
    std::complex<double>* coeff_vec = new std::complex<double>[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(std::complex<double>));

    // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
    Hermes::Hermes2D::NewtonSolver<std::complex<double> > newton(&dp);

    try{
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.printMsg();
    }
    Hermes::Hermes2D::Solution<std::complex<double> >::vector_to_solution(newton.get_sln_vector(), ref_space, &ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<std::complex<double> > ogProjection;
    ogProjection.project_global(&space, &ref_sln, &sln);

    // Calculate element errors and total error estimate.
    Adapt<std::complex<double> >* adaptivity = new Adapt<std::complex<double> >(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // Calculate exact error.
    bool solutions_for_adapt = false;
    double err_exact_rel = adaptivity->calc_err_exact(&sln, &sln_exact, solutions_for_adapt) * 100;

    // If err_est_rel too large, adapt the mesh.
    if(err_est_rel < ERR_STOP) done = true;
    else
    {
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      // Increase the counter of performed adaptivity steps.
      if(done == false)  as++;
    }
    if(space.get_num_dofs() >= NDOF_STOP) done = true;

    // Clean up.
    delete [] coeff_vec;
    delete adaptivity;
    if(done == false)
      delete ref_space->get_mesh();
    delete ref_space;
  }
  while (done == false);

  ndof = space.get_num_dofs();

  if(ndof == 382) // Tested value as of 12 Jul 2011.
  {
    printf("Success!\n");
    return 0;
  }
  else {
    printf("Failure!\n");
    return -1;
  }
}