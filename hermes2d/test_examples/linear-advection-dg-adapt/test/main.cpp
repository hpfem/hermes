#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "../definitions.h"

//  This example solves a linear advection equation using Dicontinuous Galerkin (DG) method.
//  It is intended to show how evalutation of surface matrix forms that take basis functions defined
//  on different elements work. It is the same example as linear-advection-dg, but with automatic adaptivity.
//
//  PDE: \nabla \cdot (\Beta u) = 0, where \Beta = (-x_2, x_1) / |x| represents a circular counterclockwise flow field.
//
//  Domain: Square (0, 1) x (0, 1).
//
//  BC:    Dirichlet, u = 1 where \Beta(x) \cdot n(x) < 0, that is on[0,0.5] x {0}, and g = 0 anywhere else.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF = 1;
// Initial polynomial degrees of mesh elements in vertical and horizontal directions.
const int P_INIT = 0;
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.2;
// Adaptive strategy:
// STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
//   error is processed. If more elements have similar errors, refine
//   all to keep the mesh symmetric.
// STRATEGY = 1 ... refine all elements whose error is larger
//   than THRESHOLD times maximum element error.
// STRATEGY = 2 ... refine all elements whose error is larger
//   than THRESHOLD.
const int STRATEGY = 1;
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
// H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO;
// Maximum allowed level of hanging nodes:
// MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
// MESH_REGULARITY = 1 ... at most one-level hanging nodes,
// MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
// Note that regular meshes are not supported, this is due to
// their notoriously bad performance.
const int MESH_REGULARITY = -1;
// Stopping criterion for adaptivity.
const double ERR_STOP = 5.0;
// This parameter influences the selection of
// candidates in hp-adaptivity. Default value is 1.0.
const double CONV_EXP = 1.0;
// Adaptivity process stops when the number of degrees of freedom grows
// over this limit. This is to prevent h-adaptivity to go on forever.
const int NDOF_STOP = 60000;
// Name of the iterative method employed by AztecOO (ignored by the other solvers).
// Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* iterative_method = "bicgstab";
// Name of the preconditioner employed by AztecOO (ignored by the other solvers).
// Possibilities: none, jacobi, neumann, least-squares, or a
// preconditioner from IFPACK (see solver/aztecoo.h).
const char* preconditioner = "jacobi";

int main(int argc, char* args[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("../square.mesh", &mesh);

  // Perform initial mesh refinement.
  for (int i=0; i<INIT_REF; i++) 
    mesh.refine_all_elements();

  // Create an L2 space.
  L2Space<double> space(&mesh, P_INIT);

  // Initialize refinement selector.
  L2ProjBasedSelector<double> selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Disable weighting of refinement candidates.
  selector.set_error_weights(1, 1, 1);

  Solution<double> sln;
  Solution<double> ref_sln;

  // Initialize the weak formulation.
  CustomWeakForm wf("Bdy_bottom_left", &mesh);

  // Initialize the FE problem.
  DiscreteProblemLinear<double> dp(&wf, &space);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(&dp);

  int as = 1; bool done = false;
  do
  {
    // Construct globally refined reference mesh
    // and setup reference space.
    Space<double>* ref_space = Space<double>::construct_refined_space(&space);

    dp.set_space(ref_space);

    // Solve the linear system. If successful, obtain the solution.
    try
    {
      linear_solver.solve();
      Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), ref_space, &ref_sln);
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }
    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<double> ogProjection;
    ogProjection.project_global(&space, &ref_sln, &sln, HERMES_L2_NORM);

    // Calculate element errors and total error estimate.
    Adapt<double>* adaptivity = new Adapt<double>(&space);
    double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln) * 100;

    // If err_est_rel too large, adapt the mesh.
    if(err_est_rel < ERR_STOP) done = true;
    else
    {
      done = adaptivity->adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);

      if(Space<double>::get_num_dofs(&space) >= NDOF_STOP)
      {
        done = true;
        break;
      }
    }

    // Clean up.
    delete adaptivity;
    if(done == false)
      delete ref_space->get_mesh();
    delete ref_space;

    as++;
  }
  while (done == false);

  if(done)
  {
    printf("Success!\n");
    return 0;
  }
  else
  {
    printf("Failure!\n");
    return -1;
  }
}