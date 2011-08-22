#define HERMES_REPORT_INFO
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// This test makes sure that subdomains work correctly.

// Uniform polynomial degree of mesh elements.
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;
// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;  

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh_whole_domain, mesh_bottom_left_corner, mesh_supplement;
  Hermes::vector<Mesh*> mesh_vector (&mesh_whole_domain, &mesh_bottom_left_corner, &mesh_supplement);
  MeshReaderH2DXML mloader;
  mloader.load("subdomains.xml", mesh_vector);

  // Perform initial mesh refinements (optional).
  for(int i = 0; i < INIT_REF_NUM; i++)
    for(unsigned int meshes_i = 0; meshes_i < mesh_vector.size(); meshes_i++)
      mesh_vector[meshes_i]->refine_all_elements();

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf;

  // Initialize essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential_whole_domain(Hermes::vector<std::string>("Bottom Left", "Bottom Right", "Top Left", "Top Right"), 0.0);
  EssentialBCs<double> bcs_whole_domain(&bc_essential_whole_domain);

  DefaultEssentialBCConst<double> bc_essential_bottom_left_corner(Hermes::vector<std::string>("Bottom Left", "Vertical Bottom", "Horizontal Left"), 0.0);
  EssentialBCs<double> bcs_bottom_left_corner(&bc_essential_bottom_left_corner);

  DefaultEssentialBCConst<double> bc_essential_supplement(Hermes::vector<std::string>("Bottom Right", "Top Right", "Top Left", "Horizontal Left", "Vertical Bottom"), 0.0);
  EssentialBCs<double> bcs_supplement(&bc_essential_supplement);

  // Create H1 spaces with default shapeset.
  H1Space<double> space_whole_domain(&mesh_whole_domain, &bcs_whole_domain, P_INIT);
  int ndof = space_whole_domain.get_num_dofs();
  info("Space whole domain ndof = %d", ndof);

  H1Space<double> space_bottom_left_corner(&mesh_bottom_left_corner, &bcs_bottom_left_corner, P_INIT);
  ndof = space_bottom_left_corner.get_num_dofs();
  info("Space bottom left corner ndof = %d", ndof);

  H1Space<double> space_supplement(&mesh_supplement, &bcs_supplement, P_INIT);
  ndof = space_supplement.get_num_dofs();
  info("Space supplement ndof = %d", ndof);

  // Initialize the FE problem.
  Hermes::Hermes2D::DiscreteProblem<double> dp(&wf, Hermes::vector<Space<double>*>(&space_whole_domain, &space_bottom_left_corner, &space_supplement));

  // Initial coefficient vector for the Newton's method.  
  double* coeff_vec = new double[Space<double>::get_num_dofs(Hermes::vector<Space<double>*>(&space_whole_domain, &space_bottom_left_corner, &space_supplement))];
  memset(coeff_vec, 0, ndof*sizeof(double));

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::Solution<double> sln_whole_domain, sln_bottom_left_corner, sln_supplement;
  Hermes::Hermes2D::NewtonSolver<double> newton(&dp, matrix_solver_type);
  if (!newton.solve(coeff_vec)) 
    error("Newton's iteration failed.");
  else
    Hermes::Hermes2D::Solution<double>::vector_to_solutions(newton.get_sln_vector(), Hermes::vector<Space<double>*>(&space_whole_domain, &space_bottom_left_corner, &space_supplement), Hermes::vector<Solution<double>*>(&sln_whole_domain, &sln_bottom_left_corner, &sln_supplement));

  Views::BaseView<double> b;
  b.show(&space_whole_domain);
  Views::View::wait();
  b.show(&space_bottom_left_corner);
  Views::View::wait();
  b.show(&space_supplement);
  Views::View::wait();

  return 0;
}