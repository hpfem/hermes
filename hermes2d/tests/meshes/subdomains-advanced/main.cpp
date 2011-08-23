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
  Mesh mesh_whole_domain, mesh_without_hole;
  Hermes::vector<Mesh*> meshes (&mesh_whole_domain, &mesh_without_hole);
  MeshReaderH2DXML mloader;
  mloader.load("subdomains.xml", meshes);

  // Perform initial mesh refinements (optional).
  for(int i = 0; i < INIT_REF_NUM; i++)
    for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
      meshes[meshes_i]->refine_all_elements();

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf;

  // Initialize essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential_whole_domain("b1", 0.0);
  EssentialBCs<double> bcs_whole_domain(&bc_essential_whole_domain);

  DefaultEssentialBCConst<double> bc_essential_without_hole(Hermes::vector<std::string>("b1", "b2"), 0.0);
  EssentialBCs<double> bcs_without_hole(&bc_essential_without_hole);
  
  // Create H1 spaces with default shapeset.
  H1Space<double> space_whole_domain(&mesh_whole_domain, &bcs_whole_domain, P_INIT);
  int ndof_whole_domain = space_whole_domain.get_num_dofs();
  info("Space whole domain ndof = %d", ndof_whole_domain);

  H1Space<double> space_without_hole(&mesh_without_hole, &bcs_without_hole, P_INIT);
  int ndof_without_hole = space_without_hole.get_num_dofs();
  info("Space bottom left corner ndof = %d", ndof_without_hole);

  Views::BaseView<double> b;
  b.show(&space_whole_domain);
  b.wait_for_keypress();
  b.show(&space_without_hole);
  b.wait_for_close();

  // Initialize the FE problem.
  Hermes::Hermes2D::DiscreteProblem<double> dp(&wf, Hermes::vector<Space<double>*>(&space_whole_domain, &space_without_hole));

  // Initial coefficient vector for the Newton's method.  
  double* coeff_vec = new double[ndof_whole_domain + ndof_without_hole];
  memset(coeff_vec, 0, (ndof_whole_domain + ndof_without_hole)*sizeof(double));

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::Solution<double> sln_whole_domain, sln_without_hole;
  Hermes::Hermes2D::NewtonSolver<double> newton(&dp, matrix_solver_type);
  if (!newton.solve(coeff_vec)) 
    error("Newton's iteration failed.");
  else
    Hermes::Hermes2D::Solution<double>::vector_to_solutions(newton.get_sln_vector(), Hermes::vector<Space<double>*>(&space_whole_domain, &space_without_hole), Hermes::vector<Solution<double>*>(&sln_whole_domain, &sln_without_hole));

  Views::ScalarView<double> s("Solution on the whole domain");
  s.show(&sln_whole_domain);
  s.wait_for_keypress();
  s.set_title("Solution on the complement of the hole");
  s.show(&sln_without_hole);
  s.wait_for_close();

  return 0;
}