#define HERMES_REPORT_INFO
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

// This test makes sure that subdomains work correctly.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh_whole_domain, mesh_bottom_left_corner, mesh_supplement;
  MeshReaderH2DXML mloader;
  mloader.load("subdomains.xml", Hermes::vector<Mesh*>(&mesh_whole_domain, &mesh_bottom_left_corner, &mesh_supplement));

  Views::MeshView m;
  m.show(&mesh_whole_domain);
  m.wait_for_keypress();
  m.show(&mesh_bottom_left_corner);
  m.wait_for_keypress();
  m.show(&mesh_supplement);
  m.wait_for_close();

  if(mesh_whole_domain.get_num_elements() == 4 && mesh_bottom_left_corner.get_num_elements() == 1 && mesh_supplement.get_num_elements() == 3)
  {
    info("Success!");
    return TEST_SUCCESS;
  }
  else
  {
    info("Failure!");
    return TEST_FAILURE;
  }
}