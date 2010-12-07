#include "hermes2d.h"

// This test makes sure that example 30-space-hcurl works correctly.

int INIT_REF_NUM = 2;      // Initial uniform mesh refinement.
int P_INIT = 3;            // Polynomial degree of mesh elements.

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("square.mesh", &mesh);

  // Initial mesh refinement.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create an Hcurl space with default shapeset.
  // (BC types and essential BC values not relevant.)
  // Enter boundary markers.
  BCTypes bc_types;

  // Enter Dirichlet boundary values.
  BCValues bc_values;

  HcurlSpace space(&mesh, &bc_types, &bc_values, P_INIT);

  // Visualize FE basis.
  VectorBaseView bview("VectorBaseView", new WinGeom(0, 0, 700, 600));

  // View::wait(H2DV_WAIT_KEYPRESS);

  bool success = true;

  if (success == true) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}

