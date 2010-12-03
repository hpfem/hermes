#include "hermes2d.h"

// This test makes sure that example 31-space-hdiv works correctly.

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

  // Enter boundary markers.
  // (If no markers are entered, default is a natural BC).
  BCTypes bc_types;

  // Create an Hdiv space with default shapeset.
  // (BC types and essential BC values not relevant.)
  HdivSpace space(&mesh, &bc_types, NULL, P_INIT);

  // Visualise the FE basis.
  VectorBaseView bview("VectorBaseView", new WinGeom(0, 0, 700, 600));

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

