#include "hermes2d.h"

// This example shows how to use the Hcurl space and
// visualize finite element basis functions. Note that 
// higher-order basis functions in this space comprise 
// edge functions associated with mesh edges (tangential 
// component is zero on the boundary of the element patch
// associated with the edge), and bubble functions 
// associated with elements (tangential component is 
// zero on the element boundary).
//
// The following parameters can be changed:

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
  HcurlSpace space(&mesh, (BCTypes *) NULL, NULL, P_INIT);

  // Visualize FE basis.
  VectorBaseView bview("VectorBaseView", new WinGeom(0, 0, 700, 600));
  bview.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

