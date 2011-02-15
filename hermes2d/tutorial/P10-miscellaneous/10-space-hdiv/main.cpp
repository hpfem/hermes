#include "hermes2d.h"

// This example shows how to use the Hdiv space and
// visualize finite element basis functions. Note that 
// higher-order basis functions in this space comprise 
// edge functions associated with mesh edges (normal 
// component is zero on the boundary of the element patch
// associated with the edge), and bubble functions 
// associated with elements (normal component is 
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

  // Enter boundary markers (default is Neumann boundary).
  BCTypes bc_types;

  // Enter Dirichlet boundary values (not relevant here).
  BCValues bc_values;

  // Create an Hdiv space with default shapeset.
  // (BC types and essential BC values not relevant.)
  HdivSpace space(&mesh, &bc_types, &bc_values, P_INIT);

  // Visualise the FE basis.
  VectorBaseView bview("VectorBaseView", new WinGeom(0, 0, 700, 600));
  bview.show(&space);

  // Wait for all views to be closed.
  View::wait();
  return 0;
}

