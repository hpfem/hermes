#include "hermes2d.h"

// This test makes sure that example 01-mesh works correctly.

int main(int argc, char* argv[])
{
  // load the mesh file
  Mesh mesh;
  H2DReader mloader;
  mloader.load("../domain.mesh", &mesh);

  // perform some sample initial refinements
  mesh.refine_all_elements();           // Refines all elements.
  mesh.refine_towards_vertex(3, 4);     // Refines mesh towards vertex #3 (4x).
  mesh.refine_towards_boundary("Outer", 4); // Refines all elements along boundary 2 (4x).
  mesh.refine_element_id(86, 0);        // Refines element #86 isotropically.
  mesh.refine_element_id(112, 0);       // Refines element #112 isotropically.
  mesh.refine_element_id(84, 2);        // Refines element #84 anisotropically.
  mesh.refine_element_id(114, 1);       // Refines element #114 anisotropically.

  // new code for the test
  int n_elem = mesh.get_num_elements();
  printf("n_elem = %d\n", n_elem);
  int n_active = mesh.get_num_active_elements();
  printf("n_active = %d\n", n_active);
  int n_base = mesh.get_num_base_elements();
  printf("n_base = %d\n", n_base);

  if (n_elem == 576 && n_active == 424 && n_base == 4) {
    printf("Success!\n");
    return ERR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERR_FAILURE;
  }
}
