#define HERMES_REPORT_ALL
#include "hermes2d.h"

// This example shows how to load a mesh, perform various types
// of initial refinements, and use keyboard and mouse controls.
//
// Geometry: L-Shape domain (see file domain.mesh).

static char text[] = "\
Click into the image window and:\n\
  press 'm' to show/hide element material markers,\n\
  press 'i' to show/hide element indices,\n\
  press 'b' to toggle boundary markers,\n\
  enlarge your window and press 'c' to center the mesh,\n\
  zoom into the mesh using the right mouse button\n\
  move the mesh around using the left mouse button\n\
  press 'c' to center the mesh again,\n\
  press 'h' to render high resolution image,\n\
  press 's' to save a screenshot in bmp format\n\
  press 'q' to quit.\n\
  Press F1 for help.\n";

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);

  // Optional rescaling of mesh (all vertex x- and y-coordinates are 
  // divided by x_ref and y_ref, respectively). Mesh with curved edges 
  // cannot be rescaled. So to try this feature, comment out the "curves" 
  // section in the mesh file.
  double x_ref = 2.0, y_ref = 3.0;
  if(!mesh.rescale(x_ref, y_ref)) info("Mesh was not rescaled.");
  else {
    info("Mesh scaled by the factors of %g and %g in the x- and y- direction, respectively.", 
         x_ref, y_ref);
  }

  // Conversion between triangular and quadrilateral meshes (optional). 
  //mesh.convert_quads_to_triangles();
  //mesh.convert_triangles_to_quads();

  // Refine mesh uniformly (optional).
  mesh.refine_all_elements();          

  // Refine towards a mesh vertex (optional).
  mesh.refine_towards_vertex(3, 4);    // Four refinements towards vertex no. 3.

  // Refine towards boundary (optional).
  mesh.refine_towards_boundary("Outer", 4);  // Four successive refinements towards 
                                             // boundary with marker "Outer".

  // Refine individual elements (optional).
  mesh.refine_element_id(86, 0);          // 0... isotropic refinement.
  mesh.refine_element_id(112, 0);         // 0... isotropic refinement.
  mesh.refine_element_id(84, 2);          // 2... anisotropic refinement.
  mesh.refine_element_id(114, 1);         // 1... anisotropic refinement.

  // Display the mesh.
  // (0, 0) is the upper left corner position, 
  // 350 x 350 is the window size.
  MeshView mview("Hello world!", new WinGeom(0, 0, 350, 350));
  mview.show(&mesh);

  // Practice some keyboard and mouse controls.
  printf("%s", text);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}
