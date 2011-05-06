#include "hermes2d.h"

// This example demonstrates how to create a space over a mesh
// and how to visualize the finite element basis functions
// using the BaseView class. The variable P_INIT defines the
// initial degree of all mesh elements. Currently, it is set
// to 1. After visualizing the basis functions using the hints
// printed in the terminal window, change the value of P_INIT
// to 2, 3, etc. to also see higher-order shape functions.
//
// You can use this example to visualize all shape functions
// on the reference square and reference triangle domains,
// just load the corresponding mesh at the beginning of the
// main.cpp file.
//
// Geometry: L-Shape domain (see file domain.mesh). It can be replaced
//           with "ref_square.mesh" and "ref_triangle.mesh" to visualize 
//           reference element shape functions. 

const int P_INIT = 3;

static char text[] = "\
Click into the image window and:\n\
  press 'f' to make the color scale finer/coarser,\n\
  press '3' to render 3D plot of the basis functions,\n\
  change the scale on the vertical axis using the * and / keys,\n\
  use all three mouse buttons to rotate/move/enlarge the graphs,\n\
  use the right/left arrows to browse through basis functions,\n\
  press 'l' to see adaptive linearizer output for OpenGL\n\
  use the left mouse button to drag the scale to another corner,\n\
  press 'p' to switch to greyscales and back,\n\
  press 'h' to render high-resolution image,\n\
  press 'm' to hide the mesh,\n\
  press 's' to save screenshots,\n\
  press 'q' to quit.\n\
  Press F1 for help.\n";

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  mloader.load("domain.mesh", &mesh);            // L-shape domain,

  // The following can be used to view higher-order shape functions
  // on reference domains (disable uniform mesh refinememts for that).
  //mloader.load("ref_square.mesh", &mesh);      // Reference square,
  //mloader.load("ref_triangle.mesh", &mesh);    // Reference triangle,

  // Refine all elements (optional).
  mesh.refine_all_elements();

  // Create an H1 space with default shapeset and natural BC.
  H1Space space(&mesh, P_INIT);

  // View FE basis functions.
  BaseView bview("Finite Element Space", new WinGeom(0, 0, 440, 350));
  bview.fix_scale_width(50);
  bview.show(&space);

  // Practice keyboard and mouse controls.
  printf("%s", text);

  // Wait for the view to be closed.
  View::wait();
  return 0;
}

