#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#define PI 4.0*atan(1.0)
#include "hermes2d.h"

#include <iostream>

using namespace Hermes;
using namespace Hermes::Hermes2D;
using std::cout;
using std::endl;

//******************************************************************************
//  This is a very basic test for computing the perimeter of a domain.
//  The user can experiment with existing cases or contribute a new case.
//  The latter may require an additional mesh file to be created and named
//  domain.meshX where X should be an integer following the largest already
//  used in the provided set of domain.meshX files in this directory. The
//  header in the existing mesh files should be helpful.

//  To add more domains to test, edit the CMakeLists.txt file. Note: The
//  boundary markers must be 1, 2, 3, or 4, others will be skipped.

//******************************************************************************
// Controls
//
//  The following parameters can be changed:

const int P_INIT = 1;       // Initial polynomial degree of all mesh elements.
const int INIT_REF_NUM = 2; // Number of initial uniform mesh refinements.

//******************************************************************************
// Helper functions

//------------------------------------------------------------------------------
// Compute marked boundary length
//
double CalculateBoundaryLength(Mesh* mesh, int bdryMarker)
{
  // Variables declaration.
  Element* e;
  double length = 0;
  RefMap rm;
  rm.set_quad_2d(&g_quad_2d_std);
  Quad2D * quad = rm.get_quad_2d();
  int points_location;
  double3* points;
  int np;
  double3* tangents;

  // Loop through all boundary faces of all active elements.
  for_all_active_elements(e, mesh) {
    for(int edge = 0; edge < e->get_num_surf(); ++edge) {
      if((e->en[edge]->bnd) && (e->en[edge]->marker == bdryMarker)) {
        rm.set_active_element(e);
        points_location = quad->get_edge_points(edge, quad->get_max_order(e->get_mode()), e->get_mode());
        points = quad->get_points(points_location, e->get_mode());
        np = quad->get_num_points(points_location, e->get_mode());
        tangents = rm.get_tangent(edge, points_location);
        for(int i = 0; i < np; i++) {
          // Weights sum up to two on every edge, therefore the division by two must be present.
          length +=  0.5 * points[i][2] * tangents[i][2];
        }
      }
    }
  }
  return length;
} // end of CalculateBoundaryLength()

//******************************************************************************
// Main
//
int main(int argc, char* argv[])
{
  if(argc != 3)
  {
    printf("Missing mesh filename and domain perimeter as command-line parameters.");
    return -1;
  }

  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load(argv[1], &mesh);

  /*
  //Graphics not allowed on test versions; user may uncomment for debugging; vfda
  // Display the mesh.
  // (100, 0) is the upper left corner position, 600 x 500 is the window size
  MeshView mview("Mesh", 100, 0, 600, 500);
  mview.show(&mesh);
  // Wait for the view to be closed.
  View::wait();
  */

  double bdryLengthInput = atof(argv[2]);

  // Perform initial mesh refinements.
  for (int i = 0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create H1 space with a default shapeset.
  H1Space<double> space(&mesh, P_INIT);

  // FIXME: This Solution is artificial here and it should be removed. The
  // function CalculateBoundaryLength() should only take a pointer to Mesh and
  // a boundary marker as parameters.
  //Solution sln;
  //sln.set_zero(&mesh);

  // Calculate the length of the four boundaries segments.
  double l1 = CalculateBoundaryLength(&mesh, 1);

  double l2 = CalculateBoundaryLength(&mesh, 2);

  double l3 = CalculateBoundaryLength(&mesh, 3);

  double l4 = CalculateBoundaryLength(&mesh, 4);

  double perimeter = l1 + l2 + l3 + l4;

  // Set exact value from CMakeLists.txt file

  if(fabs(perimeter - bdryLengthInput) < 1e-6) {
    printf("Success!\n");
    return 0;
  }
  else {
    printf("Failure!\n");
    return -1;
  }

  return 0;
}