#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"
#define PI 4.0*atan(1.0)
#include "hermes2d.h"

#include <iostream>
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
// Boundary condition types; this does not really affect the testing
// Set your bc_types() according to the CASE
//
BCType bc_types(int marker)
{
  if      (marker == 1) return BC_ESSENTIAL; 
  else if (marker == 2) return BC_NATURAL; 
  else if (marker == 3) return BC_ESSENTIAL; 
  else if (marker == 4) return BC_NATURAL; 
  else                  return BC_NATURAL;
} // end of bc_types()

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
    for(int edge = 0; edge < e->nvert; ++edge) {
      if ((e->en[edge]->bnd) && (e->en[edge]->marker == bdryMarker)) {
        rm.set_active_element(e);
        points_location = quad->get_edge_points(edge);
        points = quad->get_points(points_location);
        np = quad->get_num_points(points_location);
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

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

//******************************************************************************
// Main
//
int main(int argc, char* argv[])
{
  if (argc != 3)
  {
    printf("Missing mesh filename and domain perimeter as command-line parameters.");
    return ERROR_FAILURE;
  }

  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
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
  for (int i=0; i<INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Create H1 space with a default shapeset.
  H1Space space(&mesh, bc_types, NULL, P_INIT);

  // FIXME: This Solution is artificial here and it should be removed. The 
  // function CalculateBoundaryLength() should only take a pointer to Mesh and 
  // a boundary marker as parameters.
  //Solution sln;
  //sln.set_zero(&mesh);

  // Calculate the length of the four boundaries segments.
  double l1 = CalculateBoundaryLength(&mesh, 1);
  info("Length of boundary 1 = %g\n", l1);

  double l2 = CalculateBoundaryLength(&mesh, 2);
  info("Length of boundary 2 = %g\n", l2);

  double l3 = CalculateBoundaryLength(&mesh, 3);
  info("Length of boundary 3 = %g\n", l3);

  double l4 = CalculateBoundaryLength(&mesh, 4);
  info("Length of boundary 4 = %g\n", l4);
  
  double perimeter = l1 + l2 + l3 + l4;
  info("Computed perimeter = %10.15f\n", perimeter);

  // Set exact value from CMakeLists.txt file
  info("Exact perimeter = %g\n", bdryLengthInput);

  if (fabs(perimeter - bdryLengthInput) < 1e-6) {
    printf("Success!\n");
    return ERROR_SUCCESS;
  }
  else {
    printf("Failure!\n");
    return ERROR_FAILURE;
  }
 

  return 0;
}
