#define HERMES_REPORT_ALL
#include "definitions.h"

// This example shows how to solve a simple PDE that describes stationary
// heat transfer in an object consisting of two materials (aluminum and
// copper). The object is heated by constant volumetric heat sources
// (generated, for example, by a DC electric current). The temperature
// on the boundary is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - initialize matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: L-Shape domain (see file domain.mesh).
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = true;              // Set to "true" to enable VTK output.
const int P_INIT = 1;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5e2;        // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 20.0;        // Fixed temperature on the boundary.

int main(int argc, char* argv[])
{
  // Set the number of threads used in Hermes.
  Hermes::HermesCommonApi.set_param_value(Hermes::exceptionsPrintCallstack, 0);
  Hermes::Hermes2D::Hermes2DApi.set_param_value(Hermes::Hermes2D::numThreads, 8);

  // Load the mesh.
  Hermes::Hermes2D::Mesh mesh;
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  try
  {
    mloader.load("domain.xml", &mesh);
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  // Perform initial mesh refinements (optional).
  mesh.refine_in_areas(Hermes::vector<std::string>("Aluminum", "Copper"), INIT_REF_NUM);
  mesh.refine_in_area("Aluminum");

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Aluminum", new Hermes::Hermes1DFunction<double>(LAMBDA_AL), "Copper",
    new Hermes::Hermes1DFunction<double>(LAMBDA_CU), new Hermes::Hermes2DFunction<double>(-VOLUME_HEAT_SRC));

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"),
    FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  Hermes::Hermes2D::H1Space<double> space(&mesh, &bcs, P_INIT);

  Hermes::Hermes2D::Element* e;
  int i = 1;
  for_all_active_elements(e, &mesh)
  {
    space.set_element_order(e->id, i++ % 9 + 1);
  }

  // Initialize the solution.
  Hermes::Hermes2D::Solution<double> sln;

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, &space);

  // Solve the linear problem.
  try
  {
    linear_solver.solve();

    // Get the solution vector.
    double* sln_vector = linear_solver.get_sln_vector();

    // Translate the solution vector into the previously initialized Solution.
    Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, &space, &sln);

    // VTK output.
    if(VTK_VISUALIZATION)
    {
      // Output solution in VTK format.
      Hermes::Hermes2D::Views::Linearizer lin;
      bool mode_3D = false;
      lin.save_solution_vtk(&sln, "sln.vtk", "Temperature", mode_3D, 1, Hermes::Hermes2D::Views::HERMES_EPS_LOW);

      // Output mesh and element orders in VTK format.
      Hermes::Hermes2D::Views::Orderizer ord;
      ord.save_mesh_vtk(&space, "mesh.vtk");
      ord.save_orders_vtk(&space, "ord.vtk");
    }

    // Visualize the solution.
    Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));
    Hermes::Hermes2D::Views::OrderView viewO("Orders", new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));
    Hermes::Hermes2D::Views::VectorView viewV("Vectors", new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));
    Hermes::Hermes2D::Views::MeshView viewM("Mesh", new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));

    if(HERMES_VISUALIZATION)
    {
      viewS.show(&sln);
      viewS.save_screenshot("000-base.bmp");
      viewS.show_contours(1.0);
      viewS.save_screenshot("001-contours1.bmp");
      viewS.show_contours(3.0);
      viewS.save_screenshot("002-contours2.bmp");
      viewS.show_contours(5.0);
      viewS.save_screenshot("003-contours3.bmp");
      viewS.set_3d_mode(true);
      viewS.save_screenshot("004-3D.bmp");
      viewS.set_3d_mode(false);
      viewS.save_screenshot("005-backTo2D.bmp");
      viewS.hide_contours();
      viewS.save_screenshot("006-noContours.bmp");
      viewS.set_palette_filter(false);
      viewS.save_screenshot("007-paletteSmooth.bmp");
      viewS.set_palette_filter(true);
      viewS.save_screenshot("008-paletteBack.bmp");
      viewS.show_mesh(false);
      viewS.save_screenshot("009-noMesh.bmp");
      viewS.set_scale_size(40, 1230, 14);
      viewS.save_screenshot("010-scaleSize.bmp");
      
      viewO.show(&space);
      viewO.save_screenshot("011-space.bmp");
      viewO.set_b_orders(true);
      viewO.save_screenshot("012-spaceBOrders.bmp");
      viewM.show(&mesh);
      viewM.save_screenshot("013-mesh.bmp");
      viewM.set_b_elem_mrk(true);
      viewM.save_screenshot("014-meshBElemMrk.bmp");
      viewV.show(&sln, &sln, Hermes::Hermes2D::Views::HERMES_EPS_VERYHIGH);
      viewV.save_screenshot("015-vectorizer.bmp");
      viewV.set_mode(2);
      viewV.save_screenshot("016-vectorizerArrowsMode.bmp");
    }
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  return 0;
}