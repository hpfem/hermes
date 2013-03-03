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
const bool VTK_VISUALIZATION = false;             // Set to "true" to enable VTK output.
const int P_INIT = 2;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 2;                       // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5e2;        // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 20.0;        // Fixed temperature on the boundary.

int main(int argc, char* argv[])
{
  // Load the mesh.
  Hermes::Hermes2D::Mesh mesh;
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  try
  {
    mloader.load("domain.xml", &mesh);
  }
  catch(Exceptions::MeshLoadFailureException& e)
  {
    e.print_msg();
    return -1;
  }

  // Refine all elements, do it INIT_REF_NUM-times.
  for(unsigned int i = 0; i < INIT_REF_NUM; i++)
    mesh.refine_all_elements();
  
  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"),
    FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // Initialize space.
  Hermes::Hermes2D::H1Space<double> space(&mesh, &bcs, P_INIT);
  
  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Aluminum", new Hermes::Hermes1DFunction<double>(LAMBDA_AL), "Copper",
    new Hermes::Hermes1DFunction<double>(LAMBDA_CU), new Hermes::Hermes2DFunction<double>(-VOLUME_HEAT_SRC));
  
  // Illustration of setting element orders.
  Hermes::Hermes2D::Element* e;
  int i = 1;
  for_all_active_elements(e, &mesh)
  {
    space.set_element_order(e->id, i++ % 4 + 2);
  }

  // One has to call this method after any changes to DOFs.
  space.assign_dofs();

  // Output of numbers of DOFs, vertex DOFs, edge DOFs, and bubble DOFs respectively.
  std::cout << space.get_num_dofs() << std::endl;
  std::cout << space.get_vertex_functions_count() << std::endl;
  std::cout << space.get_edge_functions_count() << std::endl;
  std::cout << space.get_bubble_functions_count() << std::endl;

 if (HERMES_VISUALIZATION)
 {
   Hermes::Hermes2D::Views::BaseView<double> o;
   o.show(&space);
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

    if(HERMES_VISUALIZATION)
    {
      // Visualize the solution.
      Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));

      viewS.show(&sln, Hermes::Hermes2D::Views::HERMES_EPS_LOW);
      viewS.wait_for_close();
    }
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  return 0;
}