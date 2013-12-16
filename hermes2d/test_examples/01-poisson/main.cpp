#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

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

const bool HERMES_VISUALIZATION = true;   // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = true;     // Set to "true" to enable VTK output.
const int P_INIT = 1;                     // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;               // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5e2;        // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 20.0;        // Fixed temperature on the boundary.

class MyVolumetricIntegralCalculator : public PostProcessing::SurfaceIntegralCalculator<double>
{
public:
  MyVolumetricIntegralCalculator(MeshFunctionSharedPtr<double> source_function, int number_of_integrals) : PostProcessing::SurfaceIntegralCalculator<double>(source_function, number_of_integrals)
  {
  }

  MyVolumetricIntegralCalculator(Hermes::vector<MeshFunctionSharedPtr<double> > source_functions, int number_of_integrals) : PostProcessing::SurfaceIntegralCalculator<double>(source_functions, number_of_integrals)
  {
  }

  virtual void integral(int n, double* wt, Func<double> **fns, Geom<double> *e, double* result)
  {
    for (int i = 0; i < n; i++)
      result[0] += wt[i] * fns[0]->val[i];
  };

  virtual void order(Func<Hermes::Ord> **fns, Hermes::Ord* result) {
    result[0] = Hermes::Ord(21);
  }
};

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);

  HermesCommonApi.set_integral_param_value(numThreads, 1);

  // Refine all elements, do it INIT_REF_NUM-times.
  for(unsigned int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"),
    FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // Initialize space->
  SpaceSharedPtr<double> space( new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, P_INIT));

  std::cout << "Ndofs: " << space->get_num_dofs() << std::endl;

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Aluminum", new Hermes::Hermes1DFunction<double>(LAMBDA_AL), "Copper",
    new Hermes::Hermes1DFunction<double>(LAMBDA_CU), new Hermes::Hermes2DFunction<double>(-VOLUME_HEAT_SRC));

  // Initialize the solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);
  
  // Solve the linear problem.
  try
  {
    linear_solver.solve();

    // Get the solution vector.
    double* sln_vector = linear_solver.get_sln_vector();

    // Translate the solution vector into the previously initialized Solution.
    Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, sln);

    Hermes::Hermes2D::Views::ScalarView viewSO("Solution", new Hermes::Hermes2D::Views::WinGeom(0, 0, 500, 400));
 //   viewSO.show(sln, 1.0);
  //  viewSO.wait();

    Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(0, 0, 500, 400));
    Hermes::Hermes2D::Views::MeshView viewM("Solution", new Hermes::Hermes2D::Views::WinGeom(0, 0, 500, 400));
    Views::Linearizer lin(Views::FileExport);
    Views::Linearizer lin2;
    MeshFunctionSharedPtr<double> slns[2] = { sln, sln };
    int items[2] = { 1, 1 };
    lin.save_solution_vtk(slns, items, "asf.vtk", "as", false, 1.);
    lin2.save_solution_vtk(sln, "asf2.vtk", "asdf", 1, true, 1.);
    viewS.show(sln, 1.0);
    viewS.wait();
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }
  return 0;
}