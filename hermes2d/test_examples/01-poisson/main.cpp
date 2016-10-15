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

// Set to "false" to suppress Hermes OpenGL visualization.
const bool HERMES_VISUALIZATION = true;
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;
// Uniform polynomial degree of mesh elements.
const int P_INIT = 10;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 4;

// Problem parameters.
// Thermal cond. of Al, Cu for temperatures around 20 deg Celsius.
const double LAMBDA_AL = 236.0, LAMBDA_CU = 386.0;
// Volume heat sources generated (for example) by electric current.
const double VOLUME_HEAT_SRC = 5;
// Fixed temperature on the boundary.
const double FIXED_BDY_TEMP = 20;

int main(int argc, char* argv[])
{
	// Load the mesh.
	MeshSharedPtr mesh(new Mesh);
	Hermes::Hermes2D::MeshReaderH2DXML mloader;
	mloader.load("domain.xml", mesh);

	// Refine all elements, do it INIT_REF_NUM-times.
	for (unsigned int i = 0; i < INIT_REF_NUM; i++)
		mesh->refine_all_elements();

	// Initialize essential boundary conditions.
	Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential({ "Bottom", "Inner", "Outer", "Left" }, FIXED_BDY_TEMP);
	Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

	// Initialize space.
	SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, P_INIT));
	std::cout << "Ndofs: " << space->get_num_dofs() << std::endl;

	// Initialize the weak formulation.
	WeakFormSharedPtr<double> wf(new CustomWeakFormPoisson("Aluminum", new Hermes::Hermes1DFunction<double>(LAMBDA_AL), "Copper",
		new Hermes::Hermes1DFunction<double>(LAMBDA_CU), new Hermes::Hermes2DFunction<double>(VOLUME_HEAT_SRC)));

	// Initialize the solution.
	MeshFunctionSharedPtr<double> sln(new Solution<double>);

	// Initialize linear solver.
	Hermes::Hermes2D::LinearSolver<double> linear_solver(wf, space);

	// Solve the problem.
	linear_solver.solve();

	// Translate the solution vector into the previously initialized Solution.
	Hermes::Hermes2D::Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), space, sln);

	// Visualize the solution.
	Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(750, 50, 600, 600));
	viewS.get_linearizer()->set_criterion(Views::LinearizerCriterionFixed(3));
	viewS.show(sln);

	// Wait for view to be closed.
	Views::View::wait();
	return 0;
}