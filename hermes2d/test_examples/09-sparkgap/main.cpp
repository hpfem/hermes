#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

// Problem parameters.
const double SRC = 1000.;
const std::vector<std::string> SRC_BOUNDARY = { "0", "5" };
const std::vector<std::string> GND_BOUNDARY = { "7", "8" };

// Error calculation & adaptivity.
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.75;
// Error calculator
DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Selector with predefined list of element refinement candidates.
H1ProjBasedSelector<double> selector(H2D_HP_ANISO);
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-8;

int main(int argc, char* argv[])
{
	// Load the mesh.
	MeshSharedPtr mesh(new Mesh);
	Hermes::Hermes2D::MeshReaderH2DXML mloader;
	mloader.load("sparkgap.msh", mesh);

	// Initialize essential boundary conditions.
	Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_src(SRC_BOUNDARY, SRC);
	Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_gnd(GND_BOUNDARY, 0.);
	Hermes::Hermes2D::EssentialBCs<double> bcs({ &bc_essential_src, &bc_essential_gnd });

	// Initialize space of the 3-rd order.
	SpaceSharedPtr<double> space(new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, 3));
	OrderView Oview("FE Space", new Hermes::Hermes2D::Views::WinGeom(100, 50, 600, 600));

	// Initialize the weak formulation.
	WeakFormSharedPtr<double> wf(new WeakFormsH1::DefaultWeakFormLaplaceLinear<double>(HERMES_ANY, HERMES_AXISYM_Y));

	// Initialize the solution.
	MeshFunctionSharedPtr<double> sln(new Solution<double>);
	MeshFunctionSharedPtr<double> ref_sln(new Solution<double>);

	// Initialize the view.
	Hermes::Hermes2D::Views::ScalarView Sview("Solution", new Hermes::Hermes2D::Views::WinGeom(750, 50, 600, 600));

	// Initialize linear solver.
	Hermes::Hermes2D::LinearSolver<double> linear_solver(wf, space);

	adaptivity.set_space(space);
	while (true)
	{
		// Construct globally refined reference mesh and setup reference space.
		Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
		MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
		Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
		SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();
		linear_solver.set_space(ref_space);

		// Solve the problem.
		linear_solver.solve();

		// Translate the solution vector into the previously initialized Solution.
		Hermes::Hermes2D::Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), ref_space, ref_sln);

		// Visualize the solution.
		Oview.show(ref_space);
		Sview.show(ref_sln);

		// Project the reference solution to the coarse space in H1 norm for error calculation.
		OGProjection<double>::project_global(space, ref_sln, sln);

		// Calculate element errors and total error estimate.
		errorCalculator.calculate_errors(sln, ref_sln);
		double total_error_estimate = errorCalculator.get_total_error_squared() * 100;
		Sview.set_title("Solution, spatial error: %g%%.", total_error_estimate);

		// If error is too large, adapt the (coarse) space.
		if (total_error_estimate > ERR_STOP)
			adaptivity.adapt(&selector);
		else
			break;
	}

	// Wait for view to be closed.
	Views::View::wait();
	return 0;
}