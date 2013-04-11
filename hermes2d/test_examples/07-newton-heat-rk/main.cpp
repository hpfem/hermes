#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  This example is a continuation of the example "09-timedep-basic" and it shows how
//  to perform time integration with arbitrary Runge-Kutta methods, using Butcher's
//  tables as input parameters. Currently (as of January 2011) approx. 30 tables are
//  available by default, as can be seen below. They are taken from various sources
//  including J. Butcher's book, journal articles, and the Wikipedia page
//  http://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods. If you know
//  about some other interesting R-K method that we are missing here, please let us
//  know!
//
//  PDE: non-stationary heat transfer equation
//       HEATCAP * RHO * dT/dt - LAMBDA * Laplace T = 0.
//  This equation is, however, written in such a way that the time-derivative
//  is on the left and everything else on the right:
//
//  dT/dt = LAMBDA * Laplace T / (HEATCAP * RHO).
//
//  As opposed to the previous example where the time-discretization was hardwired
//  into the weak formulation, now we only need the weak formulation of the
//  right-hand side.
//
//  Domain: St. Vitus cathedral (file cathedral.mesh).
//
//  IC:  T = TEMP_INIT.
//  BC:  T = TEMP_INIT on the bottom edge ... Dirichlet,
//       LAMBDA * dT/dn = ALPHA*(t_exterior(time) - T) ... Newton, time-dependent.
//
//  Time-stepping: Arbitrary Runge-Kutta methods.
//

//  The following parameters can be changed:
const bool REUSE_SOLUTION = false;
const int P_INIT = 2;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 1;                   // Number of initial uniform mesh refinements towards the boundary.
double time_step = 3e+2;                    // Time step in seconds.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int max_allowed_iterations = 100;                  // Maximum allowed number of Newton iterations.
MatrixSolverType matrix_solver_type = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Choose one of the following time-integration methods, or define your own Butcher's table. The last number
// in the name of each method is its order. The one before last, if present, is the number of stages.
// Explicit methods:
//   Explicit_RK_1, Explicit_RK_2, Explicit_RK_3, Explicit_RK_4.
// Implicit methods:
//   Implicit_RK_1, Implicit_Crank_Nicolson_2_2, Implicit_SIRK_2_2, Implicit_ESIRK_2_2, Implicit_SDIRK_2_2,
//   Implicit_Lobatto_IIIA_2_2, Implicit_Lobatto_IIIB_2_2, Implicit_Lobatto_IIIC_2_2, Implicit_Lobatto_IIIA_3_4,
//   Implicit_Lobatto_IIIB_3_4, Implicit_Lobatto_IIIC_3_4, Implicit_Radau_IIA_3_5, Implicit_SDIRK_5_4.
// Embedded explicit methods:
//   Explicit_HEUN_EULER_2_12_embedded, Explicit_BOGACKI_SHAMPINE_4_23_embedded, Explicit_FEHLBERG_6_45_embedded,
//   Explicit_CASH_KARP_6_45_embedded, Explicit_DORMAND_PRINCE_7_45_embedded.
// Embedded implicit methods:
//   Implicit_SDIRK_CASH_3_23_embedded, Implicit_ESDIRK_TRBDF2_3_23_embedded, Implicit_ESDIRK_TRX2_3_23_embedded,
//   Implicit_SDIRK_BILLINGTON_3_23_embedded, Implicit_SDIRK_CASH_5_24_embedded, Implicit_SDIRK_CASH_5_34_embedded,
//   Implicit_DIRK_ISMAIL_7_45_embedded.

ButcherTableType butcher_table_type = Implicit_RK_1;

// Problem parameters.
const double TEMP_INIT = 10;       // Temperature of the ground (also initial temperature).
const double ALPHA = 10;           // Heat flux coefficient for Newton's boundary condition.
const double LAMBDA = 1e2;         // Thermal conductivity of the material.
const double HEATCAP = 1e2;        // Heat capacity.
const double RHO = 3000;           // Material density.
const double T_FINAL = 86400;      // Length of time interval (24 hours) in seconds.

int main(int argc, char* argv[])
{
	// Choose a Butcher's table or define your own.
	ButcherTable bt(butcher_table_type);

	// Class for progressive saving of solutions and all necessary classes for the computation.
	CalculationContinuity<double> continuity(CalculationContinuity<double>::onlyTime);

	// Initialize the time.
	double current_time = 0;

	// mesh->
	MeshSharedPtr mesh(new Mesh);

	// Init mesh->
	MeshReaderH2D mloader;
	mloader.load("cathedral.mesh", mesh);

	// Perform initial mesh refinements.
	for(int i = 0; i < INIT_REF_NUM; i++)
		mesh->refine_all_elements();
	mesh->refine_towards_boundary("Boundary_air", INIT_REF_NUM_BDY);
	mesh->refine_towards_boundary("Boundary_ground", INIT_REF_NUM_BDY);

	// Initialize boundary conditions.
	Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential("Boundary_ground", TEMP_INIT);
	Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

	// space->
	SpaceSharedPtr<double> space;

	// Solution pointer.
	MeshFunctionSharedPtr<double> sln_time_prev = new ConstantSolution<double>(mesh, TEMP_INIT);

	if(REUSE_SOLUTION && continuity.have_record_available())
	{
		try
		{
			continuity.get_last_record()->load_mesh(mesh);
			space = continuity.get_last_record()->load_space(mesh, &bcs);
			continuity.get_last_record()->load_solution(sln_time_prev, space);
			Views::ScalarView s;
			s.show(sln_time_prev, Views::HERMES_EPS_NORMAL, H2D_FN_VAL_0);
			s.wait_for_close();
			current_time = continuity.get_last_record()->get_time();
		}
		catch(Hermes2D::CalculationContinuityException& e)
		{
			e.print_msg();
		}
	}
	else
	{
		space = new H1Space<double>(mesh, &bcs, P_INIT);
	}


	CustomWeakFormHeatRK wf("Boundary_air", ALPHA, LAMBDA, HEATCAP, RHO,
		&current_time, TEMP_INIT, T_FINAL);

	MeshFunctionSharedPtr<double> sln_time_new = new Solution<double>(mesh);

	int ndof = space->get_num_dofs();

	// Initialize views.
	Hermes::Hermes2D::Views::ScalarView Tview("Temperature", new Hermes::Hermes2D::Views::WinGeom(0, 0, 450, 600));
	Tview.set_min_max_range(0, 20);
	Tview.fix_scale_width(30);

	// Initialize Runge-Kutta time stepping.
	RungeKutta<double> runge_kutta(&wf, space, &bt);

	runge_kutta.set_verbose_output(true);
	runge_kutta.output_matrix(1);
	runge_kutta.set_matrix_number_format("%a");
	runge_kutta.output_rhs(2);
	runge_kutta.set_global_integration_order(10);

	// Iteration number.
	int iteration = 0;

	// Time stepping loop:
	do
	{
		// Perform one Runge-Kutta time step according to the selected Butcher's table.
		try
		{
			runge_kutta.set_space(space);
			runge_kutta.set_time(current_time);
			runge_kutta.set_time_step(time_step);
			runge_kutta.rk_time_step_newton(sln_time_prev, sln_time_new);
		}
		catch(Exceptions::Exception& e)
		{
			e.print_msg();
		}

		// Show the new time level solution.
		char title[100];
		sprintf(title, "Time %3.2f s", current_time);
		Tview.set_title(title);
		Tview.show(sln_time_new);

		// Save the progress.
		continuity.add_record(current_time, mesh, space, sln_time_prev);

		// Copy solution for the new time step.
		sln_time_prev->copy(sln_time_new);

		// Increase current time and time step counter.
		current_time += time_step;
		iteration++;
	}
	while (current_time < T_FINAL);

	// Wait for the view to be closed.
	Hermes::Hermes2D::Views::View::wait();
	return 0;
}
