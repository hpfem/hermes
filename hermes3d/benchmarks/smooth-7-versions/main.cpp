#define H3D_REPORT_WARN
#define H3D_REPORT_INFO
#define H3D_REPORT_VERBOSE
#include "config.h"
#include <hermes3d.h>

//  This benchmark solves the Poisson equation with an exact solution from fichera corner. 
//  The problem geometry is a cube with missing corner. The exact solution has a singular 
//  gradient (an analogy of infinite stress) at the center. The knowledge of the exact 
//  solution makes it possible to calculate the approximation error exactly.  You can 
//  compare h- and hp-adaptivity from the point of view of both CPU time requirements and 
//  discrete problem size. Also look at the quality of the a-posteriori error estimator 
//  used by Hermes. 
//
//  PDE: -Laplace u = f.
//
//  Known exact solution, see functions fn() and fndd().
//
//  Domain: unit cube (0, 0, 1)x(0, 1, 0)x(1, 0, 0), with a missing corner, see the file "fichera-corner.mesh3d".
//
//  BC:  Essential Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

const int INIT_REF_NUM = 0;         // Number of initial uniform mesh refinements.
const int P_INIT_X = 2,
          P_INIT_Y = 2,
          P_INIT_Z = 2;             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;       // Error threshold for element refinement of the adapt(...) function 
                                    // (default) STRATEGY = 0 ... refine elements elements until sqrt(THRESHOLD) 
                                    // times total error is processed. If more elements have similar errors, 
                                    // refine all to keep the mesh symmetric.
                                    // STRATEGY = 1 ... refine all elements whose error is larger
                                    // than THRESHOLD times maximum element error.
const double ERR_STOP = 1e-4;       // Stopping criterion for adaptivity (rel. error tolerance between the
                                    // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 100000;       // Adaptivity process stops when the number of degrees of freedom grows
                                    // over this limit. This is to prevent h-adaptivity to go on forever.
bool solution_output = true;        // Generate output files (if true).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).

#undef REFERENCE_SOLUTION			// Use ref. solution to guide adaptivity. Otherwise
                                                // Exact solution is used.

#define ANISO_X						1
#define ANISO_Y						2
#define ANISO_Z						4

int init_q = 1;					// Initial order.
int init_p = 2;					// Initial order.
int ANISO_TYPE = 0;                             // Type of anisotropy.

// Exact solution
#include "exact_solution.cpp"

BCType bc_types(int marker)
{
  switch (ANISO_TYPE) {
    case ANISO_X: if (marker == 1 || marker == 2) return BC_ESSENTIAL; else return BC_NATURAL;
    case ANISO_Y: if (marker == 3 || marker == 4) return BC_ESSENTIAL; else return BC_NATURAL;
    case ANISO_Z: if (marker == 5 || marker == 6) return BC_ESSENTIAL; else return BC_NATURAL;
    case ANISO_X | ANISO_Y: if (marker == 1 || marker == 2 || marker == 3 || marker == 4) return BC_ESSENTIAL; else return BC_NATURAL;
    case ANISO_X | ANISO_Z: if (marker == 1 || marker == 2 || marker == 5 || marker == 6) return BC_ESSENTIAL; else return BC_NATURAL;
    case ANISO_Y | ANISO_Z: if (marker == 3 || marker == 4 || marker == 5 || marker == 6) return BC_ESSENTIAL; else return BC_NATURAL;
    case ANISO_X | ANISO_Y | ANISO_Z: return BC_ESSENTIAL;
  }

  return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}

void parse_aniso_type(char *str)
{
  int type = 0;
  if (strchr(str, 'x') != NULL) type |= ANISO_X;
  if (strchr(str, 'y') != NULL) type |= ANISO_Y;
  if (strchr(str, 'z') != NULL) type |= ANISO_Z;
  ANISO_TYPE = type;
}

bool check_order(const Ord3 &spord)
{
	switch (ANISO_TYPE)
	{
		case ANISO_X:
			if (spord.y == init_q && spord.z == init_q) return true;
			break;

		case ANISO_Y:
			if (spord.x == init_q && spord.z == init_q) return true;
			break;

		case ANISO_Z:
			if (spord.x == init_q && spord.y == init_q) return true;
			break;

		case ANISO_X | ANISO_Y:
			if (spord.z == init_q) return true;
			break;

		case ANISO_X | ANISO_Z:
			if (spord.y == init_q) return true;
			break;

		case ANISO_Y | ANISO_Z:
			if (spord.x == init_q) return true;
			break;

		case ANISO_X | ANISO_Y | ANISO_Z:
			return true;
	}

	return false;
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char **args)
{
	if (argc < 3) error("Not enough parameters.");

	parse_aniso_type(args[2]);

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh;
	H3DReader mesh_loader;
	mesh_loader.load(args[1], &mesh);

	Ord3 order;
	switch (ANISO_TYPE)
	{
		case ANISO_X: order = Ord3(init_p, init_q, init_q); break;
		case ANISO_Y: order = Ord3(init_q, init_p, init_q); break;
		case ANISO_Z: order = Ord3(init_q, init_q, init_p); break;
		case ANISO_X | ANISO_Y: order = Ord3(init_p, init_p, init_q); break;
		case ANISO_X | ANISO_Z: order = Ord3(init_p, init_q, init_p); break;
		case ANISO_Y | ANISO_Z: order = Ord3(init_q, init_p, init_p); break;
		case ANISO_X | ANISO_Y | ANISO_Z: order = Ord3(init_p, init_p, init_p); break;
	}

	printf("  - Setting uniform order to (%d, %d, %d)\n", order.x, order.y, order.z);
	H1Space space(&mesh, bc_types, essential_bc_values, order);

	WeakForm wf;
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM, HERMES_ANY);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY);

        bool is_linear = true;
	DiscreteProblem dp(&wf, &space, is_linear);

        // Initialize the solver in the case of SOLVER_PETSC or SOLVER_MUMPS.
        initialize_solution_environment(matrix_solver, argc, args);

	bool done = false;
	int iter = 0;
        Space* ref_space;
	do {
		printf("\n=== Iter #%d ================================================================\n", iter);

		// check if we are doing all right
		FOR_ALL_ACTIVE_ELEMENTS(eid, &mesh) {
			Ord3 spord = space.get_element_order(eid);
			printf("#%u: order = (%d, %d, %d)\n", eid, spord.x, spord.y, spord.z);
		}

		if (mesh.get_num_elements() != 1) {
			printf("failed\n");
			break;
		}

		Ord3 spord = space.get_element_order(1);
		if (!check_order(spord)) {
			printf("failed\n");
			break;
		}

		// we're good -> go ahead
		printf("\nSolution\n");

                // Set up the solver, matrix, and rhs according to the solver selection.
                SparseMatrix* matrix = create_matrix(matrix_solver);
                Vector* rhs = create_vector(matrix_solver);
                Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

		int ndofs = space.assign_dofs();
		printf("  - Number of DOFs: %d\n", ndofs);

		// assemble stiffness matrix
		printf("  - Assembling... "); fflush(stdout);
		dp.assemble(matrix, rhs);
		printf("done\n");

		// solve the stiffness matrix
		printf("  - Solving... "); fflush(stdout);
		bool solved = solver->solve();
		if (solved)
			printf("done\n");
		else {
			printf("failed\n");
			break;
		}

		Solution sln(&mesh);
		Solution::vector_to_solution(solver->get_solution(), &space, &sln);

		printf("Reference solution\n");

		Mesh ref_mesh;
		ref_mesh.copy(mesh);
		ref_mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

#ifdef REFERENCE_SOLUTION

		ref_space = space.dup(&ref_mesh);
		ref_space->copy_orders(space, 1);

		DiscreteProblem rdp(&wf, ref_space, is_linear);
		rdp.set_space(ref_space);

		int ref_ndof = ref_space->assign_dofs();
		printf("  - Number of DOFs: %d\n", ref_ndof);

		printf("  - Assembling... "); fflush(stdout);
		rdp.assemble(&mat, &rhs);
		printf("done\n");

		printf("  - Solving... "); fflush(stdout);
		bool rsolved = rsolver.solve();
		if (rsolved)
			printf("done\n");
		else {
			printf("failed\n");
			break;
		}
		Solution ref_sln(&ref_mesh);
		Solution::vector_to_solution(ref_solver.get_solution(), ref_space, ref_sln);
#else
		ExactSolution ref_sln(&ref_mesh, exact_solution);
#endif

		{
			double h1_err_norm = h1_error(&sln, &ref_sln) * 100;
			printf("  - H1 error norm:      % le\n", h1_err_norm);
		}

		printf("Adaptivity:\n");
                Adapt *adaptivity = new Adapt(&space, HERMES_H1_NORM);
		adaptivity->set_aniso(true);
                bool solutions_for_adapt = true;
                double err_est_rel = adaptivity->calc_err_est(&sln, &ref_sln, solutions_for_adapt) * 100;
		printf("  - tolerance: "); fflush(stdout);
		printf("% lf\n", err_est_rel);
		if (err_est_rel < ERR_STOP || iter == 8) {
			ExactSolution ex_sln(&mesh, exact_solution);
			printf("\nDone\n");
			// norm
			double h1_sln_norm = h1_norm(&sln);
			double h1_err_norm = h1_error(&sln, &ex_sln);
			printf("  - H1 solution norm:   % le\n", h1_sln_norm);
			printf("  - H1 error norm:      % le\n", h1_err_norm);

			double l2_sln_norm = l2_norm(&sln);
			double l2_err_norm = l2_error(&sln, &ex_sln);
			printf("  - L2 solution norm:   % le\n", l2_sln_norm);
			printf("  - L2 error norm:      % le\n", l2_err_norm);

			break;
		}

		printf("  - adapting... "); fflush(stdout);
		adaptivity->adapt(THRESHOLD);
		printf("done (refined %d element(s))\n", adaptivity->get_num_refined_elements());

                // Clean up.
#ifdef REFERENCE_SOLUTION
                delete ref_space->get_mesh();
                delete ref_space;
#endif
                delete matrix;
                delete rhs;
                delete solver;
                delete adaptivity;

		iter++;
	} while (!done && iter < 5);

        // Properly terminate the solver in the case of SOLVER_PETSC or SOLVER_MUMPS.
        finalize_solution_environment(matrix_solver);

	return 1;
}
