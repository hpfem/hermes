// This file is part of Hermes3D
//
// Copyright (c) 2010 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

// Test to verify that hp-adaptivity worsk well.
//

#include "config.h"
#include <hermes3d.h>
#include "../../../../hermes_common/trace.h"
#include "../../../../hermes_common/common_time_period.h"
#include "../../../../hermes_common/error.h"


MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_MUMPS, SOLVER_NOX, 
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_UMFPACK.
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).

#undef REFERENCE_SOLUTION			// use ref. solution for guiding the adaptivity

const double TOLERANCE = 0.000001;		// error tolerance in percent
const double THRESHOLD = 0.3;			// error threshold for element refinement

#define ANISO_X						1
#define ANISO_Y						2
#define ANISO_Z						4

int init_q = 1;							// initial order
int init_p = 2;							// initial order
int aniso_type = 0;

// error should be smaller than this epsilon
#define EPS								10e-10F

double fnc(double x, double y, double z)
{
	switch (aniso_type) {
		case ANISO_X: return sin(x);
		case ANISO_Y: return sin(y);
		case ANISO_Z: return sin(z);

		case ANISO_X | ANISO_Y: return sin(x) * sin(y);
		case ANISO_X | ANISO_Z: return sin(x) * sin(z);
		case ANISO_Y | ANISO_Z: return sin(y) * sin(z);

		case ANISO_X | ANISO_Y | ANISO_Z: return sin(x) * sin(y) * sin(z);

		default:
			return 0;
	}
}

template<typename T>
T rhs(T x, T y, T z)
{
	T ddxx = 0;
	T ddyy = 0;
	T ddzz = 0;

	switch (aniso_type) {
		case ANISO_X: ddxx = -sin(x); break;
		case ANISO_Y: ddyy = -sin(y); break;
		case ANISO_Z: ddzz = -sin(z); break;

		case ANISO_X | ANISO_Y:
			ddxx = - sin(x) * sin(y);
			ddyy = - sin(x) * sin(y);
			break;

		case ANISO_X | ANISO_Z:
			ddxx = - sin(x) * sin(z);
			ddzz = - sin(x) * sin(z);
			break;

		case ANISO_Y | ANISO_Z:
			ddyy = - sin(y) * sin(z);
			ddzz = - sin(y) * sin(z);
			break;

		case ANISO_X | ANISO_Y | ANISO_Z:
			ddxx = - sin(x) * sin(y) * sin(z);
			ddyy = - sin(x) * sin(y) * sin(z);
			ddzz = - sin(x) * sin(y) * sin(z);
			break;
	}

	return ddxx + ddyy + ddzz;
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz)
{
	switch (aniso_type) {
		case ANISO_X:
			dx = cos(x);
			dy = 0;
			dz = 0;
			break;

		case ANISO_Y:
			dx = 0;
			dy = cos(y);
			dz = 0;
			break;

		case ANISO_Z:
			dx = 0;
			dy = 0;
			dz = cos(z);
			break;

		case ANISO_X | ANISO_Y:
			dx = cos(x) * sin(y);
			dy = sin(x) * cos(y);
			dz = 0;
			break;

		case ANISO_X | ANISO_Z:
			dx = cos(x) * sin(z);
			dy = 0;
			dz = sin(x) * cos(z);
			break;

		case ANISO_Y | ANISO_Z:
			dx = 0;
			dy = cos(y) * sin(z);
			dz = sin(y) * cos(z);
			break;

		case ANISO_X | ANISO_Y | ANISO_Z:
			dx = cos(x) * sin(y) * sin(z);
			dy = sin(x) * cos(y) * sin(z);
			dz = sin(x) * sin(y) * cos(z);
			break;
	}

	return fnc(x, y, z);
}

//

BCType bc_types(int marker)
{
	switch (aniso_type) {
		case ANISO_X:
			if (marker == 1 || marker == 2) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_Y:
			if (marker == 3 || marker == 4) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_Z:
			if (marker == 5 || marker == 6) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_X | ANISO_Y:
			if (marker == 1 || marker == 2 || marker == 3 || marker == 4) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_X | ANISO_Z:
			if (marker == 1 || marker == 2 || marker == 5 || marker == 6) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_Y | ANISO_Z:
			if (marker == 3 || marker == 4 || marker == 5 || marker == 6) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_X | ANISO_Y | ANISO_Z:
			return BC_ESSENTIAL;
	}

	return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return 0;
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, Func<res_t> *u_ext[], Func<f_t> *u, Func<f_t> *v, Geom<f_t> *e, ExtData<res_t> *data)
{
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, Func<res_t> *u_ext[], Func<f_t> *u, Geom<f_t> *e, ExtData<res_t> *data)
{
	return -int_F_v<f_t, res_t>(n, wt, rhs, u, e);
}

void out_fn_vtk(MeshFunction *x, const char *name, int i) {
#ifdef OUTPUT_DIR
	char of_name[1024];
	FILE *ofile;
	// mesh out
	sprintf(of_name, "%s/%s-%d.vtk", OUTPUT_DIR, name, i);
	ofile = fopen(of_name, "w");
	if (ofile != NULL) {
		VtkOutputEngine output(ofile);
		output.out(x, name);
		fclose(ofile);
	}
	else {
		warning("Can not open '%s' for writing.", of_name);
	}
#endif
}

void parse_aniso_type(char *str)
{
	int type = 0;
	if (strchr(str, 'x') != NULL) type |= ANISO_X;
	if (strchr(str, 'y') != NULL) type |= ANISO_Y;
	if (strchr(str, 'z') != NULL) type |= ANISO_Z;
	aniso_type = type;
}

bool check_order(const Ord3 &spord)
{
	switch (aniso_type)
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

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args)
{
	int res = ERR_SUCCESS;

	if (argc < 3) error("Not enough parameters.");

	parse_aniso_type(args[2]);

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh;
	H3DReader mloader;
	if (!mloader.load(args[1], &mesh)) error("Loading mesh file '%s'\n", args[1]);

	printf("* Setting the space up\n");

	Ord3 order;
	switch (aniso_type)
	{
		case ANISO_X:
			order = Ord3(init_p, init_q, init_q);
			break;

		case ANISO_Y:
			order = Ord3(init_q, init_p, init_q);
			break;

		case ANISO_Z:
			order = Ord3(init_q, init_q, init_p);
			break;

		case ANISO_X | ANISO_Y:
			order = Ord3(init_p, init_p, init_q);
			break;

		case ANISO_X | ANISO_Z:
			order = Ord3(init_p, init_q, init_p);
			break;

		case ANISO_Y | ANISO_Z:
			order = Ord3(init_q, init_p, init_p);
			break;

		case ANISO_X | ANISO_Y | ANISO_Z:
			order = Ord3(init_p, init_p, init_p);
			break;
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
			res = ERR_FAILURE;
			printf("failed\n");
			break;
		}

		Ord3 spord = space.get_element_order(1);
		if (!check_order(spord)) {
			res = ERR_FAILURE;
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
			res = ERR_FAILURE;
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
			res = ERR_FAILURE;
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
		if (err_est_rel < TOLERANCE || iter == 8) {
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

			if (h1_err_norm > EPS || l2_err_norm > EPS) {
				// calculated solution is not enough precise
				res = ERR_FAILURE;
			}

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

	return res;
}
