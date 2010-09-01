// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
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

// Test to verify the hp-adaptivity
//
// Starts with linear elements and should find the solution (just using p refinements)
//

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

const double TOLERANCE = 0.001;		// error tolerance in percent
const double THRESHOLD = 0.5;		// error threshold for element refinement

int m = 2;
int n = 2;
int o = 2;

// error should be smaller than this epsilon
#define EPS								10e-10F

double fnc(double x, double y, double z) {
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
}

template<typename T>
T dfnc(T x, T y, T z) {
	T ddxx = m * (m - 1) * pow(x, m - 2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
	T ddyy = n * (n - 1) * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = o * (o - 1) * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);

	return -(ddxx + ddyy + ddzz);
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	// u(x, y, z) = pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4)
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 3) + 4 * pow(z, 3);

	return fnc(x, y, z);
}

//

BCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z) {
	return fnc(x, y, z);
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_F_v<f_t, res_t>(n, wt, dfnc, u, e);
}


// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &argv, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	if (argc < 5) error("Not enough parameters.");

	H1ShapesetLobattoHex shapeset;

	printf("* Loading mesh '%s'\n", argv[1]);
	Mesh mesh;
	Mesh3DReader mloader;
	if (!mloader.load(argv[1], &mesh)) error("Loading mesh file '%s'\n", argv[1]);

	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_essential_bc_values(essential_bc_values);

	int o[3] = { 0, 0, 0 };
	sscanf(argv[2], "%d", o + 0);
	sscanf(argv[3], "%d", o + 1);
	sscanf(argv[4], "%d", o + 2);
	order3_t order(o[0], o[1], o[2]);
	printf("  - Setting uniform order to (%d, %d, %d)\n", order.x, order.y, order.z);
	space.set_uniform_order(order);

	WeakForm wf;
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM, ANY);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>, ANY);

	LinearProblem lp(&wf);
	lp.set_space(&space);

	bool done = false;
	int iter = 0;
	do {
		Timer assemble_timer("Assembling stiffness matrix");
		Timer solve_timer("Solving stiffness matrix");

		printf("\n=== Iter #%d ================================================================\n", iter);

		printf("\nSolution\n");

#if defined WITH_UMFPACK
		UMFPackMatrix mat;
		UMFPackVector rhs;
		UMFPackLinearSolver solver(&mat, &rhs);
#elif defined WITH_PARDISO
		PardisoMatrix mat;
		PardisoVector rhs;
		PardisoLinearSolver solver(&mat, &rhs);
#elif defined WITH_PETSC
		PetscMatrix mat;
		PetscVector rhs;
		PetscLinearSolver solver(&mat, &rhs);
#elif defined WITH_MUMPS
		MumpsMatrix mat;
		MumpsVector rhs;
		MumpsSolver solver(&mat, &rhs);
#endif

		int ndofs = space.assign_dofs();
		printf("  - Number of DOFs: %d\n", ndofs);

		// assemble stiffness matrix
		printf("  - Assembling... "); fflush(stdout);
		assemble_timer.reset();
		assemble_timer.start();
		lp.assemble(&mat, &rhs);
		assemble_timer.stop();
		printf("done in %s (%lf secs)\n", assemble_timer.get_human_time(), assemble_timer.get_seconds());

		// solve the stiffness matrix
		printf("  - Solving... "); fflush(stdout);
		solve_timer.reset();
		solve_timer.start();
		bool solved = solver.solve();
		solve_timer.stop();
		if (solved)
			printf("done in %s (%lf secs)\n", solve_timer.get_human_time(), solve_timer.get_seconds());
		else {
			res = ERR_FAILURE;
			printf("failed\n");
			break;
		}

		printf("Reference solution\n");

#if defined WITH_UMFPACK
		UMFPackLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_PARDISO
		PardisoLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_PETSC
		PetscLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_MUMPS
		MumpsSolver rsolver(&mat, &rhs);
#endif

		Mesh rmesh;
		rmesh.copy(mesh);
		rmesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

		Space *rspace = space.dup(&rmesh);
		rspace->copy_orders(space, 1);

		LinearProblem rlp(&wf);
		rlp.set_space(rspace);

		int rndofs = rspace->assign_dofs();
		printf("  - Number of DOFs: %d\n", rndofs);

		printf("  - Assembling... "); fflush(stdout);
		assemble_timer.reset();
		assemble_timer.start();
		rlp.assemble(&mat, &rhs);
		assemble_timer.stop();
		printf("done in %s (%lf secs)\n", assemble_timer.get_human_time(), assemble_timer.get_seconds());

		printf("  - Solving... "); fflush(stdout);
		solve_timer.reset();
		solve_timer.start();
		bool rsolved = rsolver.solve();
		solve_timer.stop();
		if (rsolved)
			printf("done in %s (%lf secs)\n", solve_timer.get_human_time(), solve_timer.get_seconds());
		else {
			res = ERR_FAILURE;
			printf("failed\n");
			break;
		}

		Solution sln(&mesh);
		sln.set_fe_solution(&space, solver.get_solution());

		Solution rsln(&rmesh);
		rsln.set_fe_solution(rspace, rsolver.get_solution());

		printf("Adaptivity:\n");
		H1Adapt hp(&space);
		double tol = hp.calc_error(&sln, &rsln) * 100;
		printf("  - tolerance: "); fflush(stdout);
		printf("% lf\n", tol);
		if (tol < TOLERANCE) {
			printf("\nDone\n");
			ExactSolution ex_sln(&mesh, exact_solution);
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

		Timer t("");
		printf("  - adapting... "); fflush(stdout);
		t.start();
		hp.adapt(THRESHOLD);
		t.stop();
		printf("done in %lf secs (refined %d element(s))\n", t.get_seconds(), hp.get_num_refined_elements());

		iter++;
	} while (!done);


#ifdef WITH_PETSC
	PetscFinalize();
#endif

	return res;
}

