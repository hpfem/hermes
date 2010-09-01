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

/*
 * hex-neumann.cc
 *
 * exhaustive testing of HERMES3D with neumann boundary condition
 *
 * -\Delta u + u = f in Omega
 * du/dn = g on dOmega
 *
 * u(x,y,z) = x^m * y^n * z^o + x^2 * y^3 - x^3 * z + z^4
 */

#include "config.h"
#include <math.h>
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>


// error should be smaller than this epsilon
#define EPS								10e-10F

int m, n, o;

template<typename T>
T fnc(T x, T y, T z) {
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
}

template<typename T>
T dfnc(T x, T y, T z) {
	T ddxx = m*(m-1) * pow(x, m-2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
	T ddyy = n*(n-1) * pow(x, m) * pow(y, n-2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = o*(o-1) * pow(x, m) * pow(y, n) * pow(z, o-2) + 12 * pow(z, 2);

	return -(ddxx + ddyy + ddzz) + fnc(x, y, z);
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 3) + 4 * pow(z, 3);

	return fnc(x, y, z);
}

//

BCType bc_types(int marker) {
	return BC_NATURAL;
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return
		int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e) +
		int_u_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_F_v<f_t, res_t>(n, wt, dfnc, u, e);
}

template<typename f_t, typename res_t>
res_t linear_form_surf(int np, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data) {
	res_t result = 0;
	for (int i = 0; i < np; i++) {
		res_t dx = m * pow(e->x[i], m-1) * pow(e->y[i], n) * pow(e->z[i], o) + 2 * e->x[i] * pow(e->y[i], 3) - 3 * pow(e->x[i], 2) * e->z[i];
		res_t dy = n * pow(e->x[i], m) * pow(e->y[i], n-1) * pow(e->z[i], o) + 3 * pow(e->x[i], 2) * pow(e->y[i], 2);
		res_t dz = o * pow(e->x[i], m) * pow(e->y[i], n) * pow(e->z[i], o-1) - pow(e->x[i], 3) + 4 * pow(e->z[i], 3);

		result += wt[i] * (u->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i]));
	}
	return result;
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

	if (argc < 5) error("Not enough parameters");

	sscanf(args[2], "%d", &m);
	sscanf(args[3], "%d", &n);
	sscanf(args[4], "%d", &o);

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh;
	Mesh3DReader mloader;
	if (!mloader.load(args[1], &mesh)) error("Loading mesh file '%s'\n", args[1]);

	H1ShapesetLobattoHex shapeset;
	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);

	int mx = maxn(4, m, n, o, 4);
	order3_t order(mx, mx, mx);
//	order3_t order(1, 1, 1);
//	order3_t order(m, n, o);
	printf("  - Setting uniform order to (%d, %d, %d)\n", mx, mx, mx);
	space.set_uniform_order(order);

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	printf("* Calculating a solution\n");

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

	WeakForm wf;
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>);
	wf.add_vector_form_surf(linear_form_surf<double, scalar>, linear_form_surf<ord_t, ord_t>);

	LinearProblem lp(&wf);
	lp.set_space(&space);

	// assemble stiffness matrix
	printf("  - assembling...\n"); fflush(stdout);
	Timer assemble_timer;
	assemble_timer.start();
	lp.assemble(&mat, &rhs);
	assemble_timer.stop();
	printf("%s (%lf secs)\n", assemble_timer.get_human_time(), assemble_timer.get_seconds());

	// solve the stiffness matrix
	printf("  - solving... "); fflush(stdout);
	Timer solve_timer;
	solve_timer.start();
	bool solved = solver.solve();
	solve_timer.stop();
	printf("%s (%lf secs)\n", solve_timer.get_human_time(), solve_timer.get_seconds());

//	mat.dump(stdout, "a");
//	rhs.dump(stdout, "b");

	if (solved) {
		Solution sln(&mesh);
		sln.set_fe_solution(&space, solver.get_solution());

//		printf("* Solution:\n");
//		double *s = solver.get_solution();
//		for (int i = 1; i <= ndofs; i++) {
//			printf(" x[% 3d] = % lf\n", i, s[i]);
//		}

		ExactSolution ex_sln(&mesh, exact_solution);
		// norm
		double h1_sln_norm = h1_norm(&sln);
		double h1_err_norm = h1_error(&sln, &ex_sln);
		printf(" - H1 solution norm:   % le\n", h1_sln_norm);
		printf(" - H1 error norm:      % le\n", h1_err_norm);

		double l2_sln_norm = l2_norm(&sln);
		double l2_err_norm = l2_error(&sln, &ex_sln);
		printf(" - L2 solution norm:   % le\n", l2_sln_norm);
		printf(" - L2 error norm:      % le\n", l2_err_norm);

		if (h1_err_norm > EPS || l2_err_norm > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}

#if 0 //def OUTPUT_DIR
		printf("* Output\n");
		// output
		const char *of_name = OUTPUT_DIR "/solution.pos";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			ExactSolution ex_sln(&mesh, exact_solution);
			DiffFilter eh(&sln, &ex_sln);
//			DiffFilter eh_dx(&mesh, &sln, &ex_sln, FN_DX, FN_DX);
//			DiffFilter eh_dy(&mesh, &sln, &ex_sln, FN_DY, FN_DY);
//			DiffFilter eh_dz(&mesh, &sln, &ex_sln, FN_DZ, FN_DZ);

			GmshOutputEngine output(ofile);
			output.out(&sln, "Uh");
//			output.out(&sln, "Uh dx", FN_DX_0);
//			output.out(&sln, "Uh dy", FN_DY_0);
//			output.out(&sln, "Uh dz", FN_DZ_0);
			output.out(&eh, "Eh");
//			output.out(&eh_dx, "Eh dx");
//			output.out(&eh_dy, "Eh dy");
//			output.out(&eh_dz, "Eh dz");
			output.out(&ex_sln, "U");
//			output.out(&ex_sln, "U dx", FN_DX_0);
//			output.out(&ex_sln, "U dy", FN_DY_0);
//			output.out(&ex_sln, "U dz", FN_DZ_0);

			fclose(ofile);
		}
		else {
			warning("Can not open '%s' for writing.", of_name);
		}
#endif
	}
	else
		res = ERR_FAILURE;

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

