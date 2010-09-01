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

// Two uncoupled equations:
// -\laplace u_1 = f_1        u_1 = (1 - x^2) (1 - y^2) (1 - z^2)
// -\laplace u_2 = f_2        u_2 = (1 - x^2) x^2 (1 - y^2) y^2 (1 - z^2) z^2
//
// u_1 = 0 on \partial\Omega
// u_2 = 0 on \partial\Omega
//

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

// error should be smaller than this epsilon
#define EPS								10e-10F

// needed for calculation norms and used by visualizator
double exact_sln_fn_1(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = -2 * x * (1 - y*y) * (1 - z*z);
	dy = -2 * (1 - x*x) * y * (1 - z*z);
	dz = -2 * (1 - x*x) * (1 - y*y) * z;

	return (1 - x*x) * (1 - y*y) * (1 - z*z);
}

double exact_sln_fn_2(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = 2 * x * (1 - x*x) * y*y * (1 - y*y) * z*z * (1 - z*z) - 2 * x*x*x * y*y * (1 - y*y) * z*z * (1 - z*z);
	dy = 2 * x*x * (1 - x*x) * y * (1 - y*y) * z*z * (1 - z*z) - 2 * x*x * (1 - x*x) * y*y*y * z*z * (1 - z*z);
	dz = 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * z * (1 - z*z) - 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * z*z*z;

	return (1 - x*x) * x*x * (1 - y*y) * y*y * (1 - z*z) * z*z;
}

//

BCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

template<typename f_t, typename res_t>
res_t bilinear_form_1(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t bilinear_form_2(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename T>
T f1(T x, T y, T z) {
	T ddxx = -2 * (1 - y*y) * (1 - z*z);
	T ddyy = -2 * (1 - x*x) * (1 - z*z);
	T ddzz = -2 * (1 - x*x) * (1 - y*y);

	return -(ddxx + ddyy + ddzz);
}

template<typename f_t, typename res_t>
res_t linear_form_1(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_F_v<f_t, res_t>(n, wt, f1, u, e);
}

template<typename T>
T f2(T x, T y, T z) {
	T ddxx = 2 * (1 - x*x) * y*y * (1 - y*y) * z*z * (1 - z*z) - 10 * x*x * y*y * (1 - y*y) * z*z * (1 - z*z);
	T ddyy = 2 * x*x * (1 - x*x) * (1 - y*y) * z*z * (1 - z*z) - 10 * x*x * (1 - x*x) * y*y * z*z * (1 - z*z);
	T ddzz = 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * (1 - z*z) - 10 * x*x * (1 - x*x) * y*y * (1 - y*y) * z*z;

	return -(ddxx + ddyy + ddzz);
}

template<typename f_t, typename res_t>
res_t linear_form_2(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_F_v<f_t, res_t>(n, wt, f2, u, e);
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	if (argc < 2) error("Not enough parameters");

	H1ShapesetLobattoHex shapeset;

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh)) error("Loading mesh file '%s'\n", args[1]);

	printf("* Setup space #1\n");
	H1Space space1(&mesh, &shapeset);
	space1.set_bc_types(bc_types);

	order3_t o1(2, 2, 2);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o1.x, o1.y, o1.z);
	space1.set_uniform_order(o1);

	printf("* Setup space #2\n");
	H1Space space2(&mesh, &shapeset);
	space2.set_bc_types(bc_types);

	order3_t o2(4, 4, 4);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o2.x, o2.y, o2.z);
	space2.set_uniform_order(o2);

	int ndofs = 0;
	ndofs += space1.assign_dofs();
	ndofs += space2.assign_dofs(ndofs);
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

	WeakForm wf(2);
	wf.add_matrix_form(0, 0, bilinear_form_1<double, scalar>, bilinear_form_1<ord_t, ord_t>, SYM);
	wf.add_vector_form(0, linear_form_1<double, scalar>, linear_form_1<ord_t, ord_t>);

	wf.add_matrix_form(1, 1, bilinear_form_2<double, scalar>, bilinear_form_2<ord_t, ord_t>, SYM);
	wf.add_vector_form(1, linear_form_2<double, scalar>, linear_form_2<ord_t, ord_t>);

	LinearProblem lp(&wf);
	lp.set_spaces(Tuple<Space *>(&space1, &space2));

	// assemble stiffness matrix
	Timer assemble_timer("Assembling stiffness matrix");
	assemble_timer.start();
	lp.assemble(&mat, &rhs);
	assemble_timer.stop();

	// solve the stiffness matrix
	Timer solve_timer("Solving stiffness matrix");
	solve_timer.start();
	bool solved = solver.solve();
	solve_timer.stop();

	// output the measured values
	printf("%s: %s (%lf secs)\n", assemble_timer.get_name(), assemble_timer.get_human_time(), assemble_timer.get_seconds());
	printf("%s: %s (%lf secs)\n", solve_timer.get_name(), solve_timer.get_human_time(), solve_timer.get_seconds());

	if (solved) {
		// solution 1
		Solution sln1(&mesh);
		sln1.set_fe_solution(&space1, solver.get_solution());

		ExactSolution esln1(&mesh, exact_sln_fn_1);
		// norm
		double h1_sln_norm1 = h1_norm(&sln1);
		double h1_err_norm1 = h1_error(&sln1, &esln1);

		printf(" - H1 solution norm:   % le\n", h1_sln_norm1);
		printf(" - H1 error norm:      % le\n", h1_err_norm1);

		double l2_sln_norm1 = l2_norm(&sln1);
		double l2_err_norm1 = l2_error(&sln1, &esln1);
		printf(" - L2 solution norm:   % le\n", l2_sln_norm1);
		printf(" - L2 error norm:      % le\n", l2_err_norm1);

		if (h1_err_norm1 > EPS || l2_err_norm1 > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}

		// solution 2
		Solution sln2(&mesh);
		sln2.set_fe_solution(&space2, solver.get_solution());

		ExactSolution esln2(&mesh, exact_sln_fn_2);
		// norm
		double h1_sln_norm2 = h1_norm(&sln2);
		double h1_err_norm2 = h1_error(&sln2, &esln2);

		printf(" - H1 solution norm:   % le\n", h1_sln_norm2);
		printf(" - H1 error norm:      % le\n", h1_err_norm2);

		double l2_sln_norm2 = l2_norm(&sln2);
		double l2_err_norm2 = l2_error(&sln2, &esln2);
		printf(" - L2 solution norm:   % le\n", l2_sln_norm2);
		printf(" - L2 error norm:      % le\n", l2_err_norm2);

		if (h1_err_norm2 > EPS || l2_err_norm2 > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}

#ifdef OUTPUT_DIR
		// output
		const char *of_name = OUTPUT_DIR "/solution.pos";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			GmshOutputEngine output(ofile);
			output.out(&sln1, "Uh_1");
			output.out(&esln1, "U1");
			output.out(&sln2, "Uh_2");
			output.out(&esln2, "U2");

			fclose(ofile);
		}
		else {
			warning("Can not open '%s' for writing.", of_name);
		}
#endif
	}

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

