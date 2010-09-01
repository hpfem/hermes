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

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

// first two Lobatto shape functions
#define l0(x) ((1.0 - (x)) * 0.5)
#define l1(x) ((1.0 + (x)) * 0.5)


// error should be smaller than this epsilon
#define EPS								10e-10F

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = -0.5 * x * l0(y) * l1(y) * l0(z) * l1(z);
	dy = -0.5 * y * l0(x) * l1(x) * l0(z) * l1(z);
	dz = -0.5 * z * l0(x) * l1(x) * l0(y) * l1(y);

	return l0(x) * l1(x) * l0(y) * l1(y) * l0(z) * l1(z);
}

//

BCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename T>
T f(T x, T y, T z) {
	return
		0.5 * l0(y) * l1(y) * l0(z) * l1(z) +
		0.5 * l0(x) * l1(x) * l0(z) * l1(z) +
		0.5 * l0(x) * l1(x) * l0(y) * l1(y);
}


template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_F_v<f_t, res_t>(n, wt, f, u, e);
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	_F_

	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	if (argc < 3) error("Not enough parameters");

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh;
	Mesh3DReader mloader;
	if (!mloader.load(args[1], &mesh)) error("Loading mesh file '%s'\n", args[1]);

	FOR_ALL_ELEMENTS(idx, &mesh) {
		Element *e = mesh.elements[idx];

		printf("elem = %ld:", idx);
		for (int iface = 0; iface < 6; iface++) {
			printf(" %d,", e->get_face_orientation(iface));
		}
		printf("\n");
	}

	H1ShapesetLobattoHex shapeset;
	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);

	int o, p, q;
	sscanf(args[2], "%d", &o);
	if (argc > 3) sscanf(args[3], "%d", &p);
	else p = o;
	if (argc > 4) sscanf(args[4], "%d", &q);
	else q = o;
	order3_t order(o, p, q);
	printf("  - Setting uniform order to %s\n", order.str());
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
	wf.add_matrix_form(FORM_CB(bilinear_form), SYM);
	wf.add_vector_form(FORM_CB(linear_form));

	LinearProblem lp(&wf);
	lp.set_space(&space);

	printf("  - assembling... ");
	Timer assemble_timer;
	assemble_timer.start();
	lp.assemble(&mat, &rhs);
	assemble_timer.stop();
	printf("%s (%lf secs)\n", assemble_timer.get_human_time(), assemble_timer.get_seconds());

//	mat.dump(stdout, "a");
//	rhs.dump(stdout, "b");

	// solve the stiffness matrix
	printf("  - solving... ");
	Timer solve_timer;
	solve_timer.start();
	bool solved = solver.solve();
	solve_timer.stop();

	// output the measured values
	printf("%s (%lf secs)\n", solve_timer.get_human_time(), solve_timer.get_seconds());

	if (solved) {
		Solution sln(&mesh);
		sln.set_fe_solution(&space, solver.get_solution());

		printf("* Solution:\n");

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

#ifdef OUTPUT_DIR
		// output
		const char *of_name = OUTPUT_DIR "/solution.pos";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			ExactSolution ex_sln(&mesh, exact_solution);
//			DiffFilter eh(&sln, &ex_sln);
//			DiffFilter eh_dx(&mesh, &sln, &ex_sln, FN_DX, FN_DX);
//			DiffFilter eh_dy(&mesh, &sln, &ex_sln, FN_DY, FN_DY);
//			DiffFilter eh_dz(&mesh, &sln, &ex_sln, FN_DZ, FN_DZ);

			GmshOutputEngine output(ofile);
			output.out(&sln, "Uh");
//			output.out(&sln, "Uh dx", FN_DX_0);
//			output.out(&sln, "Uh dy", FN_DY_0);
//			output.out(&sln, "Uh dz", FN_DZ_0);
//			output.out(&eh, "Eh");
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

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	return res;
}

