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

// error should be smaller than this epsilon
#define EPS								10e-10F


#define la0(x,y,z) (((y) + 1) / 2)
#define la1(x,y,z) (-(1 + (x) + (y) + (z)) / 2)
#define la2(x,y,z) (((x) + 1) / 2)
#define la3(x,y,z) (((z) + 1) / 2)
#define ph0 (-2.0 * 1.22474487139158904909864203735)

double fnc(double x, double y, double z) {
	return la0(x, y, z) * la1(x, y, z) * la2(x, y, z) * la3(x, y, z) * ph0 * ph0 * ph0;
}

double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
	dx = 0.75 * sqrt(3.0 / 2.0) * (((x + 1) * (y + 1) * (z + 1)) - ((y + 1) * (- z - y - x - 1) * (z + 1)));
	dy = 0.75 * sqrt(3.0 / 2.0) * (((x + 1) * (y + 1) * (z + 1)) - ((x + 1) * (- z - y - x - 1) * (z + 1)));
	dz = 0.75 * sqrt(3.0 / 2.0) * (((x + 1) * (y + 1) * (z + 1)) - ((x + 1) * (y + 1) * (- z - y - x - 1)));

	return fnc(x, y, z);
}

////

BCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *user_data) {
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename T>
T dfnc(T x, T y, T z) {
	return
		- 1.5 * sqrt(3.0 / 2.0) * (y + 1) * (z + 1)
		- 1.5 * sqrt(3.0 / 2.0) * (x + 1) * (z + 1)
		- 1.5 * sqrt(3.0 / 2.0) * (x + 1) * (y + 1);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *user_data) {
	return int_F_v<f_t, res_t>(n, wt, dfnc, u, e);
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv) {
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &argv, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	if (argc < 3) error("Not enough parameters");

	printf("* Loading mesh '%s'\n", argv[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) error("Loading mesh file '%s'\n", argv[1]);

	printf("* Setting the space up\n");
	H1ShapesetLobattoTetra ss;
	H1Space space(&mesh, &ss);
	space.set_bc_types(bc_types);

	int o;
	sscanf(argv[2], "%d", &o);
	printf("  - Setting uniform order to %d\n", o);
	space.set_uniform_order(o);

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

//	mat.dump(stdout, "a");
//	rhs.dump(stdout, "b");

	if (solved) {
		Solution sln(&mesh);
		sln.set_fe_solution(&space, solver.get_solution() );

		ExactSolution ex_sln(&mesh, exact_solution);
		// norm
//		double h1_sln_norm = h1_norm(&sln);
		double h1_err_norm = h1_error(&sln, &ex_sln);

//		printf(" - H1 solution norm:   % le\n", h1_sln_norm);
		printf(" - H1 error norm:      % le\n", h1_err_norm);

//		double l2_sln_norm = l2_norm(&sln);
//		double l2_err_norm = l2_error(&sln, &ex_sln);
//		printf(" - L2 solution norm:   % le\n", l2_sln_norm);
//		printf(" - L2 error norm:      % le\n", l2_err_norm);

//		if (h1_err_norm > EPS || l2_err_norm > EPS) {
			// calculated solution is not enough precise
//			res = ERR_FAILURE;
//		}

#ifdef AOUTPUT_DIR
		// output
		const char *of_name = OUTPUT_DIR "/solution.pos";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			DiffFilter eh(&sln, &ex_sln);
//			DiffFilter eh_dx(&sln, &ex_sln, FN_DX, FN_DX);
//			DiffFilter eh_dy(&sln, &ex_sln, FN_DY, FN_DY);
//			DiffFilter eh_dz(&sln, &ex_sln, FN_DZ, FN_DZ);

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

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	return res;
}

