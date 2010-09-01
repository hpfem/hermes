// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - Pavel Kus
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

//
// Time harmonic Maxwell equations with impedance BC
//
// calculates solution of the problem:
//  curl(curl(u)) = f on Omega
//  u x n = 0 on Gamma_P
//  curl(u) x n - i (n x u) x n = g on Gamma_I
//  where n is outer normal
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

std::complex<double> img(0, 1);

// general polynomial function

scalar3 &exact_solution(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz) {
	dx[0] = 3.*x*x*y*y - 3.*y*y*y*z;
	dy[0] = 2.*x*x*x*y - 9.*x*y*y*z;
	dz[0] = -3.*x*y*y*y;

	dx[1] = 3.*x*x*y*y*y*z*z*z + 4.*x*y*z;
	dy[1] = 3.*x*x*x*y*y*z*z*z + 2.*x*x*z;
	dz[1] = 3.*x*x*x*y*y*y*z*z + 2.*x*x*y;

	dx[2] = -12.*x*x;
	dy[2] = z*z*z;
	dz[2] = 3.*y*z*z;

	static scalar3 val;
	val[0] = x*x*x*y*y - 3.*x*y*y*y*z;
	val[1] = x*x*x*y*y*y*z*z*z + 2.*x*x*y*z;
	val[2] = y*z*z*z - 4.*x*x*x;
	return val;
}

template<typename S, typename T>
void exact_sln(S x, S y, S z, T (&fn)[3], T (&dx)[3], T (&dy)[3], T (&dz)[3]) {
	fn[0] = x*x*x*y*y - 3.*x*y*y*y*z;
	fn[1] = x*x*x*y*y*y*z*z*z + 2.*x*x*y*z;
	fn[2] = y*z*z*z - 4.*x*x*x;

	dx[0] = 3.*x*x*y*y - 3.*y*y*y*z;
	dx[1] = 3.*x*x*y*y*y*z*z*z + 4.*x*y*z;
	dx[2] = -12.*x*x;

	dy[0] = 2.*x*x*x*y - 9.*x*y*y*z;
	dy[1] = 3.*x*x*x*y*y*z*z*z + 2.*x*x*z;
	dy[2] = z*z*z;

	dz[0] = -3.*x*y*y*y;
	dz[1] = 3.*x*x*x*y*y*y*z*z + 2.*x*x*y;
	dz[2] = 3.*y*z*z;
}

template<typename ct, typename res_t>
void f(ct x, ct y, ct z, res_t (&val)[3]) {
	res_t ev[3], dx[3], dy[3], dz[3];
	exact_sln(x, y, z, ev, dx, dy, dz);

	val[0] = 4*x*z + 18*x*y*z - 2*x*x*x + 9*x*x*y*y*z*z*z - ev[0];
	val[1] = -4*y*z + 3*z*z - 9*z*y*y + 6*y*x*x - 6*x*y*y*y*z*z*z - 6*z*x*x*x*y*y*y - ev[1];
	val[2] = 24*x + 2*x*x + 9*x*x*x*y*y*z*z - 3*y*y*y - ev[2];
}

/*
// TODO: this could be written in a much simpler way. Just use curl of exact solution
// and cross product defined in Scalar3D...
scalar3 &bc_values(int ess_bdy_marker, double x, double y, double z) {
	static scalar bc[3] = { 0., 0., 0. };

	switch (marker) {
		case 1:
			bc[1] = -4*x*y*z - 9*x*z*y*y + 2*y*x*x*x - 3*x*x*y*y*y*z*z*z;
			bc[2] = 12*x*x - 3*x*y*y*y;
			break;

		case 2:
			bc[1] = 4*x*y*z + 9*x*z*y*y - 2*y*x*x*x + 3*x*x*y*y*y*z*z*z;
			bc[2] = -12*x*x + 3*x*y*y*y;
			break;

		case 3:
			bc[0] = 4*x*y*z + 9*x*z*y*y - 2*y*x*x*x + 3*x*x*y*y*y*z*z*z;
			bc[2] = 2*y*x*x + 3*x*x*x*y*y*y*z*z - z*z*z;
			break;

		case 4:
			bc[0] = -4*x*y*z - 9*x*z*y*y + 2*y*x*x*x - 3*x*x*y*y*y*z*z*z;
			bc[2] = -2*y*x*x - 3*x*x*x*y*y*y*z*z + z*z*z;
			break;

		case 5:
			bc[0] = -12*x*x + 3*x*y*y*y;
			bc[1] = -2*y*x*x - 3*x*x*x*y*y*y*z*z + z*z*z;
			break;

		case 6:
			bc[0] = 12*x*x - 3*x*y*y*y;
			bc[1] = 2*y*x*x + 3*x*x*x*y*y*y*z*z - z*z*z;
			break;

		default:
			EXIT(H3D_ERR_FACE_INDEX_OUT_OF_RANGE);
	}

	switch (marker) {
		case 1:
		case 2:
			bc[1] -= img * (2*y*z*x*x + x*x*x*y*y*y*z*z*z);
			bc[2] -= img * (-4*x*x*x + y*z*z*z);
			break;

		case 3:
		case 4:
			bc[0] -= img * (x*x*x*y*y - 3*x*z*y*y*y);
			bc[2] -= img * (-4*x*x*x + y*z*z*z);
			break;

		case 5:
		case 6:
			bc[0] -= img * (x*x*x*y*y - 3*x*z*y*y*y);
			bc[1] -= img * (2*y*z*x*x + x*x*x*y*y*y*z*z*z);
			break;

		default:
			EXIT(H3D_ERR_FACE_INDEX_OUT_OF_RANGE);
	}

	return bc;
}
*/

BCType bc_types(int marker) {
	return BC_NATURAL;
}

/// definition of the forms

template<typename ct, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<ct> *u, fn_t<ct> *v, geom_t<ct> *e, user_data_t<res_t> *data) {
	return
		hcurl_int_curl_u_curl_v<ct, res_t>(n, wt, u, v, e) -
		hcurl_int_u_v<ct, res_t>(n, wt, u, v, e);
}

template<typename ct, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<ct> *v, geom_t<ct> *e, user_data_t<res_t> *data) {
	return hcurl_int_F_v<ct, res_t>(n, wt, f, v, e);
}

template<typename ct, typename res_t>
res_t bilinear_form_surf(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<ct> *u, fn_t<ct> *v, geom_t<ct> *e, user_data_t<res_t> *data) {
	return -img * hcurl_int_u_v<ct, res_t>(n, wt, u, v, e);
}

template<typename ct, typename res_t>
res_t linear_form_surf(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<ct> *v, geom_t<ct> *e, user_data_t<res_t> *data) {
	res_t result = 0.0;
	for (int i = 0; i < n; i++) {
		res_t ev[3], dx[3], dy[3], dz[3];
		exact_sln(e->x[i], e->y[i], e->z[i], ev, dx, dy, dz);

		res_t curl_e[3];
		calc_curl(dx, dy, dz, curl_e);
		res_t tpe[3];
		calc_tan_proj(e->nx[i], e->ny[i], e->nz[i], ev, tpe);

		res_t g[3] = {
			(e->nz[i] * curl_e[1] - e->ny[i] * curl_e[2]) - img * tpe[0],
			(e->nx[i] * curl_e[2] - e->nz[i] * curl_e[0]) - img * tpe[1],
			(e->ny[i] * curl_e[0] - e->nx[i] * curl_e[1]) - img * tpe[2],
		};
		result += wt[i] * (v->fn0[i] * g[0] + v->fn1[i] * g[1] + v->fn2[i] * g[2]);
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

	if (argc < 3) error("Not enough parameters");

	HcurlShapesetLobattoHex shapeset;

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh)) error("Loading mesh file '%s'\n", args[1]);

	printf("* Setting the space up\n");
	HcurlSpace space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
//	space.set_essential_bc_values(essential_bc_values);

	int order;
	sscanf(args[2], "%d", &order);
	int dir_x = order, dir_y = order, dir_z = order;
	order3_t o(dir_x, dir_y, dir_z);
	printf("  - Setting uniform order to (%d, %d, %d)\n", dir_x, dir_y, dir_z);
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
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, UNSYM);
	wf.add_matrix_form_surf(bilinear_form_surf<double, scalar>, bilinear_form_surf<ord_t, ord_t>);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>);
	wf.add_vector_form_surf(linear_form_surf<double, scalar>, linear_form_surf<ord_t, ord_t>);
//	wf.add_biform(bilinear_form, SYM);
//	wf.add_biform_surf(bilinear_form_surf);
//	wf.add_liform(linear_form);
//	wf.add_liform_surf(linear_form_surf);

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

//	mat.dump(stdout, "a");
//	rhs.dump(stdout, "b");

	if (solved) {
		Solution sln(&mesh);
		sln.set_fe_solution(&space, solver.get_solution());

		printf("* Solution:\n");
//		scalar *s = solver.get_solution();
//		for (int i = 0; i < ndofs; i++) {
//			printf(" x[% 3d] = " SCALAR_FMT "\n", i, SCALAR(s[i]));
//		}

		// output the measured values
		printf("%s: %s (%lf secs)\n", assemble_timer.get_name(), assemble_timer.get_human_time(), assemble_timer.get_seconds());
		printf("%s: %s (%lf secs)\n", solve_timer.get_name(), solve_timer.get_human_time(), solve_timer.get_seconds());

		// norm
		ExactSolution ex_sln(&mesh, exact_solution);
		double hcurl_sln_norm = hcurl_norm(&sln);
		double hcurl_err_norm = hcurl_error(&sln, &ex_sln);
		printf(" - Hcurl solution norm: % le\n", hcurl_sln_norm);
		printf(" - Hcurl error norm:    % le\n", hcurl_err_norm);

		double l2_sln_norm = l2_norm_hcurl(&sln);
		double l2_err_norm = l2_error_hcurl(&sln, &ex_sln);
		printf(" - L2 solution norm:    % le\n", l2_sln_norm);
		printf(" - L2 error norm:       % le\n", l2_err_norm);

		if (hcurl_err_norm > EPS || l2_err_norm > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}


#if 0 //def OUTPUT_DIR
		// output
		printf("starting output\n");
		const char *of_name = OUTPUT_DIR "/solution.pos";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			ExactSolution ex_sln(&mesh, exact_solution_0, exact_solution_1, exact_solution_2);

			RealPartFilter real_sln(&mesh, &sln, FN_VAL);
			ImagPartFilter imag_sln(&mesh, &sln, FN_VAL);

			GmshOutputEngine output(ofile);

			output.out(&real_sln, "real_Uh", FN_VAL);
			output.out(&imag_sln, "imag_Uh", FN_VAL);

			output.out(&real_sln, "real_Uh_0", FN_VAL_0);
			output.out(&real_sln, "real_Uh_1", FN_VAL_1);
			output.out(&real_sln, "real_Uh_2", FN_VAL_2);

			output.out(&imag_sln, "imag_Uh_0", FN_VAL_0);
			output.out(&imag_sln, "imag_Uh_1", FN_VAL_1);
			output.out(&imag_sln, "imag_Uh_2", FN_VAL_2);

			DiffFilter eh(&mesh, &sln, &ex_sln);
//			DiffFilter eh_dx(mesh, &sln, &ex_sln, FN_DX, FN_DX);
//			DiffFilter eh_dy(mesh, &sln, &ex_sln, FN_DY, FN_DY);
//			DiffFilter eh_dz(mesh, &sln, &ex_sln, FN_DZ, FN_DZ);

			RealPartFilter real_eh(&mesh, &eh, FN_VAL);
			ImagPartFilter imag_eh(&mesh, &eh, FN_VAL);

			output.out(&real_eh, "real_Eh", FN_VAL);
			output.out(&imag_eh, "imag_Eh", FN_VAL);

			fclose(ofile);
		}
		else {
			ERROR("Can not open '%s' for writing.", of_name);
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

