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

int m = 2;
int n = 2;
int o = 2;

double fnc(double x, double y, double z)
{
	_F_
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
}

// needed for calculation norms and used by visualizator
double exact_sln_fn(double x, double y, double z, double &dx, double &dy, double &dz)
{
	_F_
	dx = m * pow(x, m - 1) * pow(y, n) * pow(z, o) - 3 * pow(x, 2) * z + 2 * x * pow(y, 3);
	dy = n * pow(x, m) * pow(y, n - 1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o - 1) + 4 * pow(z, 3) - pow(x, 3);

	return fnc(x, y, z);
}

//

BCType bc_types(int marker)
{
	_F_
	return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	_F_
	return fnc(x, y, z);
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                    user_data_t<res_t> *data)
{
	_F_
	return
		int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e) +
		int_dudx_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename T>
T f(T x, T y, T z)
{
	_F_
	T ddxx = (m - 1) * m * pow(x, m - 2) * pow(y, n) * pow(z, o) - 6 * x * z + 2 * pow(y, 3);
	T ddyy = (n - 1) * n * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = (o - 1) * o * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);
	T dx = m * pow(x, m - 1) * pow(y, n) * pow(z, o) - 3 * pow(x, 2) * z + 2 * x * pow(y, 3);

	return -(ddxx + ddyy + ddzz) + dx;
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data)
{
	_F_
	return int_F_v<f_t, res_t>(n, wt, f, u, e);
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	_F_
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &argv, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	if (argc < 2) error("Not enough parameters.");

	H1ShapesetLobattoHex shapeset;

	printf("* Loading mesh '%s'\n", argv[1]);
	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(argv[1], &mesh)) error("loading mesh file '%s' failed.", argv[1]);

	printf("* Setup space\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_essential_bc_values(essential_bc_values);

	order3_t o(4, 4, 4);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o.x, o.y, o.z);
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
	wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>);

	LinearProblem lp(&wf);
	lp.set_space(&space);

	// assemble stiffness matrix
	printf("  - assembling... "); fflush(stdout);
	Timer tmr_assemble("");
	tmr_assemble.start();
	lp.assemble(&mat, &rhs);
	tmr_assemble.stop();
	printf("done in %s (%lf secs)\n", tmr_assemble.get_human_time(), tmr_assemble.get_seconds());

	// solve the stiffness matrix
	printf("  - solving... "); fflush(stdout);
	Timer tmr_solve("");
	tmr_solve.start();
	bool solved = solver.solve();
	tmr_solve.stop();
	printf("done in %s (%lf secs)\n", tmr_solve.get_human_time(), tmr_solve.get_seconds());

	if (solved) {
		// solution
		Solution sln(&mesh);
		sln.set_fe_solution(&space, solver.get_solution());

		ExactSolution esln(&mesh, exact_sln_fn);
		// norm
		double h1_sln_norm = h1_norm(&sln);
		double h1_err_norm = h1_error(&sln, &esln);

		printf("  - H1 solution norm:   % le\n", h1_sln_norm);
		printf("  - H1 error norm:      % le\n", h1_err_norm);

		double l2_sln_norm = l2_norm(&sln);
		double l2_err_norm = l2_error(&sln, &esln);
		printf("  - L2 solution norm:   % le\n", l2_sln_norm);
		printf("  - L2 error norm:      % le\n", l2_err_norm);

		if (h1_err_norm > EPS || l2_err_norm > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}

#ifdef OUTPUT_DIR
		// output
		const char *of_name = OUTPUT_DIR "/solution.pos";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			GmshOutputEngine output(ofile);
			output.out(&sln, "Uh");
			output.out(&esln, "U");

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
