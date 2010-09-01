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
// Testing projections
//
//

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>
#include <float.h>

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1

#define EPS											10e-10

//#define X2_Y2_Z2_025
//#define X2_Y2_Z2
//#define X3_Y3_Z3
#define XN_YM_ZO


int m = 2, n = 2, o = 2;

double fnc(double x, double y, double z) {
#if defined X2_Y2_Z2_025
	return pow(x*x + y*y + z*z, .25);
#elif defined X2_Y2_Z2
	return x*x + y*y + z*z;
#elif defined X3_Y3_Z3
	return x*x*x + y*y*y + z*z*z;
#elif defined XN_YM_ZO
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 2) * z + pow(z, 4);
#endif
}

template<typename T>
T dfnc(T x, T y, T z) {
#if defined X2_Y2_Z2_025
	if ((x*x + y*y + z*z) == 0.0)
		return -DBL_MAX;
	else
		return -0.75 * pow(x*x + y*y + z*z, -0.75);
#elif defined X3_Y3_Z3
	return -6.0 * x - 6.0 * y - 6.0 * z;
#elif defined X2_Y2_Z2
	return -6.0;
#elif defined XN_YM_ZO
	T ddxx = m*(m-1) * pow(x, m-2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 2 * z;
	T ddyy = n*(n-1) * pow(x, m) * pow(y, n-2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = o*(o-1) * pow(x, m) * pow(y, n) * pow(z, o-2) + 12 * pow(z, 2);
	return -(ddxx + ddyy + ddzz);
#endif
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
#if defined X2_Y2_Z2_025
	if ((x*x + y*y + z*z) != 0.0) {
		dx = 0.5 * x * pow(x*x + y*y + z*z, -.75);
		dy = 0.5 * y * pow(x*x + y*y + z*z, -.75);
		dz = 0.5 * z * pow(x*x + y*y + z*z, -.75);
	}
	else {
		// pow(x*x + y*y + z*z, -.75) is not defined
		dx = 0.0;
		dy = 0.0;
		dz = 0.0;
	}
#elif defined X2_Y2_Z2
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;
#elif defined X3_Y3_Z3
	dx = 3 * x*x;
	dy = 3 * y*y;
	dz = 3 * z*z;
#elif defined XN_YM_ZO
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 2 * x * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 2) + 4 * pow(z, 3);
#endif

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


//

#define CHECK_ERROR \
	printf("Elem #%ld: error = % lf\n", e->id, error); \
	if (error > EPS) { \
		ret = ERR_FAILURE; \
		break; \
	}

//
// main
//

int main(int argc, char *argv[])
{
	_F_
	int ret = ERROR_SUCCESS;

	if (argc < 3) {
		fprintf(stderr, "ERROR: not enough parameters\n");
		return ERR_FAILURE;
	}

	if (strcmp(argv[1], "h1") != 0 && strcmp(argv[1], "h1-ipol")) {
		fprintf(stderr, "ERROR: unknown type of the projection\n");
		return ERR_FAILURE;
	}

#ifdef WITH_PETSC
	PetscInitialize(NULL, NULL, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	H1ShapesetLobattoHex shapeset;

	Mesh mesh;
	Mesh3DReader mloader;
	if (!mloader.load(argv[2], &mesh)) {
		fprintf(stderr, "ERROR: loading mesh file '%s'\n", argv[2]);
		return ERR_FAILURE;
	}

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

	Word_t ne = mesh.elements.count();
	// make the mesh for the ref. solution
	mesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_essential_bc_values(essential_bc_values);

#if defined X2_Y2_Z2
	order3_t o(2, 2, 2);
#elif defined X3_Y3_Z3
	order3_t o(3, 3, 3);
#elif defined XN_YM_ZO
	order3_t o(2, 3, 4);
#endif
	space.set_uniform_order(o);

	WeakForm wf;
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM, ANY);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>, ANY);

	LinearProblem lp(&wf);
	lp.set_space(&space);

	space.assign_dofs();

	// assemble the stiffness matrix
	lp.assemble(&mat, &rhs);

	// solve the stiffness matrix
	solver.solve();

	Solution sln(&mesh);
	sln.set_fe_solution(&space, solver.get_solution());

	for (Word_t idx = mesh.elements.first(); idx <= ne; idx = mesh.elements.next(idx)) {
		Element *e = mesh.elements[idx];

		order3_t order(4, 4, 4);
		double error;

		Projection *proj;
		if (strcmp(argv[1], "h1") == 0) proj = new H1Projection(&sln, e, &shapeset);
		else if (strcmp(argv[1], "h1-ipol") == 0) proj = new H1ProjectionIpol(&sln, e, &shapeset);
		else return ERR_FAILURE;

		//
		error = 0.0;
		error += proj->get_error(H3D_REFT_HEX_NONE, -1, order);
		error = sqrt(error);
		CHECK_ERROR;

		//
		error = 0.0;
		error += proj->get_error(H3D_REFT_HEX_X, 20, order);
		error += proj->get_error(H3D_REFT_HEX_X, 21, order);
		error = sqrt(error);
		CHECK_ERROR;

		//
		error = 0.0;
		error += proj->get_error(H3D_REFT_HEX_Y, 22, order);
		error += proj->get_error(H3D_REFT_HEX_Y, 23, order);
		error = sqrt(error);
		CHECK_ERROR;

		//
		error = 0.0;
		error += proj->get_error(H3D_REFT_HEX_Z, 24, order);
		error += proj->get_error(H3D_REFT_HEX_Z, 25, order);
		error = sqrt(error);
		CHECK_ERROR;

		//
		error = 0.0;
		error += proj->get_error(H3D_H3D_REFT_HEX_XY,  8, order);
		error += proj->get_error(H3D_H3D_REFT_HEX_XY,  9, order);
		error += proj->get_error(H3D_H3D_REFT_HEX_XY, 10, order);
		error += proj->get_error(H3D_H3D_REFT_HEX_XY, 11, order);
		error = sqrt(error);
		CHECK_ERROR;

		//
		error = 0.0;
		error += proj->get_error(H3D_H3D_REFT_HEX_XZ, 12, order);
		error += proj->get_error(H3D_H3D_REFT_HEX_XZ, 13, order);
		error += proj->get_error(H3D_H3D_REFT_HEX_XZ, 14, order);
		error += proj->get_error(H3D_H3D_REFT_HEX_XZ, 15, order);
		error = sqrt(error);
		CHECK_ERROR;

		//
		error = 0.0;
		error += proj->get_error(H3D_H3D_REFT_HEX_YZ, 16, order);
		error += proj->get_error(H3D_H3D_REFT_HEX_YZ, 17, order);
		error += proj->get_error(H3D_H3D_REFT_HEX_YZ, 18, order);
		error += proj->get_error(H3D_H3D_REFT_HEX_YZ, 19, order);
		error = sqrt(error);
		CHECK_ERROR;

		//
		error = 0.0;
		for (int j = 0; j < 8; j++)
			error += proj->get_error(H3D_H3D_H3D_REFT_HEX_XYZ, j, order);
		error = sqrt(error);
		CHECK_ERROR;

		delete proj;
	}

#ifdef WITH_PETSC
	PetscFinalize();
#endif

	return ret;
}
