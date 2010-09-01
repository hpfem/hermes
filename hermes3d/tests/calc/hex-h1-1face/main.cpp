// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
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

/*
 * hang-nodes-continuity.cc
 *
 * usage: $0 <mesh file> <element id> <refinement id> [<element id> <refinement id>...]
 *
 */

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

#define BEGIN_BLOCK							{
#define END_BLOCK							}

//#define DIRICHLET
#define NEWTON

//#define X2_Y2_Z2
//#define XM_YN_ZO
#define XM_YN_ZO_2


// Everything is tested on the following geometry
//
//
//      7             6             11
//        +-----------+-----------+
//       /|          /|          /|
//      / |         / |         / |
//     /  |      5 /  |     10 /  |
//  4 +-----------+-----------+   |
//    |   |       |   |       |   |
//    |   +-------|---+-------|---+
//    |  / 3      |  / 2      |  / 9
//    | /         | /         | /
//    |/          |/          |/
//    +-----------+-----------+
//   0            1            8
//
//  ^ z
//  |
//  | / y
//  |/
//  +---> x

// vertices
Point3D vtcs[] = {
	{ -2, -1, -1 },
	{  0, -1, -1 },
	{  0,  1, -1 },
	{ -2,  1, -1 },
	{ -2, -1,  1 },
	{  0, -1,  1 },
	{  0,  1,  1 },
	{ -2,  1,  1 },
	{  2, -1, -1 },
	{  2,  1, -1 },
	{  2, -1,  1 },
	{  2,  1,  1 },
};

// mesh
Word_t hexs[2][48][8] = {
	// hexs 1
	{
		{  0,  1,  2,  3,  4,  5,  6,  7 },
		{  1,  2,  3,  0,  5,  6,  7,  4 },
		{  2,  3,  0,  1,  6,  7,  4,  5 },
		{  3,  0,  1,  2,  7,  4,  5,  6 },

		{  4,  7,  6,  5,  0,  3,  2,  1 },
		{  5,  4,  7,  6,  1,  0,  3,  2 },
		{  6,  5,  4,  7,  2,  1,  0,  3 },
		{  7,  6,  5,  4,  3,  2,  1,  0 },

		{  1,  5,  6,  2,  0,  4,  7,  3 },
		{  5,  6,  2,  1,  4,  7,  3,  0 },
		{  6,  2,  1,  5,  7,  3,  0,  4 },
		{  2,  1,  5,  6,  3,  0,  4,  7 },

		{  0,  3,  7,  4,  1,  2,  6,  5 },
		{  4,  0,  3,  7,  5,  1,  2,  6 },
		{  7,  4,  0,  3,  6,  5,  1,  2 },
		{  3,  7,  4,  0,  2,  6,  5,  1 },

		{  1,  0,  4,  5,  2,  3,  7,  6 },
		{  0,  4,  5,  1,  3,  7,  6,  2 },
		{  4,  5,  1,  0,  7,  6,  2,  3 },
		{  5,  1,  0,  4,  6,  2,  3,  7 },

		{  2,  6,  7,  3,  1,  5,  4,  0 },
		{  3,  2,  6,  7,  0,  1,  5,  4 },
		{  7,  3,  2,  6,  4,  0,  1,  5 },
		{  6,  7,  3,  2,  5,  4,  0,  1 },

		{  0,  3,  7,  4,  1,  2,  6,  5 },
		{  4,  0,  3,  7,  5,  1,  2,  6 },
		{  7,  4,  0,  3,  6,  5,  1,  2 },
		{  3,  7,  4,  0,  2,  6,  5,  1 },

		{  1,  5,  6,  2,  0,  4,  7,  3 },
		{  5,  6,  2,  1,  4,  7,  3,  0 },
		{  6,  2,  1,  5,  7,  3,  0,  4 },
		{  2,  1,  5,  6,  3,  0,  4,  7 },

		{  7,  6,  5,  4,  3,  2,  1,  0 },
		{  4,  7,  6,  5,  0,  3,  2,  1 },
		{  5,  4,  7,  6,  1,  0,  3,  2 },
		{  6,  5,  4,  7,  2,  1,  0,  3 },

		{  3,  0,  1,  2,  7,  4,  5,  6 },
		{  0,  1,  2,  3,  4,  5,  6,  7 },
		{  1,  2,  3,  0,  5,  6,  7,  4 },
		{  2,  3,  0,  1,  6,  7,  4,  5 },

		{  2,  6,  7,  3,  1,  5,  4,  0 },
		{  6,  7,  3,  2,  5,  4,  0,  1 },
		{  7,  3,  2,  6,  4,  0,  1,  5 },
		{  3,  2,  6,  7,  0,  1,  5,  4 },

		{  1,  0,  4,  5,  2,  3,  7,  6 },
		{  5,  1,  0,  4,  6,  2,  3,  7 },
		{  4,  5,  1,  0,  7,  6,  2,  3 },
		{  0,  4,  5,  1,  3,  7,  6,  2 }
	},
	// hexs 2
	{
		{  1,  8,  9,  2,  5, 10, 11,  6 },
		{  8,  9,  2,  1, 10, 11,  6,  5 },
		{  9,  2,  1,  8, 11,  6,  5, 10 },
		{  2,  1,  8,  9,  6,  5, 10, 11 },

		{  5,  6, 11, 10,  1,  2,  9,  8 },
		{ 10,  5,  6, 11,  8,  1,  2,  9 },
		{ 11, 10,  5,  6,  9,  8,  1,  2 },
		{  6, 11, 10,  5,  2,  9,  8,  1 },

		{  8, 10, 11,  9,  1,  5,  6,  2 },
		{ 10, 11,  9,  8,  5,  6,  2,  1 },
		{ 11,  9,  8, 10,  6,  2,  1,  5 },
		{  9,  8, 10, 11,  2,  1,  5,  6 },

		{  1,  2,  6,  5,  8,  9, 11, 10 },
		{  5,  1,  2,  6, 10,  8,  9, 11 },
		{  6,  5,  1,  2, 11, 10,  8,  9 },
		{  2,  6,  5,  1,  9, 11, 10,  8 },

		{  8,  1,  5, 10,  9,  2,  6, 11 },
		{  1,  5, 10,  8,  2,  6, 11,  9 },
		{  5, 10,  8,  1,  6, 11,  9,  2 },
		{ 10,  8,  1,  5, 11,  9,  2,  6 },

		{  9, 11,  6,  2,  8, 10,  5,  1 },
		{  2,  9, 11,  6,  1,  8, 10,  5 },
		{  6,  2,  9, 11,  5,  1,  8, 10 },
		{ 11,  6,  2,  9, 10,  5,  1,  8 },

		{  1,  2,  6,  5,  8,  9, 11, 10 },
		{  5,  1,  2,  6, 10,  8,  9, 11 },
		{  6,  5,  1,  2, 11, 10,  8,  9 },
		{  2,  6,  5,  1,  9, 11, 10,  8 },

		{  8, 10, 11,  9,  1,  5,  6,  2 },
		{ 10, 11,  9,  8,  5,  6,  2,  1 },
		{ 11,  9,  8, 10,  6,  2,  1,  5 },
		{  9,  8, 10, 11,  2,  1,  5,  6 },

		{  6, 11, 10,  5,  2,  9,  8,  1 },
		{  5,  6, 11, 10,  1,  2,  9,  8 },
		{ 10,  5,  6, 11,  8,  1,  2,  9 },
		{ 11, 10,  5,  6,  9,  8,  1,  2 },

		{  2,  1,  8,  9,  6,  5, 10, 11 },
		{  1,  8,  9,  2,  5, 10, 11,  6 },
		{  8,  9,  2,  1, 10, 11,  6,  5 },
		{  9,  2,  1,  8, 11,  6,  5,  10 },

		{  9, 11,  6,  2,  8, 10,  5,  1 },
		{ 11,  6,  2,  9, 10,  5,  1,  8 },
		{  6,  2,  9, 11,  5,  1,  8, 10 },
		{  2,  9, 11,  6,  1,  8, 10,  5 },

		{  8,  1,  5, 10,  9,  2,  6, 11 },
		{ 10,  8,  1,  5, 11,  9,  2,  6 },
		{  5, 10,  8,  1,  6, 11,  9,  2 },
		{  1,  5, 10,  8,  2,  6, 11,  9 }
	}
};

Word_t bnd[10][5] = {
	{ 0,  3,  7,  4, 1 },
	{ 8,  9, 11, 10, 2 },
	{ 0,  1,  5,  4, 3 },
	{ 1,  8, 10,  5, 3 },
	{ 3,  2,  6,  7, 4 },
	{ 2,  9, 11,  6, 4 },
	{ 0,  1,  2,  3, 5 },
	{ 1,  8,  9,  2, 5 },
	{ 4,  5,  6,  7, 6 },
	{ 5, 10, 11,  6, 6 }
};


int m = 2, n = 2, o = 2;

template<typename T>
T fnc(T x, T y, T z) {
#ifdef XM_YN_ZO
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
#elif defined XM_YN_ZO_2
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 2) * z + pow(z, 4);
#elif defined X2_Y2_Z2
	return x*x + y*y + z*z;
#endif
}

template<typename T>
T dfnc(T x, T y, T z) {
#ifdef XM_YN_ZO
	T ddxx = m * (m - 1) * pow(x, m - 2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
	T ddyy = n * (n - 1) * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = o * (o - 1) * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);
	return -(ddxx + ddyy + ddzz);

#elif defined XM_YN_ZO_2
	T ddxx = m*(m-1) * pow(x, m-2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 2 * z;
	T ddyy = n*(n-1) * pow(x, m) * pow(y, n-2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = o*(o-1) * pow(x, m) * pow(y, n) * pow(z, o-2) + 12 * pow(z, 2);
	return -(ddxx + ddyy + ddzz);
#elif defined X2_Y2_Z2
	return -6.0;
#endif
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
#ifdef XM_YN_ZO
	dx = m * pow(x, m - 1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n - 1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o - 1) - pow(x, 3) + 4 * pow(z, 3);
#elif defined XM_YN_ZO_2
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 2 * x * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 2) + 4 * pow(z, 3);
#elif defined X2_Y2_Z2
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;
#endif

	return fnc(x, y, z);
}

//

BCType bc_types(int marker) {
#ifdef DIRICHLET
	return BC_ESSENTIAL;
#elif defined NEWTON
	return BC_NATURAL;
#endif
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z) {
#ifdef DIRICHLET
	return fnc(x, y, z);
#else
	return 0;
#endif
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t bilinear_form_surf(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_u_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return int_F_v<f_t, res_t>(n, wt, dfnc, u, e);
}

template<typename f_t, typename res_t>
res_t linear_form_surf(int np, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data) {
	res_t result = 0;
	for (int i = 0; i < np; i++) {
#ifdef XM_YN_ZO
		res_t dx = m * pow(e->x[i], m - 1) * pow(e->y[i], n) * pow(e->z[i], o) + 2 * e->x[i] * pow(e->y[i], 3) - 3 * pow(e->x[i], 2) * e->z[i];
		res_t dy = n * pow(e->x[i], m) * pow(e->y[i], n - 1) * pow(e->z[i], o) + 3 * pow(e->x[i], 2) * pow(e->y[i], 2);
		res_t dz = o * pow(e->x[i], m) * pow(e->y[i], n) * pow(e->z[i], o - 1) - pow(e->x[i], 3) + 4 * pow(e->z[i], 3);
#elif defined XM_YN_ZO_2
		res_t dx = m * pow(e->x[i], m-1) * pow(e->y[i], n) * pow(e->z[i], o) + 2 * e->x[i] * pow(e->y[i], 3) - 2 * e->x[i] * e->z[i];
		res_t dy = n * pow(e->x[i], m) * pow(e->y[i], n-1) * pow(e->z[i], o) + 3 * pow(e->x[i], 2) * pow(e->y[i], 2);
		res_t dz = o * pow(e->x[i], m) * pow(e->y[i], n) * pow(e->z[i], o-1) - pow(e->x[i], 2) + 4 * pow(e->z[i], 3);
#elif defined X2_Y2_Z2
		res_t dx = 2 * e->x[i];
		res_t dy = 2 * e->y[i];
		res_t dz = 2 * e->z[i];
#endif
		result += wt[i] * (u->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc(e->x[i], e->y[i], e->z[i])));
	}
	return result;
}

// helpers ////////////////////////////////////////////////////////////////////////////////////////

int parse_reft(char *str) {
	if (strcasecmp(str, "x") == 0) return H3D_REFT_HEX_X;
	else if (strcasecmp(str, "y") == 0) return H3D_REFT_HEX_Y;
	else if (strcasecmp(str, "z") == 0) return H3D_REFT_HEX_Z;
	else if (strcasecmp(str, "xy") == 0 || strcasecmp(str, "yx") == 0) return H3D_H3D_REFT_HEX_XY;
	else if (strcasecmp(str, "xz") == 0 || strcasecmp(str, "zx") == 0) return H3D_H3D_REFT_HEX_XZ;
	else if (strcasecmp(str, "yz") == 0 || strcasecmp(str, "zy") == 0) return H3D_H3D_REFT_HEX_YZ;
	else if (strcasecmp(str, "xyz") == 0) return H3D_H3D_H3D_REFT_HEX_XYZ;
	else return H3D_REFT_HEX_NONE;
}

///

const double EPS = 10e-14;


// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	int res = ERR_SUCCESS;


#ifdef WITH_PETSC
	PetscInitialize(NULL, NULL, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

	try {
		for (int i = 0; i < 48; i++) {
			for (int j = 0; j < 48; j++) {
//		int i = 5; {
//		int j = 0; {
				printf("Config: %d, %d ", i, j);

				Mesh mesh;

				for (Word_t k = 0; k < countof(vtcs); k++)
					mesh.add_vertex(vtcs[k].x, vtcs[k].y, vtcs[k].z);
				Word_t h1[] = {
						hexs[0][i][0] + 1, hexs[0][i][1] + 1, hexs[0][i][2] + 1, hexs[0][i][3] + 1,
						hexs[0][i][4] + 1, hexs[0][i][5] + 1, hexs[0][i][6] + 1, hexs[0][i][7] + 1 };
				mesh.add_hex(h1);
				Word_t h2[] = {
						hexs[1][j][0] + 1, hexs[1][j][1] + 1, hexs[1][j][2] + 1, hexs[1][j][3] + 1,
						hexs[1][j][4] + 1, hexs[1][j][5] + 1, hexs[1][j][6] + 1, hexs[1][j][7] + 1 };
				mesh.add_hex(h2);
				// bc
				for (Word_t k = 0; k < countof(bnd); k++) {
					Word_t facet_idxs[Quad::NUM_VERTICES] = { bnd[k][0] + 1, bnd[k][1] + 1, bnd[k][2] + 1, bnd[k][3] + 1 };
					mesh.add_quad_boundary(facet_idxs, bnd[k][4]);
				}

				mesh.ugh();
//				mesh.dump();

//				Element *hx[] = { mesh.elements[1], mesh.elements[2] };
//				printf("[%d, %d]\n", hx[0]->get_face_orientation(1), hx[1]->get_face_orientation(2));

//				Word_t fidx[4];
//				hx[1]->get_face_vertices(2, fidx);
//				printf("FI: %d, %d, %d, %d\n", fidx[0], fidx[1], fidx[2], fidx[3]);
				printf("\n");

#ifdef OUTPUT_DIR
				BEGIN_BLOCK
					// output the mesh
					const char *of_name = OUTPUT_DIR "/ref.msh";
					FILE *ofile = fopen(of_name, "w");
					if (ofile != NULL) {
						GmshOutputEngine output(ofile);
						output.out(&mesh);
						fclose(ofile);
					}
					else {
						warning("Can not open '%s' for writing.", of_name);
					}
				END_BLOCK
#endif

				H1ShapesetLobattoHex shapeset;

//				printf("* Setting the space up\n");
				H1Space space(&mesh, &shapeset);
				space.set_bc_types(bc_types);
				space.set_essential_bc_values(essential_bc_values);

#ifdef XM_YN_ZO
				order3_t ord(4, 4, 4);
#elif defined XM_YN_ZO_2
				order3_t ord(4, 4, 4);
#elif defined X2_Y2_Z2
				order3_t ord(2, 2, 2);
#endif
//				printf("  - Setting uniform order to (%d, %d, %d)\n", dir_x, dir_y, dir_z);
				space.set_uniform_order(ord);

				space.assign_dofs();

//				printf("* Calculating a solution\n");

#if defined WITH_UMFPACK
				UMFPackMatrix mat;
				UMFPackVector rhs;
				UMFPackLinearSolver solver(&mat, &rhs);
#elif defined WITH_PARDISO
				PardisoLinearSolver solver;
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
#ifdef DIRICHLET
				wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM);
				wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>);
#elif defined NEWTON
				wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM);
				wf.add_matrix_form_surf(bilinear_form_surf<double, scalar>, bilinear_form_surf<ord_t, ord_t>);
				wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>);
				wf.add_vector_form_surf(linear_form_surf<double, scalar>, linear_form_surf<ord_t, ord_t>);
#endif

				LinProblem lp(&wf);
				lp.set_space(&space);

				// assemble stiffness matrix
				lp.assemble(&mat, &rhs);

				// solve the stiffness matrix
				bool solved = solver.solve();
				if (!solved) throw ERR_FAILURE;

//				{
//					char file_name[1024];
//					sprintf(file_name, "%s/matrix-%d-%d", OUTPUT_DIR, i, j);
//					FILE *file = fopen(file_name, "w");
//					if (file != NULL) {
//						solver.dump_matrix(file, "A");
//						solver.dump_rhs(file, "b");
//
//						fclose(file);
//					}
//				}

				Solution sln(&mesh);
				sln.set_fe_solution(&space, solver.get_solution());

				ExactSolution exsln(&mesh, exact_solution);
				// norm
				double h1_sln_norm = h1_norm(&sln);
				double h1_err_norm = h1_error(&sln, &exsln);
				printf(" - H1 solution norm:   % le\n", h1_sln_norm);
				printf(" - H1 error norm:      % le\n", h1_err_norm);

				double l2_sln_norm = l2_norm(&sln);
				double l2_err_norm = l2_error(&sln, &exsln);
				printf(" - L2 solution norm:   % le\n", l2_sln_norm);
				printf(" - L2 error norm:      % le\n", l2_err_norm);

				assert(h1_sln_norm > 0 && h1_err_norm > 0);
				assert(l2_sln_norm > 0 && l2_err_norm > 0);

//				// out fn
//				char fname[4096];
//				sprintf(fname, "%s/cfg-%d-%d.pos", OUTPUT_DIR, i, j);
//				FILE *fnf = fopen(fname, "w");
//				assert(fnf != NULL);
//				GmshOutputEngine out(fnf);
//				char var[64];
//				sprintf(var, "%d_%d", i, j);
//				out.out(&sln, var);
//				fclose(fnf);
//
//				char mfname[4096];
//				sprintf(mfname, "%s/mesh-%d-%d.ref", OUTPUT_DIR, i, j);
//				FILE *mfnf = fopen(mfname, "w");
//				assert(mfnf != NULL);
//				GmshOutputEngine outm(mfnf);
//				outm.out(&mesh);
//				fclose(mfnf);

				if (h1_err_norm > EPS || l2_err_norm > EPS) {
					// calculated solution is not enough precise
					printf("Solution is not precise enough.\n");
					throw ERR_FAILURE;
				}

				printf("Passed\n");
			}
		}
	}
	catch (int e) {
		res = e;
		printf("Failed\n");
	}

#ifdef WITH_PETSC
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

