// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - Pavel Kus
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
// hang-nodes-continuity.cc
//
// usage: $0 <mesh file> <element id> <refinement id> [<element id> <refinement id>...]
//
//

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

#define BEGIN_BLOCK						{
#define END_BLOCK						}


//#define X2_Y2_Z2
#define FN4

// general polynomial function satisfying perfect conductor bc
// exact solution has zero tangential component on the boundary on the domain (-1, 1)^3
// DO NOT TEST on other domains

const double alpha = 1.0;

scalar3 &exact_solution(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz) {
	_F_
	static scalar3 val;

#ifdef FN4
	dx[0] = (1 - y*y)*(1 - z*z)*(z - 6*x*x);
	dx[1] = 2*(1 - x*x)*(1 - z*z) - 2*x*(1 - z*z)*(y*y*y + 2*x);
	dx[2] = -2*x*(1 - y*y)*(z*z - 3*x*y*z) - 3*y*z*(1 - x*x)*(1 - y*y);

	dy[0] = -2*y*(1 - z*z)*(1 - 2*x*x*x + x*z);
	dy[1] = 3*y*y*(1 - x*x)*(1 - z*z);
	dy[2] = -2*y*(1 - x*x)*(z*z - 3*x*y*z) - 3*x*z*(1 - x*x)*(1 - y*y);

	dz[0] = x*(1 - y*y)*(1 - z*z) - 2*z*(1 - y*y)*(1 - 2*x*x*x + x*z);
	dz[1] = -2*z*(1 - x*x)*(y*y*y + 2*x);
	dz[2] = (1 - x*x)*(1 - y*y)*(2*z - 3*x*y);

	val[0] = (1-y*y) * (1-z*z) * (x*z - 2*x*x*x + 1);
	val[1] = (1-x*x) * (1-z*z) * (y*y*y + 2*x);
	val[2] = (1-x*x) * (1-y*y) * (z*z - 3*x*y*z);
#elif defined X2_Y2_Z2
	dx[0] = 0;
	dx[1] = 0;
	dx[2] = -2 * x * (1 - y*y) * (1 - z*z);

	dy[0] = 0;
	dy[1] = 0;
	dy[2] = -2 * (1 - x*x) * y * (1 - z*z);

	dz[0] = 0;
	dz[1] = 0;
	dz[2] = -2 * (1 - x*x) * (1 - y*y) * z;

	val[0] = 0;
	val[1] = 0;
	val[2] = (1 - x*x) * (1 - y*y) * (1 - z*z);
#endif

	return val;
}

template<typename S, typename T>
void exact_sln(S x, S y, S z, T (&val)[3], T (&dx)[3], T (&dy)[3], T (&dz)[3]) {
	_F_
#ifdef FN4
	val[0] = (1-y*y) * (1-z*z) * (x*z - 2*x*x*x + 1);
	val[1] = (1-x*x) * (1-z*z) * (y*y*y + 2*x);
	val[2] = (1-x*x) * (1-y*y) * (z*z - 3*x*y*z);

	dx[0] = (1 - y*y)*(1 - z*z)*(z - 6*x*x);
	dx[1] = 2*(1 - x*x)*(1 - z*z) - 2*x*(1 - z*z)*(y*y*y + 2*x);
	dx[2] = -2*x*(1 - y*y)*(z*z - 3*x*y*z) - 3*y*z*(1 - x*x)*(1 - y*y);

	dy[0] = -2*y*(1 - z*z)*(1 - 2*x*x*x + x*z);
	dy[1] = 3*y*y*(1 - x*x)*(1 - z*z);
	dy[2] = -2*y*(1 - x*x)*(z*z - 3*x*y*z) - 3*x*z*(1 - x*x)*(1 - y*y);

	dz[0] = x*(1 - y*y)*(1 - z*z) - 2*z*(1 - y*y)*(1 - 2*x*x*x + x*z);
	dz[1] = -2*z*(1 - x*x)*(y*y*y + 2*x);
	dz[2] = (1 - x*x)*(1 - y*y)*(2*z - 3*x*y);
#elif defined X2_Y2_Z2
	dx[0] = 0;
	dx[1] = 0;
	dx[2] = -2 * x * (1 - y*y) * (1 - z*z);

	dy[0] = 0;
	dy[1] = 0;
	dy[2] = -2 * (1 - x*x) * y * (1 - z*z);

	dz[0] = 0;
	dz[1] = 0;
	dz[2] = -2 * (1 - x*x) * (1 - y*y) * z;

	val[0] = 0;
	val[1] = 0;
	val[2] = (1 - x*x) * (1 - y*y) * (1 - z*z);
#endif
}

template<typename S, typename T>
void f(S x, S y, S z, T (&val)[3]) {
	_F_
	T ev[3], dx[3], dy[3], dz[3];
	exact_sln(x, y, z, ev, dx, dy, dz);

#ifdef FN4
	T curlpart[3] = {
		2*(1 - y*y)*(1 - 2*x*x*x + x*z) + 2*(1 - z*z)*(1 - 2*x*x*x + x*z) - 6*x*y*y*(1 - z*z) - 3*y*(1 - x*x)*(1 - y*y) - 2*x*(1 - y*y)*(2*z - 3*x*y) + 4*x*z*(1 - y*y),
		2*(1 - x*x)*(y*y*y + 2*x) + 2*(1 - z*z)*(y*y*y + 2*x) + 8*x*(1 - z*z) - 3*x*(1 - x*x)*(1 - y*y) - 2*y*(1 - x*x)*(2*z - 3*x*y) - 2*y*(1 - z*z)*(z - 6*x*x),
		(1 - y*y)*(1 - z*z) + 2*(1 - x*x)*(z*z - 3*x*y*z) + 2*(1 - y*y)*(z*z - 3*x*y*z) - 6*z*y*y*(1 - x*x) - 2*z*(1 - y*y)*(z - 6*x*x) - 12*x*y*z*(1 - x*x) - 12*x*y*z*(1 - y*y)
	};
#elif defined X2_Y2_Z2
	// \nabla x \nabla x val
	T curlpart[3] = {
		4 * x * (1 - y*y) * z,
		4 * (1 - x*x) * y * z,
		2 * (1 - y*y) * (1 - z*z) + 2 * (1 - x*x) * (1 - z*z)
	};
#endif

	val[0] = curlpart[0] - alpha * ev[0];
	val[1] = curlpart[1] - alpha * ev[1];
	val[2] = curlpart[2] - alpha * ev[2];
}

BCType bc_types(int marker) {
	return BC_ESSENTIAL;
}

// definition of the forms

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return
		hcurl_int_curl_u_curl_v<f_t, res_t>(n, wt, u, v, e) -
		alpha * hcurl_int_u_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data) {
	return hcurl_int_F_v<f_t, res_t>(n, wt, f<f_t, res_t>, u, e);
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


////////////////////////////////////////////////////////////////////////////////

// maximal level of refinement considered (= 0 .. # of ref)
#ifdef DEV_TESTS
	#define MAX_LEVEL					4
#else
	#define MAX_LEVEL					2
#endif
#define NUM_RULES						((MAX_LEVEL + 1) * (MAX_LEVEL + 1))

// special quadrature used in this test to define set of points, where
// continuity is tested. Quadrature is used in order to use RefMap abilities
// calculate physical coordinates of points from the reference domain.
//
// points are chosen in this way:
// - only face points are used
// - there are more levels of points
// - 0-level consist of 16 points chosen symmetricaly
// - other levels (or orders) simulate division of adjacent elements and map
//   points in such way, that some points from divided and some from nondivided
//   face will still match in the physical domain
// - this is of course only possible up to some level of division, controled by
//   constant LEVEL. for example, for hex1 :
//   LEVEL = 1  divisions 0 x 1 y or 0 x 1 y 3 z are ok, but 0 x 1 y 3 y not
//   LEVEL = 2  division 0 x 1 y 3 y is ok, but 0 x 1 y 3 y 5 y not
// - if LEVEL is not sufficient, there will be some faces, that will not be tested,
//   because no points from the face will match to points from the constraining face
class ContQuad : public Quad3D {
public:
	ContQuad() {
		_F_
//		max_order = max_edge_order = max_face_order = NUM_RULES;

		int my_np_1d[MAX_LEVEL + 1];
		double my_tables_1d[MAX_LEVEL + 1][1000];

		my_np_1d[0] = 4;
		my_tables_1d[0][0] = -0.71;  my_tables_1d[0][1] = -0.59;
		my_tables_1d[0][2] =  0.59;  my_tables_1d[0][3] =  0.71;

		for(int order = 1; order <= MAX_LEVEL; order++) {
			my_np_1d[order] = 2 * my_np_1d[order - 1];
			for(int i = 0; i < my_np_1d[order - 1]; i++) {
				my_tables_1d[order][i] = (my_tables_1d[order - 1][i] - 1.) / 2.;
				my_tables_1d[order][i + my_np_1d[order - 1]] = (my_tables_1d[order - 1][i] + 1.) / 2.;
			}
		}

		face_tables = new Array<QuadPt3D *>[Hex::NUM_FACES];
		for (int order = 0; order < NUM_RULES; order++) {
			int ord1 = order / (MAX_LEVEL + 1);
			int ord2 = order % (MAX_LEVEL + 1);
			np_face[order] = my_np_1d[ord1] * my_np_1d[ord2];
			for (int face = 0; face < Hex::NUM_FACES; face++)
				face_tables[face][order] = new QuadPt3D[np_face[order]];

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[0][order][n].x = -1;
					face_tables[0][order][n].y = my_tables_1d[ord1][k];
					face_tables[0][order][n].z = my_tables_1d[ord2][l];
				}
			}

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[1][order][n].x = 1;
					face_tables[1][order][n].y = my_tables_1d[ord1][k];
					face_tables[1][order][n].z = my_tables_1d[ord2][l];
				}
			}

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[2][order][n].x = my_tables_1d[ord1][k];
					face_tables[2][order][n].y = -1;
					face_tables[2][order][n].z = my_tables_1d[ord2][l];
				}
			}

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[3][order][n].x = my_tables_1d[ord1][k];
					face_tables[3][order][n].y = 1;
					face_tables[3][order][n].z = my_tables_1d[ord2][l];
				}
			}

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[4][order][n].x = my_tables_1d[ord1][k];
					face_tables[4][order][n].y = my_tables_1d[ord2][l];
					face_tables[4][order][n].z = -1;
				}
			}

			for (int k = 0, n = 0; k < my_np_1d[ord1]; k++) {
				for (int l = 0; l < my_np_1d[ord2]; l++, n++) {
					assert(n < np_face[order]);
					face_tables[5][order][n].x = my_tables_1d[ord1][k];
					face_tables[5][order][n].y = my_tables_1d[ord2][l];
					face_tables[5][order][n].z = 1;
				}
			}
		}
	}

	~ContQuad() {
		_F_
		for (int face = 0; face < Hex::NUM_FACES; face++)
			for (int order = 0; order < NUM_RULES; order++)
				delete [] face_tables[face][order];
		delete [] face_tables;
	}
};

struct Point {
	double ref_x, ref_y, ref_z;
	double phys_x, phys_y, phys_z;
	Word_t elm_idx;
	int iface;

	Point(int idx, int ifa, double rx, double ry, double rz, double px, double py, double pz) {
		elm_idx = idx; iface = ifa;
		ref_x = rx;  ref_y = ry;  ref_z = rz;
		phys_x = px; phys_y = py; phys_z = pz;
	}
};

typedef
	int (*compfn)(const void*, const void*);

int compare(Point **pt1, Point **pt2) {
	_F_
	double val1 = 1000000. * (*pt1)->phys_x + 1000. * (*pt1)->phys_y + (*pt1)->phys_z;
	double val2 = 1000000. * (*pt2)->phys_x + 1000. * (*pt2)->phys_y + (*pt2)->phys_z;

	if (val1 < val2) return -1;
	else if (val1 > val2) return 1;
	else return 0;
}

const double EPS = 1e-13;
const double TOLERANCE = 1e-10;

bool equal(Point *pt1, Point *pt2) {
	_F_
	if (fabs(pt1->phys_x - pt2->phys_x) > TOLERANCE) return false;
	if (fabs(pt1->phys_y - pt2->phys_y) > TOLERANCE) return false;
	if (fabs(pt1->phys_z - pt2->phys_z) > TOLERANCE) return false;
	return true;
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args) {
	_F_
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(NULL, NULL, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	if (argc < 2) error("Not enough parameters");

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh)) error("Loading mesh file '%s'\n", args[1]);

	// apply refinements
	for (int i = 2; i < argc; i += 2) {
		int elem_id, reft_id;
		sscanf(args[i], "%d", &elem_id);
		reft_id = parse_reft(args[i + 1]);
		mesh.refine_element(elem_id + 1, reft_id);
	}
//	mesh.dump();


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

	HcurlShapesetLobattoHex shapeset;

	printf("* Setting the space up\n");
	HcurlSpace space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
//	space.set_essential_bc_values(essential_bc_values);

#if defined FN4
	order3_t o(4, 4, 4);
#elif defined X2_Y2_Z2
	order3_t o(2, 2, 2);
#endif
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

	WeakForm wf(1);
	wf.add_matrix_form(0, 0, bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM);
	wf.add_vector_form(0, linear_form<double, scalar>, linear_form<ord_t, ord_t>);

	// assemble stiffness matrix
	LinearProblem lp(&wf);
	lp.set_space(&space);

	lp.assemble(&mat, &rhs);

#if 0 //def OUTPUT_DIR
	{
		char file_name[1024];
		sprintf(file_name, "%s/matrix", OUTPUT_DIR);
		FILE *file = fopen(file_name, "w");
		if (file != NULL) {
			mat.dump(file, "A");
			rhs.dump(file, "b");

			fclose(file);
		}
	}
#endif

	try {
		Solution sln(&mesh);
#if 1
		// solve the stiffness matrix
		bool solved = solver.solve();
		if (!solved) throw ERR_FAILURE;

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

		if (h1_err_norm > EPS || l2_err_norm > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
			printf("Solution is not precise enough.\n");
		}
#endif
		//
		//
		// the main code starts here
		//
		//

#if 0
		RefMap rm(&mesh);
		{
			// prepare solution correspondig to basis function with dof dof
			double sln_vector[ndofs + 1];
			memset(sln_vector, 0, (ndofs + 1) * sizeof(double));
			sln_vector[0] = 1.0;
			sln.set_fe_solution(&space, sln_vector);

			rm.set_active_element(mesh.elements[2]);

			int np = 1;
			QuadPt3D pt[1];
			pt[0].x =  1.0;
			pt[0].y = -0.5;
			pt[0].z =  0.0;
			double *phys_x = rm.get_phys_x(np, pt);
			double *phys_y = rm.get_phys_y(np, pt);
			double *phys_z = rm.get_phys_z(np, pt);

			printf("------------------------------------------------\n");
			printf("% lf, % lf, % lf\n", phys_x[0], phys_y[0], phys_z[0]);
			sln.set_active_element(mesh.elements[2]);
			sln.precalculate(1, pt, FN_VAL);
			double *v0 = sln.get_fn_values(0);
			double *v1 = sln.get_fn_values(1);
			double *v2 = sln.get_fn_values(2);
			printf("val = % lf, % lf, % lf\n", v0[0], v1[0], v2[0]);

			///
			{
				rm.set_active_element(mesh.elements[4]);


				pt[0].x = -1.0;
				pt[0].y =  0.0;
				pt[0].z =  0.0;
				double *phys_x = rm.get_phys_x(np, pt);
				double *phys_y = rm.get_phys_y(np, pt);
				double *phys_z = rm.get_phys_z(np, pt);

				printf("------------------------------------------------\n");
				printf("% lf, % lf, % lf\n", phys_x[0], phys_y[0], phys_z[0]);
				sln.set_active_element(mesh.elements[4]);
				sln.precalculate(1, pt, FN_VAL);
				double *v0 = sln.get_fn_values(0);
				double *v1 = sln.get_fn_values(1);
				double *v2 = sln.get_fn_values(2);
				printf("val = % lf, % lf, % lf\n", v0[0], v1[0], v2[0]);

			}
		}
#endif

#if 1
		ContQuad my_quad;
		RefMap ref_map(&mesh);

		int num_points = 0;
		for (int order = 0; order < NUM_RULES; order++) {
			FOR_ALL_ACTIVE_ELEMENTS(idx, &mesh) {
				for (int iface = 0; iface < Hex::NUM_FACES; iface++)
					num_points += my_quad.get_face_num_points(iface, order);
			}
		}

		Point **points = new Point *[num_points];
		int ipt = 0;

		// find points
		for (int order = 0; order < NUM_RULES; order++) {
			FOR_ALL_ACTIVE_ELEMENTS(idx, &mesh) {
				Element *e = mesh.elements[idx];
				ref_map.set_active_element(e);
//				ref_map.set_quad(&my_quad);
				for (int iface = 0; iface < Hex::NUM_FACES; iface++) {
//					Word_t fac_idx = mesh.get_facet_id(e, iface);

					QuadPt3D *quad_pts = my_quad.get_face_points(iface, order);
					int np = my_quad.get_face_num_points(iface, order);
//					double *phys_x = ref_map.get_face_phys_x(iface, order);
//					double *phys_y = ref_map.get_face_phys_y(iface, order);
//					double *phys_z = ref_map.get_face_phys_z(iface, order);
					double *phys_x = ref_map.get_phys_x(np, quad_pts);
					double *phys_y = ref_map.get_phys_y(np, quad_pts);
					double *phys_z = ref_map.get_phys_z(np, quad_pts);

					// for each face and each integration point store reference and physical coordinates
					for (int pt = 0; pt < np; pt++) {
						points[ipt++] = new Point(idx, iface,
							quad_pts[pt].x, quad_pts[pt].y, quad_pts[pt].z,
							phys_x[pt], phys_y[pt], phys_z[pt]);
					}
				}
			}
		}

		// sort points according to first phys_x, then phys_y and phys_z
		// it means, that two points, with almost identical physical coordinates
		// (even though from different elements) will be next to each other in the array
		qsort((void *) points, num_points, sizeof(Point *), (compfn) compare);

		int *pairs = new int [num_points];
		int num_pairs = 0;

		// choose those indicies, that correspond to pairs with identical physical coordinates
		// and store them in field pairs
		for (int i = 0; i < num_points - 1; i++) {
			if (equal(points[i], points[i+1])) {
				pairs[num_pairs++] = i;
			}
		}

		// check, whether we tested points from all inner active facets
		// this is done only for testing of correctness of the test itself
		int nonchecked_faces = 0;
		FOR_ALL_FACETS(fid, &mesh) {
			bool ok = false;
			Facet *fac = mesh.facets[fid];
			if (fac->type == Facet::OUTER) continue;
			if (!(fac->ractive || fac->lactive)) continue;
			for (int i = 0; i < num_pairs - 1; i++) {
				Word_t fac_idx1 = mesh.get_facet_id(points[pairs[i]]->elm_idx, points[pairs[i]]->iface);
				Word_t fac_idx2 = mesh.get_facet_id(points[pairs[i + 1]]->elm_idx, points[pairs[i + 1]]->iface);
				if ((fac_idx1 == fid) || (fac_idx2 == fid)) {
					ok = true;
					break;
				}
			}

			if (!ok) nonchecked_faces++;
		}


		// loop over all basis functions
		for (int dof = 0; dof < ndofs; dof++) {
			printf("processing dof %d...\n", dof);

			// prepare solution correspondig to basis function with dof dof
			double sln_vector[ndofs];
			memset(sln_vector, 0, ndofs * sizeof(double));
			sln_vector[dof] = 1.0;
			sln.set_fe_solution(&space, sln_vector);

			double max_difference[3] = { 0., 0.0, 0.0 };
			double max_pt_x, max_pt_y, max_pt_z, max_val_1[3], max_val_2[3];
			Word_t max_elm_1, max_elm_2;

			RefMap *rm;
			double *nx, *ny, *nz;

			// loop over all pairs of points, that correspond to one point in the physical domain
			for(int pair = 0; pair < num_pairs; pair++) {
				int i = pairs[pair];

				Element *e1 = mesh.elements[points[i]->elm_idx];
				sln.set_active_element(e1);

				rm = sln.get_refmap();
				// FIXME: !!!
				QuadPt3D pt1(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z, 1.0);
				rm->calc_face_normal(points[i]->iface, 1, &pt1, nx, ny, nz);

				// TODO: improve me!
				double3 val1 = {
					sln.get_pt_value(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z, 0),
					sln.get_pt_value(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z, 1),
					sln.get_pt_value(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z, 2),
				};

				double3 tp1;
				calc_tan_proj<double, double>(nx[0], ny[0], nz[0], val1, tp1);

				delete [] nx;
				delete [] ny;
				delete [] nz;


				Element *e2 = mesh.elements[points[i + 1]->elm_idx];
				sln.set_active_element(e2);

				rm = sln.get_refmap();
				QuadPt3D pt2(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z, 1.0);
				rm->calc_face_normal(points[i]->iface, 1, &pt2, nx, ny, nz);

				// TODO: improve me!
				double val2[3] = {
					sln.get_pt_value(points[i + 1]->ref_x, points[i + 1]->ref_y, points[i + 1]->ref_z, 0),
					sln.get_pt_value(points[i + 1]->ref_x, points[i + 1]->ref_y, points[i + 1]->ref_z, 1),
					sln.get_pt_value(points[i + 1]->ref_x, points[i + 1]->ref_y, points[i + 1]->ref_z, 2),
				};

				double3 tp2;
				calc_tan_proj(nx[0], ny[0], nz[0], val2, tp2);

				delete [] nx;
				delete [] ny;
				delete [] nz;

				//value in this point should be the same, no matter from which element we go
				double difference[3] = {
					fabs(tp1[0] - tp2[0]),
					fabs(tp1[1] - tp2[1]),
					fabs(tp1[2] - tp2[2])
				};

//				printf("(% lf, % lf, % lf): [% lf, % lf, % lf] <=> [% lf, % lf, % lf]\n",
//					points[i]->phys_x, points[i]->phys_y, points[i]->phys_z,
//				    val1[0], val1[1], val1[2], val2[0], val2[1], val2[2]);

				double norm = sqrt(sqr(difference[0]) + sqr(difference[1] + sqr(difference[2])));
				double md_norm = sqrt(sqr(max_difference[0]) + sqr(max_difference[1] + sqr(max_difference[2])));
				if (norm > md_norm) {
					max_difference[0] = difference[0];
					max_difference[1] = difference[1];
					max_difference[2] = difference[2];

					max_pt_x = points[i]->phys_x;
					max_pt_y = points[i]->phys_y;
					max_pt_z = points[i]->phys_z;
					memcpy(max_val_1, val1, 3 * sizeof(double));
					memcpy(max_val_2, val2, 3 * sizeof(double));
//					max_val_2 = val2;
					max_elm_1 = points[i]->elm_idx;
					max_elm_2 = points[i + 1]->elm_idx;
				}
			}

			double md_norm = sqrt(sqr(max_difference[0]) + sqr(max_difference[1] + sqr(max_difference[2])));
			if (md_norm > TOLERANCE) {
				printf("base fn %d NOT continuous between elements %ld and %ld @ (% lf, % lf, % lf), "
					"max difference [%g, %g, %g] ([%.15g, %.15g, %.15g] <-> [%.15g, %.15g, %.15g])\n",
					 dof, max_elm_1, max_elm_2 , max_pt_x, max_pt_y, max_pt_z, max_difference[0], max_difference[1], max_difference[2],
					 max_val_1[0], max_val_1[1], max_val_1[2], max_val_2[0], max_val_2[1], max_val_2[2]);
				res = ERR_FAILURE;
			}

//			printf("--\n");
		}

		delete [] pairs;
		delete [] points;

		printf("continuity tested in %d points and %d inner faces with at least one active adjacent element were not tested\n", num_pairs, nonchecked_faces);
#endif

#if 0
#ifdef OUTPUT_DIR
		{
			char of_name[512];
			sprintf(of_name, "%s/sln.gmsh", OUTPUT_DIR);
			FILE *ofile = fopen(of_name, "w");
			if (ofile != NULL) {
				GmshOutputEngine output(ofile);
				output.out(&sln, "Uh");
				output.out(&exsln, "U");

				fclose(ofile);
			}
			else {
				ERROR("Cannot open '%s' for writing.", of_name);
			}
		}
#endif
#endif

#if 0
		// loop over all basis functions
//		int idof = 2; {
//		for (int idof = 0; idof < ndofs; idof++) {
		for (int idof = 0; idof < 4; idof++) {
//		for (int idof = 40; idof <= 41; idof++) {
			printf("dumping fn %d...\n", idof);

			// prepare solution corresponding to basis function with dof 'dof'
			double sln_vector[ndofs];
			memset(sln_vector, 0, ndofs * sizeof(double));
			sln_vector[idof] = 1.0;
			sln.set_fe_solution(&space, sln_vector);

//			sln.enable_transform(false);

#ifdef OUTPUT_DIR
			char of_name[512];
			sprintf(of_name, "%s/f%d.gmsh", OUTPUT_DIR, idof);
//			sprintf(of_name, "%s/f%d.vtk", OUTPUT_DIR, idof);
			FILE *ofile = fopen(of_name, "w");
			if (ofile != NULL) {
				GmshOutputEngine output(ofile);
//				VtkOutputEngine output(ofile, 7);
//				MagFilter mag(&sln);

				int mask = FN_VAL;
				if (idof == 0 || idof == 1) mask = FN_VAL_1;
				else if (idof == 2 || idof == 3) mask = FN_VAL_2;
//				mask = FN_VAL_0;

				output.out(&sln, "U", mask);
//				output.out(&sln, "U_dx", FN_DX);
//				output.out(&sln, "U_dy", FN_DY);
//				output.out(&sln, "U_dz", FN_DZ);

				fclose(ofile);
			}
			else {
				ERROR("Cannot open '%s' for writing.", of_name);
			}
#endif
		}
#endif

		if (res != ERR_SUCCESS) throw res;

		printf("Passed\n");
	}
	catch (int e) {
		res = e;
		printf("Failed\n");
	}

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	TRACE_END;

	return res;
}

