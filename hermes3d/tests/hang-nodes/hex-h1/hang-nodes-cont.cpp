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

#define BEGIN_BLOCK						{
#define END_BLOCK						}

//#define DIRICHLET
//#define NEWTON

//#define X2_Y2_Z2
//#define X3_Y3_Z3
//#define XM_YN_ZO
#define XM_YN_ZO_2

int m = 2, n = 2, o = 2;

template<typename S, typename T>
T fnc(S x, S y, S z) {
#ifdef XM_YN_ZO
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
#elif defined XM_YN_ZO_2
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 2) * z + pow(z, 4);
#elif defined X2_Y2_Z2
	return x*x + y*y + z*z;
#elif defined X3_Y3_Z3
	return x*x*x + y*y*y + z*z*z;
#endif
}

template<typename S, typename T>
T dfnc(S x, S y, S z) {
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
#elif defined X3_Y3_Z3
	return -6.0 * x - 6.0 * y - 6.0 * z;
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
#elif defined X3_Y3_Z3
	dx = 3 * x*x;
	dy = 3 * y*y;
	dz = 3 * z*z;
#endif

	return fnc<double, double>(x, y, z);
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
	return fnc<double, scalar>(x, y, z);
#else
	return 0;
#endif
}

template<typename f_t, typename res_t>
res_t bilinear_form(int np, double *jwt, fn_t<res_t> *u_ext[], fn_t<f_t> *fu, fn_t<f_t> *fv, geom_t<f_t> *e, user_data_t<res_t> *ud) {
	return int_grad_u_grad_v<f_t, res_t>(np, jwt, fu, fv, e);
}

template<typename f_t, typename res_t>
res_t bilinear_form_surf(int np, double *jwt, fn_t<res_t> *u_ext[], fn_t<f_t> *fu, fn_t<f_t> *fv, geom_t<f_t> *e, user_data_t<res_t> *ud) {
	return int_u_v<f_t, res_t>(np, jwt, fu, fv, e);
}

template<typename f_t, typename res_t>
res_t linear_form(int np, double *jwt, fn_t<res_t> *u_ext[], fn_t<f_t> *fv, geom_t<f_t> *e, user_data_t<res_t> *ud) {
	return int_F_v<f_t, res_t>(np, jwt, dfnc, fv, e);
}

template<typename f_t, typename res_t>
res_t linear_form_surf(int np, double *jwt, fn_t<res_t> *u_ext[], fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *ud) {
	res_t result = 0;
#ifdef XM_YN_ZO
	for (int i = 0; i < np; i++) {
		res_t dx = m * pow(e->x[i], m - 1) * pow(e->y[i], n) * pow(e->z[i], o) + 2 * e->x[i] * pow(e->y[i], 3) - 3 * pow(e->x[i], 2) * e->z[i];
		res_t dy = n * pow(e->x[i], m) * pow(e->y[i], n - 1) * pow(e->z[i], o) + 3 * pow(e->x[i], 2) * pow(e->y[i], 2);
		res_t dz = o * pow(e->x[i], m) * pow(e->y[i], n) * pow(e->z[i], o - 1) - pow(e->x[i], 3) + 4 * pow(e->z[i], 3);
		result += jwt[i] * (v->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc<f_t, res_t>(e->x[i], e->y[i], e->z[i])));
	}
#elif defined XM_YN_ZO_2
	for (int i = 0; i < np; i++) {
		res_t dx = m * pow(e->x[i], m-1) * pow(e->y[i], n) * pow(e->z[i], o) + 2 * e->x[i] * pow(e->y[i], 3) - 2 * e->x[i] * e->z[i];
		res_t dy = n * pow(e->x[i], m) * pow(e->y[i], n-1) * pow(e->z[i], o) + 3 * pow(e->x[i], 2) * pow(e->y[i], 2);
		res_t dz = o * pow(e->x[i], m) * pow(e->y[i], n) * pow(e->z[i], o-1) - pow(e->x[i], 2) + 4 * pow(e->z[i], 3);
		result += jwt[i] * (v->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc<f_t, res_t>(e->x[i], e->y[i], e->z[i])));
	}
#elif defined X2_Y2_Z2
	for (int i = 0; i < np; i++) {
		res_t dx = 2 * e->x[i];
		res_t dy = 2 * e->y[i];
		res_t dz = 2 * e->z[i];
		result += jwt[i] * (v->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc<f_t, res_t>(e->x[i], e->y[i], e->z[i])));
	}
#elif defined X3_Y3_Z3
	for (int i = 0; i < np; i++) {
		res_t dx = 3 * e->x[i] * e->x[i];
		res_t dy = 3 * e->y[i] * e->y[i];
		res_t dz = 3 * e->z[i] * e->z[i];
		result += jwt[i] * (v->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc<f_t, res_t>(e->x[i], e->y[i], e->z[i])));
	}
#endif
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
		for (int face = 0; face < Hex::NUM_FACES; face++)
			for (int order = 0; order < NUM_RULES; order++)
				delete [] face_tables[face][order];
		delete [] face_tables;
	}
};

struct Point {
	double ref_x, ref_y, ref_z;
	double phys_x, phys_y, phys_z;
	Word_t elm_idx, fac_idx;

	Point(int idx, int f_idx, double rx, double ry, double rz, double px, double py, double pz) {
		elm_idx = idx; fac_idx = f_idx;
		ref_x = rx;  ref_y = ry;  ref_z = rz;
		phys_x = px; phys_y = py; phys_z = pz;
	}
};

typedef
	int (*compfn)(const void*, const void*);

int compare(Point **pt1, Point **pt2) {
	double val1 = 1000000. * (*pt1)->phys_x + 1000. * (*pt1)->phys_y + (*pt1)->phys_z;
	double val2 = 1000000. * (*pt2)->phys_x + 1000. * (*pt2)->phys_y + (*pt2)->phys_z;

	if (val1 < val2) return -1;
	else if (val1 > val2) return 1;
	else return 0;
}

const double EPS = 1e-13;
const double TOLERANCE = 1e-10;

bool equal(Point *pt1, Point *pt2) {
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

	TRACE_START("trace.txt");
	DEBUG_OUTPUT_ON;
	SET_VERBOSE_LEVEL(0);

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
		const char *of_name = OUTPUT_DIR "/mesh.vtk";
		FILE *ofile = fopen(of_name, "w");
		if (ofile != NULL) {
			VtkOutputEngine output(ofile);
			output.out(&mesh);
			fclose(ofile);
		}
		else {
			warning("Can not open '%s' for writing.", of_name);
		}
	END_BLOCK
#endif

	H1ShapesetLobattoHex shapeset;
//	shapeset.preload_products();

	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_essential_bc_values(essential_bc_values);

#ifdef XM_YN_ZO
	order3_t o(3, 3, 4);
#elif defined XM_YN_ZO_2
	order3_t o(2, 3, 4);
#elif defined X2_Y2_Z2
	order3_t o(2, 2, 2);
#elif defined X3_Y3_Z3
	order3_t o(3, 3, 3);
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

	WeakForm wf(1);
#ifdef DIRICHLET
	wf.add_matrix_form(0, 0, bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM);
	wf.add_vector_form(0, linear_form<double, scalar>, linear_form<ord_t, ord_t>);
#elif defined NEWTON
	wf.add_matrix_form(0, 0, bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM);
	wf.add_matrix_form_surf(0, 0, bilinear_form_surf<double, scalar>, bilinear_form_surf<ord_t, ord_t>);
	wf.add_vector_form(0, linear_form<double, scalar>, linear_form<ord_t, ord_t>);
	wf.add_vector_form_surf(0, linear_form_surf<double, scalar>, linear_form_surf<ord_t, ord_t>);
#endif

	// assemble stiffness matrix
	LinearProblem lp(&wf);
	lp.set_space(&space);

	lp.assemble(&mat, &rhs);

#ifdef OUTPUT_DIR
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
		// solve the stiffness matrix
		bool solved = solver.solve();

		if (!solved) throw ERR_FAILURE;

		Solution sln(&mesh);
		sln.set_fe_solution(&space, solver.get_solution());

		{
			double *s = solver.get_solution();
			for (int i = 0; i < ndofs; i++)
				printf("x[% 3d] = % lf\n", i, s[i]);
		}

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

		//
		//
		// the main code starts here
		//
		//
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
					Word_t fac_idx = mesh.get_facet_id(e, iface);

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
						points[ipt++] = new Point(idx, fac_idx,
							quad_pts[pt].x, quad_pts[pt].y, quad_pts[pt].z,
							phys_x[pt], phys_y[pt], phys_z[pt]);
					}
				}
			}
		}

		// sort points according to first phys_x, then phys_y and phys_z
		// it means, that two points, with almost identical physical coordinates
		// (even though from different elements) will be next to each other in the array
		qsort((void *) points, num_points, sizeof(Point*), (compfn)compare);

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
			for (int i = 0; i < num_pairs; i++) {
				if ((points[pairs[i]]->fac_idx == fid) || (points[pairs[i] + 1]->fac_idx == fid)) {
					ok = true;
					break;
				}
			}

			if (!ok) nonchecked_faces++;
		}


		// loop over all basis functions
		for (int dof = 0; dof < ndofs; dof++) {
			printf("processing dof %d...", dof); fflush(stdout);

			// prepare solution correspondig to basis function with dof dof
			double sln_vector[ndofs];
			memset(sln_vector, 0, ndofs * sizeof(double));
			sln_vector[dof] = 1.0;
			sln.set_fe_solution(&space, sln_vector);

			double max_difference = 0.;
			double max_pt_x, max_pt_y, max_pt_z, max_val_1, max_val_2;
			Word_t max_elm_1, max_elm_2;

			// loop over all pairs of points, that correspond to one point in the physical domain
			for (int pair = 0; pair < num_pairs; pair++) {
				int i = pairs[pair];
				Element *e1 = mesh.elements[points[i]->elm_idx];
				sln.set_active_element(e1);
				// TODO: improve me!
				double val1 = sln.get_pt_value(points[i]->ref_x, points[i]->ref_y, points[i]->ref_z);

				Element *e2 = mesh.elements[points[i + 1]->elm_idx];
				sln.set_active_element(e2);
				// TODO: improve me!
				double val2 = sln.get_pt_value(points[i + 1]->ref_x, points[i + 1]->ref_y, points[i + 1]->ref_z);

				//value in this point should be the same, no matter from which element we go
				double difference = fabs(val1 - val2);
				if (difference > max_difference) {
					max_difference = difference;
					max_pt_x = points[i]->phys_x;
					max_pt_y = points[i]->phys_y;
					max_pt_z = points[i]->phys_z;
					max_val_1 = val1;
					max_val_2 = val2;
					max_elm_1 = points[i]->elm_idx;
					max_elm_2 = points[i + 1]->elm_idx;
				}
			}

			if (max_difference > TOLERANCE) {
				printf("failed\nbase fn %d NOT continuous between elements %ld and %ld @ (% lf, % lf, % lf), max difference %g (%.15g <-> %.15g)\n",
						 dof, max_elm_1, max_elm_2 , max_pt_x, max_pt_y, max_pt_z, max_difference, max_val_1, max_val_2);
				res = ERR_FAILURE;
			}
			else {
				printf("ok\n");
			}

		}

		for (int i = 0; i < num_points; i++) delete points[i];
		delete [] points;
		delete [] pairs;

		printf("continuity tested in %d points and %d inner faces with at least one active adjacent element were not tested\n", num_pairs, nonchecked_faces);
#endif

#if 0
		// loop over all basis functions
		for (int dof = 0; dof <= ndofs; dof++) {
			printf("dumping fn %d...\n", dof);

			// prepare solution corresponding to basis function with dof 'dof'
			double sln_vector[ndofs];
			memset(sln_vector, 0, ndofs * sizeof(double));
			sln_vector[dof] = 1.0;
			sln.set_fe_solution(&space, sln_vector);

#ifdef OUTPUT_DIR
			char of_name[512];
			sprintf(of_name, "%s/f%d.gmsh", OUTPUT_DIR, dof);
			FILE *ofile = fopen(of_name, "w");
			if (ofile != NULL) {
				GmshOutputEngine output(ofile);
				output.out(&sln, "U");

				fclose(ofile);
			}
			else {
				warning("Can not open '%s' for writing.", of_name);
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

