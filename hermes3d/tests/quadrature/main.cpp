// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
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
// main.cc
//
// Testing numerical quadrature
//
//

//
// TODO: test numerical quadrature on faces, edges
//
//

#include "config.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1

#define EPSILON										10e-12

#define countof(a) 									(sizeof(a)/sizeof(a[0]))


bool testPrint(bool value, const char *msg, bool correct) {
	printf("%s...", msg);
	if (value == correct) {
		printf("OK\n");
		return true;
	}
	else {
		printf("failed\n");
		return false;
	}
}

//
// Test quadrature
//

//
// 3D
//

typedef
	double (*fn3d_t)(double x, double y, double z);

struct TC3D {
	double exact;		// exact value of the integral
	int min_h_order;	// minimal horz order
	int min_v_order;	// minimal vert order
	int min_u_order;	// minimal vert order
	fn3d_t fn;			// function
	const char *fn_name;		// string representation of a function

	TC3D(fn3d_t f, double e, int min_h, int min_v, int min_u, const char *n) {
		exact = e;
		fn = f;
		fn_name = n;
		min_h_order = min_h;
		min_v_order = min_v;
		min_u_order = min_u;
	}
};

// function: f(x, y) = x^2 + y^2 + z^2 + x*y*z + x + y + z + 1
double fn_3d_1(double x, double y, double z) {
	return x*x + y*y + z*z + x*y*z + x + y + z + 1;
}

double fn_3d_2(double x, double y, double z) {
	return pow(x, 2) + pow(y, 3) + pow(z, 4) + x*y*z + x + y + z + 1;
}

double fn_3d_3(double x, double y, double z) {
	return pow(x, 3) + pow(y, 4) + pow(z, 2) + x*y*z + x + y + z + 1;
}

double fn_3d_4(double x, double y, double z) {
	return pow(x, 4) + pow(y, 2) + pow(z, 3) + x*y*z + x + y + z + 1;
}

// test 3d quadrature
int test_quadrature_3d_tetra(fn3d_t fn, double exact, int min_order, const char *fn_name) {
	printf("  * f(x,y) = %s", fn_name);

	// !!! std. quadrature works on std. reference domain !!!
	QuadStdTetra quad;
	for (int order = min_order; order <= H3D_MAX_QUAD_ORDER_TETRA; order++) {
		int np = quad.get_num_points(order);
		QuadPt3D *pt = quad.get_points(order);

		double integral = 0;
		for (int i = 0; i < np; i++) {
			integral += fn(pt[i].x, pt[i].y, pt[i].z) * pt[i].w;
		}

		double err = fabs(exact - integral);
		if (err >= EPSILON) {
			printf(" ... failed for order %d, integral = %lf, expected = %lf\n", order, integral, exact);
			return ERROR_FAILURE;
		}
	}

	printf(" ... OK\n");

	return ERROR_SUCCESS;
}

int test_quadrature_3d_hex(fn3d_t fn, double exact, int min_h, int min_v, int min_u, const char *fn_name) {
	printf("  * f(x,y) = %s", fn_name);

	// !!! std. quadrature works on std. reference domain !!!
	QuadStdHex quad;
	for (int horder = min_h; horder <= H3D_MAX_QUAD_ORDER; horder++) {
		for (int vorder = min_v; vorder <= H3D_MAX_QUAD_ORDER; vorder++) {
			for (int uorder = min_u; uorder <= H3D_MAX_QUAD_ORDER; uorder++) {
				order3_t order(horder, vorder, uorder);

				int np = quad.get_num_points(order);
				QuadPt3D *pt = quad.get_points(order);

				double integral = 0;
				for (int i = 0; i < np; i++) {
					integral += fn(pt[i].x, pt[i].y, pt[i].z) * pt[i].w;
				}

				double err = fabs(exact - integral);
//				printf("  * order (h = %d, v = %d, u = %d)", horder, vorder, uorder);
				if (err >= EPSILON) {
					printf(" ... failed for order (h = %d, v = %d, u = %d), integral = %lf, expected = %lf (diff = %e)\n",
						horder, vorder, uorder, integral, exact, fabs(integral - exact));
					return ERROR_FAILURE;
				}
			}
		}
	}

	printf(" ... OK\n");

	return ERROR_SUCCESS;
}

int test_quadrature_3d_hex_surf(fn3d_t fn, double exact, int min_h, int min_v, int min_u, const char *fn_name) {
	printf("  * f(x,y) = %s", fn_name);

	// !!! std. quadrature works on std. reference domain !!!
	QuadStdHex quad;
	for (int horder = min_h; horder <= H3D_MAX_QUAD_ORDER; horder++) {
		for (int vorder = min_v; vorder <= H3D_MAX_QUAD_ORDER; vorder++) {
			for (int uorder = min_u; uorder <= H3D_MAX_QUAD_ORDER; uorder++) {
				order2_t face_order[] = {
					order2_t(vorder, uorder),
					order2_t(vorder, uorder),
					order2_t(horder, uorder),
					order2_t(horder, uorder),
					order2_t(horder, vorder),
					order2_t(horder, vorder)
				};

				double integral = 0;
				for (int face = 0; face < Hex::NUM_FACES; face++) {
					order2_t order = face_order[face];
					int np = quad.get_face_num_points(face, order);
					QuadPt3D *pt = quad.get_face_points(face, order);

					for (int i = 0; i < np; i++)
						integral += fn(pt[i].x, pt[i].y, pt[i].z) * pt[i].w;
				}

				double err = fabs(exact - integral);
//				printf("  * order (h = %d, v = %d, u = %d)", horder, vorder, uorder);
				if (err >= EPSILON) {
					printf(" ... failed for order (h = %d, v = %d, u = %d), integral = %lf, expected = %lf (diff = %e)\n",
						horder, vorder, uorder, integral, exact, fabs(integral - exact));
					return ERROR_FAILURE;
				}
			}
		}
	}

	printf(" ... OK\n");

	return ERROR_SUCCESS;
}

int test_quadrature_3d() {
	int ret = ERROR_SUCCESS;

	// hexs ///
	if (get_quadrature(MODE_HEXAHEDRON) != NULL) {
		printf("\n");
		printf("- Testing 3D quadrature (hex) -----\n");

		TC3D fn_hex[] = {
			TC3D(fn_3d_1, 16.0, 2, 2, 2, "x^2 + y^2 + z^2 + x*y*z + x + y + z + 1"),
			TC3D(fn_3d_2, 184.0/15.0, 2, 3, 4, "x^2 + y^3 + z^4 + x*y*z + x + y + z + 1"),
			TC3D(fn_3d_3, 184.0/15.0, 3, 4, 2, "x^3 + y^4 + z^2 + x*y*z + x + y + z + 1"),
			TC3D(fn_3d_4, 184.0/15.0, 4, 2, 3, "x^4 + y^2 + z^3 + x*y*z + x + y + z + 1")
		};

		for (unsigned int i = 0; i < countof(fn_hex); i++) {
			if ((ret = test_quadrature_3d_hex(fn_hex[i].fn, fn_hex[i].exact, fn_hex[i].min_h_order, fn_hex[i].min_v_order, fn_hex[i].min_u_order, fn_hex[i].fn_name)) != ERROR_SUCCESS)
				return ret;
		}

		printf("\n");
		printf("- Testing 3D quadrature (hex) - surf -----\n");

		TC3D fn_hex_surf[] = {
			TC3D(fn_3d_1, 64.0, 2, 2, 2, "x^2 + y^2 + z^2 + x*y*z + x + y + z + 1")
		};

		for (unsigned int i = 0; i < countof(fn_hex_surf); i++) {
			if ((ret = test_quadrature_3d_hex_surf(fn_hex_surf[i].fn, fn_hex_surf[i].exact, fn_hex_surf[i].min_h_order, fn_hex_surf[i].min_v_order, fn_hex_surf[i].min_u_order, fn_hex_surf[i].fn_name)) != ERROR_SUCCESS)
				return ret;
		}
	}

	// tetras ///
	if (get_quadrature(MODE_TETRAHEDRON) != NULL) {
		printf("\n");
		printf("- Testing 3D quadrature (tetra) -----\n");

		TC3D fn_tetra[] = {
			TC3D(fn_3d_1, 8.0/9.0, 3, 0, 0, "x^2 + y^2 + z^2 + x*y*z + x + y + z + 1")
		};

		for (unsigned int i = 0; i < countof(fn_tetra); i++) {
			if ((ret = test_quadrature_3d_tetra(fn_tetra[i].fn, fn_tetra[i].exact, fn_tetra[i].min_h_order, fn_tetra[i].fn_name)) != ERROR_SUCCESS)
				return ret;
		}
	}

	// TODO: prisms

	return ERROR_SUCCESS;
}

//
// main
//

int main() {
	set_verbose(false);

//	TRACE_START("trace.txt");
	DEBUG_OUTPUT_OFF;
	SET_VERBOSE_LEVEL(0);

	int ret = ERROR_SUCCESS;

	// test 3D quadrature
	if ((ret = test_quadrature_3d()) != ERROR_SUCCESS)
		return ret;

//	TRACE_END;

	return ret;
}
