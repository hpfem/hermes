#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// Testing numerical quadrature.
// TODO: Test numerical quadrature on faces, edges.

// The error should be smaller than this EPS.
#define EPS								10e-10F

#define countof(a) (sizeof(a)/sizeof(a[0]))


bool testPrint(bool value, const char *msg, bool correct) {
	info("%s...", msg);
	if (value == correct) {
		info("OK.");
		return true;
	}
	else {
		info("failed.");
		return false;
	}
}

// Test quadrature.
// 3D.

typedef
	double (*fn3d_t)(double x, double y, double z);

struct TC3D {
	double exact;		// exact value of the integral
	int min_h_order;	// minimal horz order
	int min_v_order;	// minimal vert order
	int min_u_order;	// minimal vert order
	fn3d_t fn;		// function
	const char *fn_name;	// string representation of a function

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

// Test 3d quadrature.
int test_quadrature_3d_tetra(fn3d_t fn, double exact, int min_order, const char *fn_name) {
	info("  * f(x,y) = %s", fn_name);

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
		if (err >= EPS) {
			info(" ... failed for order %d, integral = %lf, expected = %lf.", order, integral, exact);
			return ERR_FAILURE;
		}
	}

	info(" ... OK.");

	return ERR_SUCCESS;
}

int test_quadrature_3d_hex(fn3d_t fn, double exact, int min_h, int min_v, int min_u, const char *fn_name) {
	info("  * f(x,y) = %s", fn_name);

	// !!! std. quadrature works on std. reference domain !!!
	QuadStdHex quad;
	for (int horder = min_h; horder <= H3D_MAX_QUAD_ORDER; horder++) {
		for (int vorder = min_v; vorder <= H3D_MAX_QUAD_ORDER; vorder++) {
			for (int uorder = min_u; uorder <= H3D_MAX_QUAD_ORDER; uorder++) {
				Ord3 order(horder, vorder, uorder);

				int np = quad.get_num_points(order);
				QuadPt3D *pt = quad.get_points(order);

				double integral = 0;
				for (int i = 0; i < np; i++) {
					integral += fn(pt[i].x, pt[i].y, pt[i].z) * pt[i].w;
				}

				double err = fabs(exact - integral);
				if (err >= EPS) {
					info(" ... failed for order (h = %d, v = %d, u = %d), integral = %lf, expected = %lf (diff = %e).",
						horder, vorder, uorder, integral, exact, fabs(integral - exact));
					return ERR_FAILURE;
				}
			}
		}
	}

	info(" ... OK.");

	return ERR_SUCCESS;
}

int test_quadrature_3d_hex_surf(fn3d_t fn, double exact, int min_h, int min_v, int min_u, const char *fn_name) {
	info("  * f(x,y) = %s", fn_name);

	QuadStdHex quad;
	for (int horder = min_h; horder <= H3D_MAX_QUAD_ORDER; horder++) {
		for (int vorder = min_v; vorder <= H3D_MAX_QUAD_ORDER; vorder++) {
			for (int uorder = min_u; uorder <= H3D_MAX_QUAD_ORDER; uorder++) {
				Ord2 face_order[] = {
					Ord2(vorder, uorder),
					Ord2(vorder, uorder),
					Ord2(horder, uorder),
					Ord2(horder, uorder),
					Ord2(horder, vorder),
					Ord2(horder, vorder)
				};

				double integral = 0;
				for (int face = 0; face < Hex::NUM_FACES; face++) {
					Ord2 order = face_order[face];
					int np = quad.get_face_num_points(face, order);
					QuadPt3D *pt = quad.get_face_points(face, order);

					for (int i = 0; i < np; i++)
						integral += fn(pt[i].x, pt[i].y, pt[i].z) * pt[i].w;
				}

				double err = fabs(exact - integral);
				if (err >= EPS) {
					info(" ... failed for order (h = %d, v = %d, u = %d), integral = %lf, expected = %lf (diff = %e).",
						horder, vorder, uorder, integral, exact, fabs(integral - exact));
					return ERR_FAILURE;
				}
			}
		}
	}

	info(" ... OK.");

	return ERR_SUCCESS;
}

int test_quadrature_3d() {
	int ret = ERR_SUCCESS;

	// Hexs.
	if (get_quadrature(HERMES_MODE_HEX) != NULL) {
		info("- Testing 3D quadrature (hex) -----");

		TC3D fn_hex[] = {
			TC3D(fn_3d_1, 16.0, 2, 2, 2, "x^2 + y^2 + z^2 + x*y*z + x + y + z + 1"),
			TC3D(fn_3d_2, 184.0/15.0, 2, 3, 4, "x^2 + y^3 + z^4 + x*y*z + x + y + z + 1"),
			TC3D(fn_3d_3, 184.0/15.0, 3, 4, 2, "x^3 + y^4 + z^2 + x*y*z + x + y + z + 1"),
			TC3D(fn_3d_4, 184.0/15.0, 4, 2, 3, "x^4 + y^2 + z^3 + x*y*z + x + y + z + 1")
		};

		for (unsigned int i = 0; i < countof(fn_hex); i++) {
			if ((ret = test_quadrature_3d_hex(fn_hex[i].fn, fn_hex[i].exact, fn_hex[i].min_h_order, fn_hex[i].min_v_order, fn_hex[i].min_u_order, fn_hex[i].fn_name)) != ERR_SUCCESS)
				return ret;
		}

		info("- Testing 3D quadrature (hex) - surf -----");

		TC3D fn_hex_surf[] = {
			TC3D(fn_3d_1, 64.0, 2, 2, 2, "x^2 + y^2 + z^2 + x*y*z + x + y + z + 1")
		};

		for (unsigned int i = 0; i < countof(fn_hex_surf); i++) {
			if ((ret = test_quadrature_3d_hex_surf(fn_hex_surf[i].fn, fn_hex_surf[i].exact, fn_hex_surf[i].min_h_order, fn_hex_surf[i].min_v_order, fn_hex_surf[i].min_u_order, fn_hex_surf[i].fn_name)) != ERR_SUCCESS)
				return ret;
		}
	}

	// Tetras.
	if (get_quadrature(HERMES_MODE_TET) != NULL) {
		info("- Testing 3D quadrature (tetra) -----");

		TC3D Funcetra[] = {
			TC3D(fn_3d_1, 8.0/9.0, 3, 0, 0, "x^2 + y^2 + z^2 + x*y*z + x + y + z + 1")
		};

		for (unsigned int i = 0; i < countof(Funcetra); i++) {
			if ((ret = test_quadrature_3d_tetra(Funcetra[i].fn, Funcetra[i].exact, Funcetra[i].min_h_order, Funcetra[i].fn_name)) != ERR_SUCCESS)
				return ret;
		}
	}

	// TODO: Prisms.

	return ERR_SUCCESS;
}

int main() {
	int ret = ERR_SUCCESS;

	// Test 3D quadrature.
	ret = test_quadrature_3d();
	return ret;
if (ret == ERR_SUCCESS) {
    info("Success!");
    return ERR_SUCCESS;
  }
  else {
    info("Failure!");
    return ERR_FAILURE;
  }
}
