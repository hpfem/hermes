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

#include "config.h"
#include "common.h"
#include <hermes3d.h>
#include <common/trace.h>
#include <common/error.h>

double hh = 1e-4;
double tol = hh * hh * 20.;

int np = 5;
double ptx[5] = { -0.1, 0.6, -0.4, 0.5, 0.3 };
double pty[5] = { 0., -0.94, 0.42, 0.13, 0.7 };
double ptz[5] = { 0.8, -0.4, 0.2, -0.7, -0.1 };

double dif, maxdif;

bool test(Shapeset *ss, int fnum) {
	_F_
	scalar val, valph, valder;

	bool passed = true;

	for (int point = 0; point < np; point++) {
		for (int comp = 0; comp < ss->get_num_components(); comp++) {
			// dx
			val = ss->get_fn_value(fnum, ptx[point], pty[point], ptz[point], comp);
			valph = ss->get_fn_value(fnum, ptx[point] + hh, pty[point], ptz[point], comp);
			valder = ss->get_dx_value(fnum, ptx[point] + hh / 2., pty[point], ptz[point], comp);
			if ((dif = fabs(valph - val - hh * valder)) > tol) {
				passed = false;
				printf("GRADTEST ERROR @(%lf, %lf, %lf), function %d, component %d dx is %lf, should be %lf\n",
				       ptx[point], pty[point], ptz[point], fnum, comp, valder, (valph - val) / hh);
			}
			if (dif > maxdif) maxdif = dif;

			// dy
			val = ss->get_fn_value(fnum, ptx[point], pty[point], ptz[point], comp);
			valph = ss->get_fn_value(fnum, ptx[point], pty[point] + hh, ptz[point], comp);
			valder = ss->get_dy_value(fnum, ptx[point], pty[point] + hh / 2., ptz[point], comp);
			if ((dif = fabs(valph - val - hh * valder)) > tol) {
				passed = false;
				printf("GRADTEST ERROR @(%lf, %lf, %lf), function %d, component %d dy is %lf, should be %lf\n",
				       ptx[point], pty[point], ptz[point], fnum, comp, valder, (valph - val) / hh);
			}
			if (dif > maxdif) maxdif = dif;

			// dz
			val = ss->get_fn_value(fnum, ptx[point], pty[point], ptz[point], comp);
			valph = ss->get_fn_value(fnum, ptx[point], pty[point], ptz[point] + hh, comp);
			valder = ss->get_dz_value(fnum, ptx[point], pty[point], ptz[point] + hh / 2., comp);
			if ((dif = fabs(valph - val - hh * valder)) > tol) {
				passed = false;
				printf("GRADTEST ERROR @(%lf, %lf, %lf), function %d, component %d dz is %lf, should be %lf\n",
				       ptx[point], pty[point], ptz[point], fnum, comp, valder, (valph - val) / hh);
			}
			if (dif > maxdif) maxdif = dif;

		}
	}
	return passed;
}

bool test_gradients_directly(Shapeset *ss) {
	_F_
	printf("V. direct check of the gradient values\n");

	maxdif = 0.;
	bool passed = true;
	int *indices, ii;

	order1_t order = H3D_MAX_ELEMENT_ORDER;
	for (int ie = 0; ie < 12; ie++) {
		for (int ori = 0; ori < 2; ori++) {
			indices = ss->get_edge_indices(ie, ori, order);
			for (ii = 0; ii < ss->get_num_edge_fns(order); ii++) {
				printf("."); fflush(stdout);
				if (!test(ss, indices[ii])) passed = false;
			}
		}
	}

	order2_t order2(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER);
	for (int ic = 0; ic < 6; ic++) {
		for (int ori = 0; ori < 8; ori++) {
			indices = ss->get_face_indices(ic, ori, order2);
			for (ii = 0; ii < ss->get_num_face_fns(order2); ii++) {
				printf("."); fflush(stdout);
				if (!test(ss, indices[ii])) passed = false;
			}
		}
	}

	order3_t order3(H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER, H3D_MAX_ELEMENT_ORDER);
	indices = ss->get_bubble_indices(order3);
	for (ii = 0; ii < ss->get_num_bubble_fns(order3); ii++) {
		printf("."); fflush(stdout);
		if (!test(ss, indices[ii])) passed = false;
	}

	printf("\n");
	printf("maximal difference is %g, which is %g * h^2\n", maxdif, maxdif / hh / hh);

	return passed;
}

