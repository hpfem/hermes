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

#include "h3dconfig.h"
#include "forms.h"
#include <common/callstack.h>
#include "integrals/hcurl.h"

geom_t<ord_t> init_geom(int marker) {
	_F_

	geom_t<ord_t> e;
	e.marker = marker;

	static ord_t x[] = { ord_t(1) };
	static ord_t y[] = { ord_t(1) };
	static ord_t z[] = { ord_t(1) };

	static ord_t nx[] = { ord_t(1) };
	static ord_t ny[] = { ord_t(1) };
	static ord_t nz[] = { ord_t(1) };

	static ord_t tx[] = { ord_t(1) };
	static ord_t ty[] = { ord_t(1) };
	static ord_t tz[] = { ord_t(1) };

	e.x = x; e.y = y; e.z = z;
	e.nx = nx; e.ny = ny; e.nz = nz;
	e.tx = tx; e.ty = ty; e.tz = tz;
	return e;
}

geom_t<double> init_geom(int marker, RefMap *rm, const int np, const QuadPt3D *pt) {
	_F_

	geom_t<double> e;
	e.marker = marker;
	e.x = rm->get_phys_x(np, pt);
	e.y = rm->get_phys_y(np, pt);
	e.z = rm->get_phys_z(np, pt);
	return e;
}

geom_t<double> init_geom(int marker, RefMap *rm, int iface, const int np, const QuadPt3D *pt) {
	_F_

	geom_t<double> e;
	e.marker = marker;
	e.x = rm->get_phys_x(np, pt);
	e.y = rm->get_phys_y(np, pt);
	e.z = rm->get_phys_z(np, pt);
	rm->calc_face_normal(iface, np, pt, e.nx, e.ny, e.nz);
	return e;
}

void free_geom(geom_t<double> *e) {
	delete [] e->x;
	delete [] e->y;
	delete [] e->z;

	delete [] e->nx;
	delete [] e->ny;
	delete [] e->nz;
}

fn_t<ord_t> init_fn(const order3_t &order) {
	int o = order.get_ord();
	ord_t *d = new ord_t(o);

	fn_t<ord_t> f;
	f.fn = d;
	f.dx = f.dy = f.dz = d;
	f.fn0 = f.fn1 = f.fn2 = d;
	f.dx0 = f.dx1 = f.dx2 = d;
	f.dy0 = f.dy1 = f.dy2 = d;
	f.dz0 = f.dz1 = f.dz2 = d;
	f.curl0 = f.curl1 = f.curl2 = d;
	return f;
}

sfn_t *init_fn(ShapeFunction *shfn, RefMap *rm, const int np, const QuadPt3D *pt) {
	_F_

	sfn_t *u = new sfn_t; MEM_CHECK(u);
	u->nc = shfn->get_num_components();
	shfn->precalculate(np, pt, FN_DEFAULT);
	if (u->nc == 1) {
		u->fn = new double [np]; MEM_CHECK(u->fn);
		u->dx = new double [np]; MEM_CHECK(u->dx);
		u->dy = new double [np]; MEM_CHECK(u->dy);
		u->dz = new double [np]; MEM_CHECK(u->dz);

		double *fn = shfn->get_fn_values();
		double *dx = shfn->get_dx_values();
		double *dy = shfn->get_dy_values();
		double *dz = shfn->get_dz_values();
		double3x3 *m = rm->get_inv_ref_map(np, pt);
		for (int i = 0; i < np; i++) {
			u->fn[i] = fn[i];
			u->dx[i] = (dx[i] * m[i][0][0] + dy[i] * m[i][0][1] + dz[i] * m[i][0][2]);
			u->dy[i] = (dx[i] * m[i][1][0] + dy[i] * m[i][1][1] + dz[i] * m[i][1][2]);
			u->dz[i] = (dx[i] * m[i][2][0] + dy[i] * m[i][2][1] + dz[i] * m[i][2][2]);
		}
		delete [] m;
	}
	else if (u->nc == 3) {
		u->fn0 = new double [np]; MEM_CHECK(u->fn0);
		u->fn1 = new double [np]; MEM_CHECK(u->fn1);
		u->fn2 = new double [np]; MEM_CHECK(u->fn2);

		double *fn[3];
		for (int c = 0; c < 3; c++)
			fn[c] = shfn->get_fn_values(c);

		double3x3 *irm = rm->get_inv_ref_map(np, pt);
		for (int i = 0; i < np; i++) {
			u->fn0[i] = fn[0][i] * irm[i][0][0] + fn[1][i] * irm[i][0][1] + fn[2][i] * irm[i][0][2];
			u->fn1[i] = fn[0][i] * irm[i][1][0] + fn[1][i] * irm[i][1][1] + fn[2][i] * irm[i][1][2];
			u->fn2[i] = fn[0][i] * irm[i][2][0] + fn[1][i] * irm[i][2][1] + fn[2][i] * irm[i][2][2];
		}
		delete [] irm;
	}

	if (shfn->get_type() == Hcurl) {
		u->curl0 = new double [np]; MEM_CHECK(u->curl0);
		u->curl1 = new double [np]; MEM_CHECK(u->curl1);
		u->curl2 = new double [np]; MEM_CHECK(u->curl2);

		double *dx[3], *dy[3], *dz[3];
		for (int c = 0; c < 3; c++) {
			dx[c] = shfn->get_dx_values(c);
			dy[c] = shfn->get_dy_values(c);
			dz[c] = shfn->get_dz_values(c);
		}

		// NOTE: are we able to work with transformed jacobian here?
		double *jac = rm->get_jacobian(np, pt, false);
		double3x3 *m = rm->get_ref_map(np, pt);
		for (int i = 0; i < np; i++) {
			double curl[3] = { dy[2][i] - dz[1][i], dz[0][i] - dx[2][i], dx[1][i] - dy[0][i] };
			u->curl0[i] = (curl[0] * m[i][0][0] + curl[1] * m[i][0][1] + curl[2] * m[i][0][2]) / jac[i];
			u->curl1[i] = (curl[0] * m[i][1][0] + curl[1] * m[i][1][1] + curl[2] * m[i][1][2]) / jac[i];
			u->curl2[i] = (curl[0] * m[i][2][0] + curl[1] * m[i][2][1] + curl[2] * m[i][2][2]) / jac[i];
		}

		delete [] m;
		delete [] jac;
	}

	return u;
}


sfn_t *init_fn(ShapeFunction *shfn, RefMap *rm, int iface, const int np, const QuadPt3D *pt) {
	_F_

	sfn_t *u = new sfn_t; MEM_CHECK(u);
	u->nc = shfn->get_num_components();
	shfn->precalculate(np, pt, FN_DEFAULT);
	if (u->nc == 1) {
		u->fn = new double [np]; MEM_CHECK(u->fn);
		u->dx = new double [np]; MEM_CHECK(u->dx);
		u->dy = new double [np]; MEM_CHECK(u->dy);
		u->dz = new double [np]; MEM_CHECK(u->dz);

		double *fn = shfn->get_fn_values();
		double *dx = shfn->get_dx_values();
		double *dy = shfn->get_dy_values();
		double *dz = shfn->get_dz_values();
		double3x3 *m = rm->get_inv_ref_map(np, pt);
		for (int i = 0; i < np; i++) {
			u->fn[i] = fn[i];
			u->dx[i] = (dx[i] * m[i][0][0] + dy[i] * m[i][0][1] + dz[i] * m[i][0][2]);
			u->dy[i] = (dx[i] * m[i][1][0] + dy[i] * m[i][1][1] + dz[i] * m[i][1][2]);
			u->dz[i] = (dx[i] * m[i][2][0] + dy[i] * m[i][2][1] + dz[i] * m[i][2][2]);
		}
		delete [] m;
	}

	if (shfn->get_type() == Hcurl) {
		double *nx, *ny, *nz;
		rm->calc_face_normal(iface, np, pt, nx, ny, nz);

		u->fn0 = new double [np]; MEM_CHECK(u->fn0);
		u->fn1 = new double [np]; MEM_CHECK(u->fn1);
		u->fn2 = new double [np]; MEM_CHECK(u->fn2);

		double *fn[3];
		for (int c = 0; c < 3; c++)
			fn[c] = shfn->get_fn_values(c);

		double3x3 *m = rm->get_inv_ref_map(np, pt);
		for (int i = 0; i < np; i++) {
			double ev[3] = {
				fn[0][i] * m[i][0][0] + fn[1][i] * m[i][0][1] + fn[2][i] * m[i][0][2],
				fn[0][i] * m[i][1][0] + fn[1][i] * m[i][1][1] + fn[2][i] * m[i][1][2],
				fn[0][i] * m[i][2][0] + fn[1][i] * m[i][2][1] + fn[2][i] * m[i][2][2]
			};
			double tpe[3];
			calc_tan_proj(nx[i], ny[i], nz[i], ev, tpe);
			u->fn0[i] = tpe[0];
			u->fn1[i] = tpe[1];
			u->fn2[i] = tpe[2];
		}

		delete [] m;
		delete [] nx;
		delete [] ny;
		delete [] nz;
	}

	return u;
}

mfn_t *init_fn(MeshFunction *f, RefMap *rm, const int np, const QuadPt3D *pt) {
	_F_

	mfn_t *u = new mfn_t;
	u->nc = f->get_num_components();
	f->precalculate(np, pt, FN_DEFAULT);
	if (u->nc == 1) {
		u->fn = new scalar [np]; MEM_CHECK(u->fn);
		u->dx = new scalar [np]; MEM_CHECK(u->dx);
		u->dy = new scalar [np]; MEM_CHECK(u->dy);
		u->dz = new scalar [np]; MEM_CHECK(u->dz);

		memcpy(u->fn, f->get_fn_values(), np * sizeof(scalar));
		memcpy(u->dx, f->get_dx_values(), np * sizeof(scalar));
		memcpy(u->dy, f->get_dy_values(), np * sizeof(scalar));
		memcpy(u->dz, f->get_dz_values(), np * sizeof(scalar));
	}
	else if (u->nc == 3) {
		// FN
		u->fn0 = new scalar [np]; MEM_CHECK(u->fn0);
		u->fn1 = new scalar [np]; MEM_CHECK(u->fn1);
		u->fn2 = new scalar [np]; MEM_CHECK(u->fn2);

		memcpy(u->fn0, f->get_fn_values(0), np * sizeof(scalar));
		memcpy(u->fn1, f->get_fn_values(1), np * sizeof(scalar));
		memcpy(u->fn2, f->get_fn_values(2), np * sizeof(scalar));

		// DX
		u->dx0 = new scalar [np]; MEM_CHECK(u->dx0);
		u->dx1 = new scalar [np]; MEM_CHECK(u->dx1);
		u->dx2 = new scalar [np]; MEM_CHECK(u->dx2);

		memcpy(u->dx0, f->get_dx_values(0), np * sizeof(scalar));
		memcpy(u->dx1, f->get_dx_values(1), np * sizeof(scalar));
		memcpy(u->dx2, f->get_dx_values(2), np * sizeof(scalar));

		// DY
		u->dy0 = new scalar [np]; MEM_CHECK(u->dy0);
		u->dy1 = new scalar [np]; MEM_CHECK(u->dy1);
		u->dy2 = new scalar [np]; MEM_CHECK(u->dy2);

		memcpy(u->dy0, f->get_dy_values(0), np * sizeof(scalar));
		memcpy(u->dy1, f->get_dy_values(1), np * sizeof(scalar));
		memcpy(u->dy2, f->get_dy_values(2), np * sizeof(scalar));

		// DZ
		u->dz0 = new scalar [np]; MEM_CHECK(u->dz0);
		u->dz1 = new scalar [np]; MEM_CHECK(u->dz1);
		u->dz2 = new scalar [np]; MEM_CHECK(u->dz2);

		memcpy(u->dz0, f->get_dz_values(0), np * sizeof(scalar));
		memcpy(u->dz1, f->get_dz_values(1), np * sizeof(scalar));
		memcpy(u->dz2, f->get_dz_values(2), np * sizeof(scalar));
	}
	else {
		EXIT(H3D_ERR_NOT_IMPLEMENTED);
	}

	return u;
}

void free_fn(fn_t<ord_t> *f) {
	delete f->fn;
}

template<typename T>
void free_fn_tpl(fn_t<T> *f) {
	delete [] f->fn;
	delete [] f->dx;
	delete [] f->dy;
	delete [] f->dz;

	delete [] f->fn0; delete [] f->fn1; delete [] f->fn2;
	delete [] f->dx0; delete [] f->dx1; delete [] f->dx2;
	delete [] f->dy0; delete [] f->dy1; delete [] f->dy2;
	delete [] f->dz0; delete [] f->dz1; delete [] f->dz2;
	delete [] f->curl0; delete [] f->curl1; delete [] f->curl2;

	delete f;
}

void free_fn(sfn_t *f) { free_fn_tpl<double>(f); }
#ifdef H3D_COMPLEX
void free_fn(mfn_t *f) { free_fn_tpl<scalar>(f); }
#endif
