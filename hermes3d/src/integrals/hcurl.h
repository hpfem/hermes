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

#ifndef _INTEGRALS_HCURL_H_
#define _INTEGRALS_HCURL_H_

#include "refmap.h"
#include <common/trace.h>
#include <common/error.h>
#include <common/callstack.h>

/// Calculate curl of a vector
///
/// @param dx[in] - vector of x-derivatives
/// @param dy[in] - vector of y-derivatives
/// @param dz[in] - vector of z-derivatives
/// @param curl[out] - will contain the calculated curl
template<typename T>
void calc_curl(T (&dx)[3], T (&dy)[3], T (&dz)[3], T (&curl)[3]) {
	curl[0] = dy[2] - dz[1];
	curl[1] = dz[0] - dx[2];
	curl[2] = dx[1] - dy[0];
}

/// calculates N x EV x N
///
/// @param nx[in] - x-component of the normal
/// @param ny[in] - y-component of the normal
/// @param nz[in] - z-component of the normal
/// @param ev[in] - vector to work with
/// @param tpe[out] - will contain the calculated value
template<typename S, typename T>
void calc_tan_proj(S nx, S ny, S nz, T (&ev)[3], T (&tpe)[3]) {
	T ne[3] = {
		ny * ev[2] - nz * ev[1],
		nz * ev[0] - nx * ev[2],
		nx * ev[1] - ny * ev[0]
	};

	tpe[0] = nz * ne[1] - ny * ne[2];
	tpe[1] = nx * ne[2] - nz * ne[0];
	tpe[2] = ny * ne[0] - nx * ne[1];
}

/// @ingroup hcurlintegrals

template<typename f_t, typename res_t>
res_t hcurl_int_u_v(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->fn0[i] * v->fn0[i] + u->fn1[i] * v->fn1[i] + u->fn2[i] * v->fn2[i]);
	return result;
}

template<typename f_t, typename res_t>
res_t hcurl_int_curl_u_curl_v(int n, double *wt, fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->curl0[i] * v->curl0[i] + u->curl1[i] * v->curl1[i] + u->curl2[i] * v->curl2[i]);
	return result;
}

/// Integral \F \v
///
template<typename f_t, typename res_t>
res_t hcurl_int_F_v(int n, double *wt, void (*F)(f_t, f_t, f_t, res_t (&)[3]), fn_t<f_t> *v, geom_t<f_t> *e) {
	_F_
	res_t result = 0;
	res_t f[3];
	for (int i = 0; i < n; i++) {
		F(e->x[i], e->y[i], e->z[i], f);
		result += wt[i] * (v->fn0[i] * f[0] + v->fn1[i] * f[1] + v->fn2[i] * f[2]);
	}
	return result;
}

#endif
