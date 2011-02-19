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

#ifndef _INTEGRALS_H1_H_
#define _INTEGRALS_H1_H_

#include "refmap.h"
#include "order.h"
#include "../../hermes_common/error.h"
#include "../../hermes_common/trace.h"
#include "../../hermes_common/callstack.h"

/// @defgroup h1integrals H1 integrals

/// Integral \u
///
/// @ingroup h1integrals
template<typename f_t, typename res_t>
res_t int_u(int n, double *wt, Func<f_t> *u, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->val[i]);
	return result;
}

/// Integral \v
///
/// @ingroup h1integrals
template<typename f_t, typename res_t>
res_t int_v(int n, double *wt, Func<f_t> *v, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (v->val[i]);
	return result;
}


/// Integral \u \v
///
/// @ingroup h1integrals
template<typename f_t, typename res_t>
res_t int_u_v(int n, double *wt, Func<f_t> *u, Func<f_t> *v, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->val[i] * v->val[i]);
	return result;
}

/// Integral \F \u
///
/// @ingroup h1integrals
template<typename f_t, typename res_t>
res_t int_F_v(int n, double *wt, res_t (*F)(f_t x, f_t y, f_t z), Func<f_t> *v, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (v->val[i] * (*F)(e->x[i], e->y[i], e->z[i]));
	return result;
}

#define int_F_v(n, wt, F, v, e) int_F_u(n, wt, F, v, e)

/// Integral \grad u \grad v
///
/// @ingroup h1integrals
template<typename f_t, typename res_t>
res_t int_grad_u_grad_v(int n, double *wt, Func<f_t> *u, Func<f_t> *v, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] + u->dz[i] * v->dz[i]);
	return result;
}

/// Integral u u \dd{v}{x}
///
/// @ingroup h1integrals
template<typename f_t, typename res_t>
res_t int_u_dvdx(int n, double *wt, Func<f_t> *u, Func<f_t> *v, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->val[i] * v->dx[i]);
	return result;
}

/// Integral u u \dd{v}{y}
///
/// @ingroup h1integrals
template<typename f_t, typename res_t>
res_t int_u_dvdy(int n, double *wt, Func<f_t> *u, Func<f_t> *v, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->val[i] * v->dy[i]);
	return result;
}

/// Integral u \dd{v}{z}
///
/// @ingroup h1integrals
template<typename f_t, typename res_t>
res_t int_u_dvdz(int n, double *wt, Func<f_t> *u, Func<f_t> *v, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->val[i] * v->dz[i]);
	return result;
}

/// Integral u u \dd{v}{x}
///
/// @ingroup h1integrals
template<typename f_t, typename res_t>
res_t int_dudx_v(int n, double *wt, Func<f_t> *u, Func<f_t> *v, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->dx[i] * v->val[i]);
	return result;
}

/// Integral u u \dd{v}{y}
///
/// @ingroup h1integrals
template<typename f_t, typename res_t>
res_t int_dudy_v(int n, double *wt, Func<f_t> *u, Func<f_t> *v, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->dy[i] * v->val[i]);
	return result;
}

/// Integral u \dd{v}{z}
///
/// @ingroup h1integrals
template<typename f_t, typename res_t>
res_t int_dudz_v(int n, double *wt, Func<f_t> *u, Func<f_t> *v, Geom<f_t> *e) {
	_F_
	res_t result = 0;
	for (int i = 0; i < n; i++)
		result += wt[i] * (u->dz[i] * v->val[i]);
	return result;
}

#endif
