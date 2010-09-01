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

#ifndef _DETERMINANT_H_
#define _DETERMINANT_H_

inline double det2(
	double a11, double a12,
	double a21, double a22)
{
	return a11*a22-a12*a21;
}

inline double det3(
	double a11, double a12, double a13,
	double a21, double a22, double a23,
	double a31, double a32, double a33)
{
	return a11*det2(a22, a23, a32, a33)
	      -a12*det2(a21, a23, a31, a33)
	      +a13*det2(a21, a22, a31, a32);
}

inline double det4(
	double a11, double a12, double a13, double a14,
	double a21, double a22, double a23, double a24,
	double a31, double a32, double a33, double a34,
	double a41, double a42, double a43, double a44)
{
	return a11*det3(a22, a23, a24, a32, a33, a34, a42, a43, a44)
	      -a12*det3(a21, a23, a24, a31, a33, a34, a41, a43, a44)
	      +a13*det3(a21, a22, a24, a31, a32, a34, a41, a42, a44)
	      -a14*det3(a21, a22, a23, a31, a32, a33, a41, a42, a43);
}

inline double det(const double3x3 &m) {
	return
		m[0][0] * m[1][1] * m[2][2] + m[0][1] * m[1][2] * m[2][0] + m[0][2] * m[1][0] * m[2][1] -
		m[2][0] * m[1][1] * m[0][2] - m[2][1] * m[1][2] * m[0][0] - m[2][2] * m[1][0] * m[0][1];
}

#endif // _DETERMINANT_H_
