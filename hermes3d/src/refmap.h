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

#ifndef _REFMAP_H_
#define _REFMAP_H_

#include "common.h"
#include "transform.h"
#include "quad.h"
#include "mesh.h"
#include "shapefn.h"

/// @defgroup refmap Reference mapping
///
/// TODO: description


/// Reference mapping (for evaluating integrals on a physical domain)
///
/// @ingroup refmap
class RefMap : public Transformable {
public:
	/// NOTE: Use this constructor only if constructing an array of RefMaps, then call set_mesh immediately!
	RefMap();
	/// NOTE: Always use this constructor whenever it is possible
	RefMap(Mesh *mesh);
	virtual ~RefMap();

	/// Sets the mesh where the reference map is working
	/// @param mesh [in] Pointer to the mesh.
	void set_mesh(Mesh *mesh) { this->mesh = mesh; }

	/// Initializes the reference map for the specified element.
	/// Must be called prior to using all other functions in the class.
	/// @param[in] e - The element we want to work with
	virtual void set_active_element(Element *e);

	/// @return The increase in the integration order due to the reference map.
	order3_t get_ref_order() const { return ref_order; }

	/// @return The increase in the integration order due to the inverse reference map.
	order3_t get_inv_ref_order() const { return inv_ref_order; }

	/// @return The array of jacobians of the reference map precalculated at the points 'pt'.
	/// @param[in] np - The number of points
	/// @param[in] pt - Points where the jacobian should be calculated
	/// @param[in] trans - set to true if you want transformed values
	double *get_jacobian(const int np, const QuadPt3D *pt, bool trans = true);

	/// @return The jacobi matrices of the reference map precalculated at the points 'pt'.
	/// @param[in] np - The number of points
	/// @param[in] pt - Points for which we want the jacobi matrices
	double3x3 *get_ref_map(const int np, const QuadPt3D *pt);

	/// @return The transposed inverse matrices of the reference map precalculated at the points 'pt'.
	/// @param[in] np - The number of points
	/// @param[in] pt - Points for which we want the transposed jacobi matrices
	double3x3 *get_inv_ref_map(const int np, const QuadPt3D *pt);

	/// @return The x-coordinates of the points transformed to the physical domain of the element
	/// @param[in] np - The number of points
	/// @param[in] pt - Points for which we want the x-coord
	double *get_phys_x(const int np, const QuadPt3D *pt);

	/// @return The y-coordinates of the points transformed to the physical domain of the element
	/// @param[in] np - The number of points
	/// @param[in] pt - Points for which we want the y-coord
	double *get_phys_y(const int np, const QuadPt3D *pt);

	/// @return The z-coordinates of the points transformed to the physical domain of the element
	/// @param[in] np - The number of points
	/// @param[in] pt - Points for which we want the z-coord
	double *get_phys_z(const int np, const QuadPt3D *pt);

	/// @return The array of 'face jacobians' at points 'pt'
	/// @param[in] trans - set to true if you want transformed values
	double *get_face_jacobian(int face, const int np, const QuadPt3D *pt, bool trans = true);

	/// Calculate normals on the face
	void calc_face_normal(int iface, const int np, const QuadPt3D *pt, double *&nx, double *&ny, double *&nz);

	/// See Transformable::push_transform()
	virtual void push_transform(int son);

	/// See Transformable::pop_transform()
	virtual void pop_transform();

	void force_transform(uint64 sub_idx, Trf *ctm);

protected:
	Mesh *mesh;
	ShapeFunction *pss;

	bool      is_const_jacobian;
	double    const_jacobian;
	double3x3 const_inv_ref_map;
	double3x3 const_ref_map;

	order3_t ref_order;
	order3_t inv_ref_order;

	int n_coefs;					// # of coeffs in 'indices' array
	int indices[70];				// FIXME: magic number

	Vertex *coefs;
	Vertex vertex[8];				// max number of vertices (hex has 8 vertices, other elements have less)

	void calc_const_inv_ref_map();
	double calc_face_const_jacobian(int face);

	friend class ExactSolution;
};

#endif
