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

#ifndef _SHAPESET_H_
#define _SHAPESET_H_

#include "common.h"
#include "function.h"

/// @defgroup shapesets Shapesets
///
/// TODO: description

enum EPartOri {
	PART_ORI_HORZ = 0,
	PART_ORI_VERT = 1
};

struct Part {
	// quad part
	union {
		struct {					// face variant
			unsigned horz:16;		// horizontal part (stripe)
			unsigned vert:16;		// vertical part (stripe)
		};
		struct {					// edge variant
			unsigned part:32;		// (stripe)
		};
	};
};

void get_interval_part(int part, double &lo, double &hi);
void get_edge_part(int part, double &x);
double get_edge_coef(int part);

Part transform_edge_part(int ori, Part part);
Part transform_face_part(int ori, Part part);

inline int get_lower_part(int part) { 	return (part * 2) + 1; }
inline int get_higher_part(int part) { return (part + 1) * 2; }
inline int face_to_edge_part(int part) { return part + 2; }

int combine_face_part(int part, int finer_part);
int opposite_part(int part);
inline int opposite_edge_part(int part) { return opposite_part(part - 2) + 2; }

enum EShapeFnType {
	SHFN_CONSTRAINED = -1,
	SHFN_VERTEX = 0,
	SHFN_EDGE = 1,
	SHFN_FACE = 2,
	SHFN_BUBBLE = 3
};

enum ECedKeyType {
	CED_KEY_TYPE_EDGE = 0,
	CED_KEY_TYPE_FACE = 1,
	CED_KEY_TYPE_EDGE_FACE = 2
};

/// constrained key
struct CEDKey {
	unsigned type:2;			// see enum ECedKeyType
	unsigned ori:4;				// orientation of a edge/face
	union {
		unsigned face:4;		// index of a face
		unsigned edge:4;		// index of an edge
	};
	unsigned dir:1;
	unsigned variant:1;			// variant of the function (in Hcurl we need to distinguish between two variants of face functions)
	int order;					// order of an edge/face function
	Part part;					// part of the edge/face

	CEDKey() {
	}

	CEDKey(unsigned type, unsigned edge, order1_t order, unsigned ori, Part part) {
		this->type = type;
		this->ori = ori;
		this->edge = edge;
		this->order = order;
		this->part = part;
		this->dir = dir;
	}

	CEDKey(unsigned type, unsigned num, order2_t order, unsigned ori, Part part, int dir = 0, int variant = 0) {
		this->type = type;
		this->ori = ori;
		this->face = num;
		this->order = order.get_idx();
		this->part = part;
		this->dir = dir;
		this->variant = variant;
	}
};

///
struct CEDComb {
	int n;								// number of coefficients
	double *coef;						// coeficients

	CEDComb(int n, double *coefs) {
		this->coef = coefs;
		this->n = n;
	}

	~CEDComb() {
		delete [] coef;
	}
};


/// Base class for all shapesets
///
/// @ingroup shapesets
class Shapeset {
public:
	Shapeset(int id);
	virtual ~Shapeset();

	ESpaceType get_type() { return type; }
	int get_mode() const { return mode; }
	int get_num_components() const { return num_components; }

	// @return index of a vertex shape function for a vertex
	// @param [in] vertex - index of the vertex
	virtual int get_vertex_index(int vertex) const = 0;

	/// @return indices of edge shape functions
	/// @param [in] edge - edge number (local)
	/// @param [in] ori - orientation of the edge (0 or 1)
	/// @param [in] order - order on the edge
	virtual int *get_edge_indices(int edge, int ori, order1_t order) = 0;

	/// @return indices of face shape functions
	/// @param [in] face - face number (local)
	/// @param [in] ori - orientation of the face
	/// @param [in] order - order on the face
	virtual int *get_face_indices(int face, int ori, order2_t order) = 0;

	/// @return indices of bubble functions
	/// @param order - order of the bubble function
	virtual int *get_bubble_indices(order3_t order) = 0;

	virtual int get_num_edge_fns(order1_t order) const = 0;

	virtual int get_num_face_fns(order2_t order) const = 0;

	virtual int get_num_bubble_fns(order3_t order) const = 0;

	/// Get the number of possible orientations on a face
	///
	/// @return The number of possbile orientations on a face
	/// @param[in] face The number of the face for which the orientations are evaluated
	virtual int get_face_orientations(int face) const = 0;

	/// Get the number of possible orientations on an edge
	///
	/// @return The number of possbile orientations on an edge
	virtual int get_edge_orientations() const = 0;

	/// Returns the variant of a face function
	/// It is needed only in Hcurl-derived shapesets where we need to distinguish between to variants
	/// of shape functions, other shapesets are returning 0
	virtual int get_face_fn_variant(int index) const { return 0; }

	virtual order3_t get_order(int index) const = 0;

	/// Get function decomposition for product shapesets
	virtual order3_t get_dcmp(int index) const = 0;

	/// Get index of a constrained edge function.
	/// @return The index of a constrained edge function.
	/// @param[in] edge The local number of an edge.
	/// @param[in] order The polynomial order on the edge.
	/// @param[in] ori The orientation of the edge function.
	/// @param[in] part The 'part' of an edge
	virtual int get_constrained_edge_index(int edge, int ori, order1_t order, Part part);

	/// Get index of a edge function constrained by a face function.
	/// @return The index of a constrained edge function.
	/// @param[in] edge The local number of an edge.
	/// @param[in] order The polynomial order on the edge.
	/// @param[in] ori The orientation of the edge function.
	/// @param[in] part The 'part' of an edge
	virtual int get_constrained_edge_face_index(int edge, int ori, order2_t order, Part part, int dir, int variant = 0);

	/// Get index of a constrained face function.
	/// @return The index of a constrained face function.
	/// @param[in] face The local number of a face.
	/// @param[in] order The polynomial order on the face.
	/// @param[in] ori The orientation of the face function.
	/// @param[in] part The 'part' of an face
	virtual int get_constrained_face_index(int face, int ori, order2_t order, Part part, int variant = 0);

	virtual int get_shape_type(int index) const = 0;

	/// Evaluate function in the set of points
	/// @param[in] n
	/// @param[in] index - Index of the function being evaluate
	/// @param[in] np - The number of points in array 'pt'
	/// @param[in] pt - Points where the function is evaluated
	/// @param[in] component - The number of component of the evaluated function
	/// @param[out] vals - The array of vakues (caller is responsible for freeing this memory)
	virtual void get_values(int n, int index, int np, QuadPt3D *pt, int component, double *vals) = 0;

	/// Evaluate function 'index' in points 'pt'
	virtual double get_value(int n, int index, double x, double y, double z, int component) = 0;

	inline void get_fn_values (int index, int np, QuadPt3D *pt, int component, double *vals) { get_values(FN,  index, np, pt, component, vals); }
	inline void get_dx_values (int index, int np, QuadPt3D *pt, int component, double *vals) { get_values(DX,  index, np, pt, component, vals); }
	inline void get_dy_values (int index, int np, QuadPt3D *pt, int component, double *vals) { get_values(DY,  index, np, pt, component, vals); }
	inline void get_dz_values (int index, int np, QuadPt3D *pt, int component, double *vals) { get_values(DZ,  index, np, pt, component, vals); }
	inline void get_dxx_values(int index, int np, QuadPt3D *pt, int component, double *vals) { get_values(DXX, index, np, pt, component, vals); }
	inline void get_dyy_values(int index, int np, QuadPt3D *pt, int component, double *vals) { get_values(DYY, index, np, pt, component, vals); }
	inline void get_dzz_values(int index, int np, QuadPt3D *pt, int component, double *vals) { get_values(DZZ, index, np, pt, component, vals); }
	inline void get_dxy_values(int index, int np, QuadPt3D *pt, int component, double *vals) { get_values(DXY, index, np, pt, component, vals); }
	inline void get_dxz_values(int index, int np, QuadPt3D *pt, int component, double *vals) { get_values(DXZ, index, np, pt, component, vals); }
	inline void get_dyz_values(int index, int np, QuadPt3D *pt, int component, double *vals) { get_values(DYZ, index, np, pt, component, vals); }

	inline double get_fn_value (int index, double x, double y, double z, int component) { return get_value(FN,  index, x, y, z, component); }
	inline double get_dx_value (int index, double x, double y, double z, int component) { return get_value(DX,  index, x, y, z, component); }
	inline double get_dy_value (int index, double x, double y, double z, int component) { return get_value(DY,  index, x, y, z, component); }
	inline double get_dz_value (int index, double x, double y, double z, int component) { return get_value(DZ,  index, x, y, z, component); }
	inline double get_dxx_value(int index, double x, double y, double z, int component) { return get_value(DXX, index, x, y, z, component); }
	inline double get_dyy_value(int index, double x, double y, double z, int component) { return get_value(DYY, index, x, y, z, component); }
	inline double get_dzz_value(int index, double x, double y, double z, int component) { return get_value(DZZ, index, x, y, z, component); }
	inline double get_dxy_value(int index, double x, double y, double z, int component) { return get_value(DXZ, index, x, y, z, component); }
	inline double get_dxz_value(int index, double x, double y, double z, int component) { return get_value(DXZ, index, x, y, z, component); }
	inline double get_dyz_value(int index, double x, double y, double z, int component) { return get_value(DYZ, index, x, y, z, component); }

public:
	const int id;
protected:
	ESpaceType type;
	int mode;
	int num_components;

	// CED
	void   get_constrained_values(int n, int index, int np, QuadPt3D *pt, int component, double *vals);
	double get_constrained_value (int n, int index, double x, double y, double z, int component);

	order3_t get_ced_order(int index) const;
	virtual CEDComb *calc_constrained_edge_combination(int ori, const order1_t &order, Part part) { return NULL; }
	virtual CEDComb *calc_constrained_edge_face_combination(int ori, const order2_t &order, Part part, int dir, int variant = 0) { return NULL; }
	virtual CEDComb *calc_constrained_face_combination(int ori, const order2_t &order, Part part, int variant = 0) { return NULL; }
	void free_constrained_combinations();

	Map<CEDKey, CEDComb *> ced_comb;			// mapping: CEDKey => CEDComb
	Map<CEDKey, int> ced_id;					// mapping: CEDKey => ced function index
	Array<CEDKey> ced_key;						// indexing: index => CEDKey
	int ced_idx;								// ced index to assing

	CEDComb *get_ced_comb(const CEDKey &key);
	int *get_ced_indices(const CEDKey &key);

#ifdef PRELOADING
public:
	bool load_prods(const char *file_name, double *&mat);
	bool preload_products();

	scalar get_product_val(int idx1, int idx2, double *vals);

	inline scalar get_fn_product(int idx1, int idx2) { return get_product_val(idx1, idx2, fn_prods); }
	inline scalar get_dx_product(int idx1, int idx2) { return get_product_val(idx1, idx2, dx_prods); }
	inline scalar get_dy_product(int idx1, int idx2) { return get_product_val(idx1, idx2, dy_prods); }
	inline scalar get_dz_product(int idx1, int idx2) { return get_product_val(idx1, idx2, dz_prods); }

protected:
	Array<int> fnidx2idx;				// mapping from fn index to a array index
	double *fn_prods;					// products of fn. values, dx, dy, dz
	double *dx_prods;
	double *dy_prods;
	double *dz_prods;
	int num_fns;

#endif
};


#endif
