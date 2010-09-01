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

#ifndef _FUNCTION_H_
#define _FUNCTION_H_

#include "common.h"
#include "mesh.h"
#include "quad.h"
#include "order.h"
#include "transform.h"
#include <common/error.h>

// Plenty of checking stuff for the debug version
#ifdef DEBUG

#define CHECK_PARAMS \
	if (component < 0 || component > num_components) \
		EXIT("Invalid component. You are probably using scalar-valued shapeset for an Hcurl problem."); \
	if (cur_node == NULL) EXIT("Invalid node. Did you call precalculate()?");

#define CHECK_TABLE(n, msg) \
	if (cur_node->values[component][n] == NULL) \
		EXIT(msg " not precalculated for component %d. Did you call precalculate() with correct mask?", component);

#else
  #define CHECK_PARAMS
  #define CHECK_TABLE(n, msg)
#endif

/// 3 components (we are in 3D)
static const int COMPONENTS = 3;

/// number of types of values in EValueType
static const int VALUE_TYPES = 10;

enum EValueType {
	FN = 0,
	DX = 1,
	DY = 2,
	DZ = 3,
	DXX = 4,
	DYY = 5,
	DZZ = 6,
	DXY = 7,
	DYZ = 8,
	DXZ = 9
};

// Precalculation masks
enum {
	FN_VAL_0  = 0x00000001, /// Function values of the 1st component required
	FN_DX_0   = 0x00000002, /// First derivative in x of the 1st component required
	FN_DY_0   = 0x00000004, /// First derivative in y of the 1st component required
	FN_DZ_0   = 0x00000008, /// First derivative in z of the 1st component required

	FN_DXX_0  = 0x00000010, /// Second derivative in x of the 1st component required
	FN_DYY_0  = 0x00000020, /// Second derivative in y of the 1st component required
	FN_DZZ_0  = 0x00000040, /// Second derivative in y of the 1st component required

	FN_DXY_0  = 0x00000080, /// Second mixed derivative of the 1st component required
	FN_DYZ_0  = 0x00000100, /// Second mixed derivative of the 1st component required
	FN_DXZ_0  = 0x00000200, /// Second mixed derivative of the 1st component required

	FN_VAL_1  = 0x00000400, /// Function values of the 2nd component required
	FN_DX_1   = 0x00000800, /// First derivative in x of the 2nd component required
	FN_DY_1   = 0x00001000, /// First derivative in y of the 2nd component required
	FN_DZ_1   = 0x00002000, /// First derivative in y of the 2nd component required

	FN_DXX_1  = 0x00004000, /// Second derivative in x of the 2nd component required
	FN_DYY_1  = 0x00008000, /// Second derivative in y of the 2nd component required
	FN_DZZ_1  = 0x00010000, /// Second derivative in y of the 2nd component required

	FN_DXY_1  = 0x00020000, /// Second mixed derivative of the 2nd component required
	FN_DYZ_1  = 0x00040000, /// Second mixed derivative of the 2nd component required
	FN_DXZ_1  = 0x00080000, /// Second mixed derivative of the 2nd component required

	FN_VAL_2  = 0x00100000, /// Function values of the 3rd component required
	FN_DX_2   = 0x00200000, /// First derivative in x of the 3rd component required
	FN_DY_2   = 0x00400000, /// First derivative in y of the 3rd component required
	FN_DZ_2   = 0x00800000, /// First derivative in y of the 3rd component required

	FN_DXX_2  = 0x01000000, /// Second derivative in x of the 3rd component required
	FN_DYY_2  = 0x02000000, /// Second derivative in y of the 3rd component required
	FN_DZZ_2  = 0x04000000, /// Second derivative in y of the 3rd component required

	FN_DXY_2  = 0x08000000, /// Second mixed derivative of the 3rd component required
	FN_DYZ_2  = 0x10000000, /// Second mixed derivative of the 3rd component required
	FN_DXZ_2  = 0x20000000, /// Second mixed derivative of the 3rd component required
};

// All components are usually requested together...
const int FN_VAL = FN_VAL_0 | FN_VAL_1 | FN_VAL_2;
const int FN_DX  = FN_DX_0  | FN_DX_1  | FN_DX_2;
const int FN_DY  = FN_DY_0  | FN_DY_1  | FN_DY_2;
const int FN_DZ  = FN_DZ_0  | FN_DZ_1  | FN_DZ_2;
const int FN_DXX = FN_DXX_0 | FN_DXX_1 | FN_DXX_2;
const int FN_DYY = FN_DYY_0 | FN_DYY_1 | FN_DYY_2;
const int FN_DZZ = FN_DZZ_0 | FN_DZZ_1 | FN_DZZ_2;

const int FN_DXY = FN_DXY_0 | FN_DXY_1 | FN_DXY_2;
const int FN_DYZ = FN_DYZ_0 | FN_DYZ_1 | FN_DYZ_2;
const int FN_DXZ = FN_DXZ_0 | FN_DXZ_1 | FN_DXZ_2;

// The most common case is FN_DEFAULT, ie. values and the gradient.
const int FN_DEFAULT = FN_VAL | FN_DX | FN_DY | FN_DZ;            /// Default precalculation mask
/// Precalculate everything
const int FN_ALL = FN_DEFAULT | FN_DXX | FN_DYY | FN_DZZ | FN_DXY | FN_DYZ | FN_DXZ;

// First, second, third component
const int FN_COMPONENT_0 = FN_VAL_0 | FN_DX_0 | FN_DY_0 | FN_DZ_0 | FN_DXX_0 | FN_DYY_0 | FN_DZZ_0 | FN_DXY_0 | FN_DYZ_0 | FN_DXZ_0;
const int FN_COMPONENT_1 = FN_VAL_1 | FN_DX_1 | FN_DY_1 | FN_DZ_1 | FN_DXX_1 | FN_DYY_1 | FN_DZZ_1 | FN_DXY_1 | FN_DYZ_1 | FN_DXZ_1;
const int FN_COMPONENT_2 = FN_VAL_2 | FN_DX_2 | FN_DY_2 | FN_DZ_2 | FN_DXX_2 | FN_DYY_2 | FN_DZZ_2 | FN_DXY_2 | FN_DYZ_2 | FN_DXZ_2;

inline void mask_to_comp_val(int mask, int &a, int &b) {
	a = 0, b = 0;

	if (mask >= 0x00100000) {			// where the component 1 changes to component 2
		a = 2;
		mask >>= 20;					// number of bits used for comp[0] and comp[1]
	}
	else if (mask >= 0x00000400) {		// where the component 0 changes to component 1
		a = 1;
		mask >>= 10;					// number of bits used for comp[0]
	}
	while (!(mask & 1)) {
		mask >>= 1;
		b++;
	}
}


template<typename TYPE>
class Function : public Transformable {
public:
	/// Default constructor.
	Function();
	/// Dectructor.
	virtual ~Function();

	/// @return The polynomial degree of the function currently being represented by the class.
	order3_t get_fn_order() const { return order; }

	/// @return The number of vector components of the function being represented by the class.
	int get_num_components() const { return num_components; }

	/// @param component [in] The component of the function (0-2).
	/// @return The values of the function at all points of the current integration rule.
	TYPE *get_fn_values(int component = 0) {
		CHECK_PARAMS;
		CHECK_TABLE(FN, "Function values");
		return cur_node->values[component][FN];
	}

	/// @param component [in] The component of the function (0-2).
	/// @return The x partial derivative of the function at all points of the current integration rule.
	TYPE *get_dx_values(int component = 0) {
		CHECK_PARAMS;
		CHECK_TABLE(DX, "DX values");
		return cur_node->values[component][DX];
	}

	/// @param component [in] The component of the function (0-2).
	/// @return The y partial derivative of the function at all points of the current integration rule.
	TYPE *get_dy_values(int component = 0) {
		CHECK_PARAMS;
		CHECK_TABLE(DY, "DY values");
		return cur_node->values[component][DY];
	}

	/// @param component [in] The component of the function (0-2).
	/// @return The z partial derivative of the function at all points of the current integration rule.
	TYPE *get_dz_values(int component = 0) {
		CHECK_PARAMS;
		CHECK_TABLE(DZ, "DZ values");
		return cur_node->values[component][DZ];
	}

	/// This function provides the both often-used dx and dy values in one call.
	/// @param dx [out] Variable which receives the pointer to the first partial derivatives by x
	/// @param dy [out] Variable which receives the pointer to the first partial derivatives by y
	/// @param dz [out] Variable which receives the pointer to the first partial derivatives by z
	/// @param component [in] The component of the function (0 or 1).
	void get_dx_dy_dz_values(TYPE *&dx, TYPE *&dy, TYPE *&dz, int component = 0) {
		CHECK_PARAMS;
		CHECK_TABLE(DX, "DX values");
		CHECK_TABLE(DY, "DY values");
		CHECK_TABLE(DZ, "DZ values");

		dx = cur_node->values[component][DX];
		dy = cur_node->values[component][DY];
		dz = cur_node->values[component][DZ];
	}

	/// @param component [in] The component of the function (0 or 1).
	/// @return The x second partial derivative of the function at all points of the current integration rule.
	TYPE *get_dxx_values(int component = 0) {
		CHECK_PARAMS;
		CHECK_TABLE(DXX, "DXX values");
		return cur_node->values[component][DXX];
	}

	/// @param component [in] The component of the function (0 or 1).
	/// @return The y second partial derivative of the function at all points of the current integration rule.
	TYPE *get_dyy_values(int component = 0) {
		CHECK_PARAMS;
		CHECK_TABLE(DYY, "DYY values");
		return cur_node->values[component][DYY];
	}

	/// @param component [in] The component of the function (0 or 1).
	/// @return The z second partial derivative of the function at all points of the current integration rule.
	TYPE *get_dzz_values(int component = 0) {
		CHECK_PARAMS;
		CHECK_TABLE(DZZ, "DZZ values");
		return cur_node->values[component][DZZ];
	}

	/// @param component [in] The component of the function (0 or 1).
	/// @return The second mixed derivative of the function at all points of the current integration rule.
	TYPE *get_dxy_values(int component = 0) {
		CHECK_PARAMS;
		CHECK_TABLE(DXY, "DXY values");
		return cur_node->values[component][DXY];
	}

	/// @param component [in] The component of the function (0 or 1).
	/// @return The second mixed derivative of the function at all points of the current integration rule.
	TYPE *get_dyz_values(int component = 0) {
		CHECK_PARAMS;
		CHECK_TABLE(DYZ, "DYZ values");
		return cur_node->values[component][DYZ];
	}

	/// @param component [in] The component of the function (0 or 1).
	/// @return The second mixed derivative of the function at all points of the current integration rule.
	TYPE *get_dxz_values(int component = 0) {
		CHECK_PARAMS;
		CHECK_TABLE(DXZ, "DXZ values");
		return cur_node->values[component][DXZ];
	}

	/// For internal use.
	TYPE *get_values(int component, int b) {
		CHECK_PARAMS;
		return cur_node->values[component][b];
	}

	/// Frees all precalculated tables.
	virtual void free() = 0;

	/// precalculates the current function at the current integration points.
	virtual void precalculate(const int np, const QuadPt3D *pt, int mask) = 0;

protected:
	static const int QUAD_COUNT = 8;

	order3_t order;			/// current function polynomial order
	int num_components; 	/// number of vector components

	struct Node {
		int mask;                               /// a combination of FN_XXX: specifies which tables are present
		int size;                               /// size in bytes of this struct (for maintaining total_mem)
		TYPE *values[COMPONENTS][VALUE_TYPES];  /// pointers to 'data'
		TYPE data[0];                           /// value tables
	};

	Node *cur_node;

	Quad3D *quads[QUAD_COUNT];	/// list of available quadratures
	int cur_quad;				/// active quadrature (index into 'quads')
	int total_mem;				/// total memory in bytes used by the tables
	int max_mem;				/// peak memory usage

	int find_quad(Quad3D *quad); /// searches 'quads' for the given quadrature
	int register_quad(Quad3D *quad); ///
	virtual Node *new_node(int mask, int num_points); /// allocates a new Node structure

	void free_cur_node() {
	    if (cur_node != NULL) {
	    	total_mem -= cur_node->size;
	    	::free(cur_node);
	    	cur_node = NULL;
	    }
	}

	void replace_cur_node(Node *node) {
		free_cur_node();
		cur_node = node;
	}

	/// For internal use only.
	void force_transform(uint64 sub_idx, Trf *ctm) {
		this->sub_idx = sub_idx;
		this->ctm = ctm;
		free_cur_node();
	}

	static int idx2mask[][COMPONENTS];  /// index to mask table (3 dimensions)
};


typedef
	Function<double> RealFunction;

typedef
	Function<scalar> ScalarFunction;


#endif
