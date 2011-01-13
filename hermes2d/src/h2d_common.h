// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 2 of the License,
// or (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_COMMON_H_
#define __H2D_COMMON_H_

#include "../../hermes_common/common.h"

#include "auto_local_array.h"
//#include "vector.h"


// H2D-specific error codes.
#define H2D_ERR_EDGE_INDEX_OUT_OF_RANGE         "Edge index out of range."

#define H2D_NUM_MODES 2 ///< A number of modes, see enum ElementMode2D.

// how many bits the order number takes
const int H2D_ORDER_BITS = 5;
const int H2D_ORDER_MASK = (1 << H2D_ORDER_BITS) - 1;

// macros for combining quad horizontal and vertical orders
#define H2D_MAKE_QUAD_ORDER(h_order, v_order) (((v_order) << H2D_ORDER_BITS) + (h_order))
#define H2D_GET_H_ORDER(order) ((order) & H2D_ORDER_MASK)
#define H2D_GET_V_ORDER(order) ((order) >> H2D_ORDER_BITS)

extern HERMES_API const std::string h2d_get_quad_order_str(const int quad_order); ///< Returns string representation of the quad order: used for debugging purposses.
extern HERMES_API int h2d_make_edge_order(int edge, int encoded_order, int mode); ///< Returns the correct axial order for given edge.

// Enabling second derivatives in weak forms. Turned off by default. Second
// derivatives are employed, among others, by stabilization methods for
// transport equations. For usage see the example linear-convection-diffusion.
#define H2D_SECOND_DERIVATIVES_ENABLED

/* Uncomment this line to disable internal mesh compatibility 
   tests in Traverse:begin(). */ 
//#define H2D_DISABLE_MULTIMESH_TESTS

#endif

