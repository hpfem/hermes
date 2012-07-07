// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __H2D_SHAPESET_H1_QUAD_H
#define __H2D_SHAPESET_H1_QUAD_H

extern Shapeset::shape_fn_t* simple_quad_shape_fn_table[1];
extern Shapeset::shape_fn_t* simple_quad_shape_fn_table_dx[1];
extern Shapeset::shape_fn_t* simple_quad_shape_fn_table_dy[1];
extern Shapeset::shape_fn_t* simple_quad_shape_fn_table_dxx[1];
extern Shapeset::shape_fn_t* simple_quad_shape_fn_table_dxy[1];
extern Shapeset::shape_fn_t* simple_quad_shape_fn_table_dyy[1];

extern int simple_quad_vertex_indices[4];
extern int* simple_quad_edge_indices[4];
extern int* simple_quad_bubble_indices[];
extern int simple_quad_bubble_count[];
extern int simple_quad_index_to_order[];

#endif