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
#include "common.h"
#include "function.h"
#include <common/callstack.h>


// the order of items must match values of EValueType
template<typename TYPE>
int Function<TYPE>::idx2mask[][COMPONENTS] = {
	{ FN_VAL_0, FN_VAL_1, FN_VAL_2 },
	{ FN_DX_0,  FN_DX_1,  FN_DX_2  },
	{ FN_DY_0,  FN_DY_1,  FN_DY_2  },
	{ FN_DZ_0,  FN_DZ_1,  FN_DZ_2  },
	{ FN_DXX_0, FN_DXX_1, FN_DXX_2 },
	{ FN_DYY_0, FN_DYY_1, FN_DYY_2 },
	{ FN_DZZ_0, FN_DZZ_1, FN_DZZ_2 },
	{ FN_DXY_0, FN_DXY_1, FN_DXY_2 },
	{ FN_DYZ_0, FN_DYZ_1, FN_DYZ_2 },
	{ FN_DXZ_0, FN_DXZ_1, FN_DXZ_2 }
};


template<typename TYPE>
Function<TYPE>::Function() {
	_F_
	order = 0;
	max_mem = total_mem = 0;

	cur_node = NULL;
	memset(quads, 0, sizeof(quads));
	cur_quad = 0;
	sub_idx = 0;
}


template<typename TYPE>
Function<TYPE>::~Function() {
	_F_
}

template<typename TYPE>
typename Function<TYPE>::Node *Function<TYPE>::new_node(int mask, int num_points) {
	_F_
	// get the number of tables
	int nt = 0, m = mask;
	if (num_components < 3) m &= FN_VAL_0 | FN_DX_0 | FN_DY_0 | FN_DZ_0 | FN_DXX_0 | FN_DYY_0 | FN_DZZ_0 | FN_DXY_0 | FN_DXZ_0 | FN_DYZ_0;
	while (m) {
		nt += m & 1;
		m >>= 1;
	}

	// allocate a node including its data part, init table pointers
	int size = sizeof(Node) + sizeof(TYPE) * num_points * nt;
	Node *node = (Node *) malloc(size);
	node->mask = mask;
	node->size = size;
	memset(node->values, 0, sizeof(node->values));
	TYPE *data = node->data;
	for (int j = 0; j < num_components; j++) {
		for (int i = 0; i < VALUE_TYPES; i++)
			if (mask & idx2mask[i][j]) {
				node->values[j][i] = data;
				data += num_points;
			}
	}

	total_mem += size;
	if (max_mem < total_mem) max_mem = total_mem;

	return node;
}

#undef CHECK_PARAMS
#undef CHECK_TABLE
