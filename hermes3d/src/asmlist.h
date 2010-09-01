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

#ifndef _ASMLIST_H_
#define _ASMLIST_H_

#include <common/array.h>
#include <common/error.h>

/// AsmList is a simple container for the element assembly arrays idx, dof, coef. and ori
/// These arrays are filled by Space::get_element_assembly_list() and used by the
/// assembling procedure and the Solution class. The arrays are allocated and deallocated
/// automatically by the class.
///
/// @ingroup assembling
class AsmList {
public:
	long int *idx;
	int *dof;
	scalar *coef;

	int cnt, cap;

	AsmList() {
		idx = NULL;
		dof = NULL;
		coef = NULL;
		cnt = cap = 0;
	}

	~AsmList() {
    	free(idx);
    	free(dof);
    	free(coef);
	}

	void clear() {
		cnt = 0;
	}

	inline void add(long int idx, int dof, scalar coef) {
		if (coef == 0.0) return;
		if (cnt >= cap) enlarge();

		this->idx[cnt] = idx;
		this->dof[cnt] = dof;
		this->coef[cnt] = coef;
		cnt++;
	}

	void dump(FILE *stream = stdout) {
		fprintf(stream, "\nasmlist:\n");
		for (int i = 0; i < cnt; i++)
			fprintf(stream, " [% 3d] ind = %ld, dof = %d, coef = " SCALAR_FMT "\n", i, idx[i], dof[i], SCALAR(coef[i]));
		fprintf(stream, "\n");
	}

protected:
	void enlarge() {
		cap = !cap ? 256 : cap * 2;

		idx = (long int *) realloc(idx, sizeof(long int) * cap); MEM_CHECK(idx);
		dof = (int *) realloc(dof, sizeof(int) * cap); MEM_CHECK(dof);
		coef = (scalar *) realloc(coef, sizeof(scalar) * cap); MEM_CHECK(coef);
	}
};


#endif
