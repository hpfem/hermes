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
#include "matrix.h"

#include <common/error.h>
#include <common/callstack.h>


#define H3D_TINY 1.0e-20


// ludcmp, lubksb - LU decomposition and back-substitution routines from
// the book Numerical Recipes in C, adjusted to zero-based indexing

void ludcmp(double **a, int n, int *indx, double *d)
{
	_F_
	int i, imax = 0, j, k;
	double big, dum, sum, temp;
	double *vv = new double[n];
	MEM_CHECK(vv);

	*d = 1.0;
	for (i = 0; i < n; i++) {
		big=0.0;
		for (j = 0; j < n; j++)
			if ((temp = fabs(a[i][j])) > big)
				big = temp;
		if (big == 0.0) EXIT("Singular matrix in routine LUDCMP!");
		vv[i] = 1.0 / big;
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
			if ((dum = vv[i]*fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) a[j][j] = H3D_TINY;
		if (j != n-1) {
			dum = 1.0 / (a[j][j]);
			for (i = j+1; i < n; i++) a[i][j] *= dum;
		}
	}
	delete [] vv;
}

// choldc, cholsl - Cholesky decomposition and solution routines from
// the book Numerical Recipes in C, adjusted to zero-based indexing

void choldc(double **a, int n, double p[])
{
	_F_
	int i, j, k;
	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			double sum = a[i][j];
			k = i;
			while (--k >= 0)
				sum -= a[i][k] * a[j][k];

			if (i == j) {
				if (sum <= 0.0) EXIT("CHOLDC failed!");
				else p[i] = sqrt(sum);
			}
			else
				a[j][i] = sum / p[i];
		}
	}
}

// SparseMatrix ////////////////////////////////////////////////////////////////////////////////////

void qsort_int(int* pbase, size_t total_elems); // defined in qsort.cpp

SparseMatrix::SparseMatrix()
{
	_F_
	size = 0;
	pages = NULL;

	row_storage = false;
	col_storage = false;
}

SparseMatrix::~SparseMatrix()
{
	_F_
	delete [] pages;
}

void SparseMatrix::prealloc(int n)
{
	_F_
	this->size = n;

	pages = new Page *[n];
	MEM_CHECK(pages);
	memset(pages, 0, n * sizeof(Page *));
}

void SparseMatrix::pre_add_ij(int row, int col)
{
	_F_
	if (pages[col] == NULL || pages[col]->count >= PAGE_SIZE) {
		Page *new_page = new Page;
		MEM_CHECK(new_page);
		new_page->count = 0;
		new_page->next = pages[col];
		pages[col] = new_page;
	}
	pages[col]->idx[pages[col]->count++] = row;
}

int SparseMatrix::sort_and_store_indices(Page *page, int *buffer, int *max)
{
	_F_
	// gather all pages in the buffer, deleting them along the way
	int *end = buffer;
	while (page != NULL) {
		memcpy(end, page->idx, sizeof(int) * page->count);
		end += page->count;
		Page *tmp = page;
		page = page->next;
		delete tmp;
	}

	// sort the indices and remove duplicities
	qsort_int(buffer, end - buffer);
	int *q = buffer;
	for (int *p = buffer, last = -1; p < end; p++)
		if (*p != last)
			*q++ = last = *p;

	return q - buffer;
}

int SparseMatrix::get_num_indices()
{
	_F_
	int total = 0;
	for (int i = 0; i < size; i++)
		for (Page *page = pages[i]; page != NULL; page = page->next)
			total += page->count;

	return total;
}
