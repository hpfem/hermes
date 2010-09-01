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

/*
 * main.cpp
 *
 * Lightweight solver that:
 *  - reads the file with matrix and rhs (in hermes2d binary format)
 *  - solves it (using UMFPack)
 *  - outputs the result
 *
 * Usage: h-solve <infile> <outfile>
 *   <infile> input file with matrix and rhs (in hermes2d binary format)
 *   <outfile> output file with result (in hermes2d binary format)
 *
 */

#include <stdio.h>
#include <string.h>

extern "C" {
#include <umfpack.h>
}

int *Ap = NULL;
int *Ai = NULL;
double *Ax = NULL;
double *rhs = NULL;
double *sln = NULL;

int ndofs;

int readFromBinaryFile(const char *file) {
	FILE *fin = fopen(file, "rb");
	if (fin == NULL) {
		printf("Can not open file '%s'.\n", file);
		return -1;
	}

	char buf[8];
	int version;
	int sizeofscalar;
	int nnz;

	// matrix
	fread(buf, sizeof(char), 8, fin);
	if (strncmp(buf, "H2DX", 4) != 0) {
		printf("Matrix file corrupted.\n");
		return -1;
	}

	fread(&sizeofscalar, sizeof(int), 1, fin);
	fread(&ndofs, sizeof(int), 1, fin);
	fread(&nnz, sizeof(int), 1, fin);

	Ap = new int [ndofs + 1];
	Ai = new int [nnz];
	Ax = new double [nnz];

	fread(Ap, sizeof(int), ndofs + 1, fin);
	fread(Ai, sizeof(int), nnz, fin);
	fread(Ax, sizeofscalar, nnz, fin);

	// rhs
	fread(buf, sizeof(char), 8, fin);
	if (strncmp(buf, "H2DR", 4) != 0) {
		printf("Matrix file corrupted.\n");
		return -1;
	}
	int ssize;
	fread(&ssize, sizeof(int), 1, fin);
	fread(&ndofs, sizeof(int), 1, fin);

	rhs = new double [ndofs];
	fread(rhs, ssize, ndofs, fin);

	fclose(fin);

	return 0;
}

void dump_rhs(const char *file_name) {
	FILE *file = fopen(file_name, "wb");
	if (file == NULL) {
		printf("Can not open file '%s' for writing.\n", file_name);
		return;
	}

	fwrite("H2DR\001\000\000\000", 1, 8, file);
	int ssize = sizeof(double);
	fwrite(&ssize, sizeof(int), 1, file);
	fwrite(&ndofs, sizeof(int), 1, file);
	fwrite(sln, ssize, ndofs, file);
}

int main(int argc, char *argv[]) {
	if (argc < 3) {
		printf("Not enough parameters.\n");
		return -1;
	}

	if (readFromBinaryFile(argv[1]) < 0) {
		return -1;
	}

	printf("symbolic analysis\n");
	void *symbolic, *numeric;
	if (umfpack_di_symbolic(ndofs, ndofs, Ap, Ai, Ax, &symbolic, NULL, NULL) < 0) {
		printf("umfpack_di_symbolic failed\n");
		return -1;
	}
	if (symbolic == NULL) {
		printf("umfpack_di_symbolic error: symbolic == NULL\n");
		return -1;
	}

	printf("numeric analysis\n");
	if (umfpack_di_numeric(Ap, Ai, Ax, symbolic, &numeric, NULL, NULL) < 0) {
		printf("umfpack_di_numeric failed\n");
		return -1;
	}
	if (numeric == NULL) {
		printf("umfpack_di_numeric error: numeric == NULL\n");
		return -1;
	}

	printf("solving\n");
	sln = new double [ndofs];

	if (umfpack_di_solve(UMFPACK_A, Ap, Ai, Ax, sln, rhs, numeric, NULL, NULL) < 0) {
		printf("umfpack_di_solve failed\n");
		return -1;
	}

	umfpack_di_free_symbolic(&symbolic);
	umfpack_di_free_numeric(&numeric);

	dump_rhs(argv[2]);

	printf("done\n");

	return 0;
}

