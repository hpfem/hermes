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


// TODO: improve the way how to compare two meshes (compare outputs from dump() is not enough)

#include <hermes3d.h>
#include <common/trace.h>
#include <common/callstack.h>

#define ERROR_SUCCESS								0
#define ERROR_FAILURE								-1
#define ERROR_NOT_ENOUGH_PARAMS						-2

bool testPrint(bool value, const char *msg, bool correct)
{
	_F_
	printf("%s...", msg);
	if (value == correct) {
		printf("OK\n");
		return true;
	}
	else {
		printf("failed\n");
		return false;
	}
}

//
// tests themselves
//

int test_mesh3d_loader(char *file_name)
{
	_F_
	Mesh mesh;
	Mesh3DReader mloader;
	if (mloader.load(file_name, &mesh)) {
		mesh.dump();
		return ERROR_SUCCESS;
	}
	else {
		printf("failed\n");
		return ERROR_FAILURE;
	}
}

int test_hdf5_loader(char *file_name)
{
	_F_
	Mesh mesh;
	HDF5Reader mloader;
	if (mloader.load(file_name, &mesh)) {
		mesh.dump();
		return ERROR_SUCCESS;
	}
	else {
		printf("failed\n");
		return ERROR_FAILURE;
	}
}

int test_exodusii_loader(char *file_name)
{
	_F_
	Mesh mesh;
	ExodusIIReader mloader;
	if (mloader.load(file_name, &mesh)) {
		mesh.dump();
		return ERROR_SUCCESS;
	}
	else {
		printf("failed\n");
		return ERROR_FAILURE;
	}
}

int main(int argc, char *argv[])
{
	_F_
	int ret = ERROR_SUCCESS;
	set_verbose(false);

	if (argc < 3)
		return ERROR_NOT_ENOUGH_PARAMS;

	if (strcmp(argv[1], "m3d") == 0) {
		ret = test_mesh3d_loader(argv[2]);
	}
	else if (strcmp(argv[1], "hdf5") == 0) {
		ret = test_hdf5_loader(argv[2]);
	}
	else if (strcmp(argv[1], "exoii") == 0) {
		ret = test_exodusii_loader(argv[2]);
	}

	return ret;
}
