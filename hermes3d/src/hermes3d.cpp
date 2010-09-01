// This file is part of Hermes3D
//
// Copyright (c) 2010 hp-FEM group at the University of Nevada, Reno (UNR).
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

#include <stdio.h>

// verbose output by default
bool verbose = true;

// print out the banner
void banner() {
	if (verbose) {
		printf("----------------------------------------------\n");
		printf("  This is Hermes3D - a C++ library for rapid\n");
		printf("development of adaptive FEM and hp-FEM solvers\n");
		printf("      developed by the hp-FEM group at UNR\n");
		printf("     and distributed under the GPL license.\n");
		printf("    For more details visit http://hpfem.org/.\n");
		printf("----------------------------------------------\n");
		printf("       To turn off verbose output, call\n");
		printf("       set_verbose(false) in your code.\n");
	}
}

void set_verbose(bool verb) {
	verbose = verb;
}
