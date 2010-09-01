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

#include "config.h"
#include <math.h>
#include <hermes3d.h>

#if defined GMSH
	GmshOutputEngine output(stdout);
#elif defined VTK
	VtkOutputEngine output(stdout, 1);
#endif


int m, n, o;

// error should be smaller than this epsilon
#define EPSILON								10e-10F

double fnc(double x, double y, double z)
{
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
}

template<typename T>
T dfnc(T x, T y, T z)
{
	T ddxx = m * (m - 1) * pow(x, m - 2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
	T ddyy = n * (n - 1) * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = o * (o - 1) * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);

	return -(ddxx + ddyy + ddzz);
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz)
{
	// u(x, y, z) = x^m * y^n * z^o + x^2 * y^3 - x^3 * z + z^4
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 3) + 4 * pow(z, 3);

	return fnc(x, y, z);
}

// exact solution for vector functions
double exact_solution0(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = dy = dz = 0;
	return 1.0;
}

double exact_solution1(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = dy = dz = 0;
	return 0.0;
}

double exact_solution2(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = dy = dz = 0;
	return 0.0;
}

//
scalar3 &exact_vec_solution(double x, double y, double z, scalar3 &dx, scalar3 &dy, scalar3 &dz)
{
	static scalar3 val;

	dx[0] = dx[1] = dx[2] = 0;
	dy[0] = dy[1] = dy[2] = 0;
	dz[0] = dz[1] = dz[2] = 0;

	val[0] = 0.0;
	val[1] = 1.0;
	val[2] = 0.0;

	return val;
}


//

BCType bc_types(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return fnc(x, y, z);
}

template<typename f_t, typename res_t>
res_t bilinear_form(int np, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                    user_data_t<res_t> *data)
{
	return int_grad_u_grad_v<f_t, res_t>(np, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data)
{
	return int_F_v<f_t, res_t>(n, wt, dfnc, u, e);
}

#if defined WITH_UMFPACK
#define StiffMatrix		UMFPackMatrix
#elif defined WITH_PARDISO
#define StiffMatrix		PardisoMatrix
#elif defined WITH_PETSC
#define StiffMatrix		PetscMatrix
#elif defined WITH_MUMPS
#define StiffMatrix		MumpsMatrix
#elif defined WITH_TRILINOS
#define StiffMatrix		EpetraMatrix
#endif

void test_mat(Mesh *mesh, StiffMatrix &mat)
{
#if defined WITH_UMFPACK
	UMFPackVector rhs;
	UMFPackLinearSolver solver(&mat, &rhs);
#elif defined WITH_PARDISO
	PardisoVector rhs;
	PardisoLinearSolver solver(&mat, &rhs);
#elif defined WITH_PETSC
	PetscVector rhs;
	PetscLinearSolver solver(&mat, &rhs);
#elif defined WITH_MUMPS
	MumpsVector rhs;
	MumpsSolver solver(&mat, &rhs);
#elif defined WITH_TRILINOS
	EpetraVector rhs;
	AztecOOSolver solver(&mat, &rhs);
#endif

	H1ShapesetLobattoHex shapeset;
	H1Space space(mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_essential_bc_values(essential_bc_values);

	m = n = o = 2;
	int mx = maxn(4, m, n, o, 4);
	order3_t order(mx, mx, mx);
	space.set_uniform_order(order);

	space.assign_dofs();

	WeakForm wf(1);
	wf.add_matrix_form(0, 0, bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM);
	wf.add_vector_form(0, linear_form<double, scalar>, linear_form<ord_t, ord_t>);

	LinearProblem lp(&wf);
	lp.set_space(&space);

	// assemble stiffness matrix
	lp.assemble(&mat, &rhs);
	solver.solve();
}

void test_mm(Mesh *mesh)
{
	_F_
	// testing a visualization of a vector-valued solution, where each component is on
	// a different mesh
	Mesh mesh0, mesh1, mesh2;
	mesh0.copy(*mesh);
	mesh1.copy(*mesh);
	mesh2.copy(*mesh);

	mesh0.refine_element(1, H3D_REFT_HEX_X);
	mesh1.refine_element(1, H3D_REFT_HEX_Y);
	mesh2.refine_element(1, H3D_REFT_HEX_Z);

	ExactSolution ex_sln0(&mesh0, exact_solution0);
	ExactSolution ex_sln1(&mesh1, exact_solution1);
	ExactSolution ex_sln2(&mesh2, exact_solution2);
	output.out(&ex_sln0, &ex_sln1, &ex_sln2, "U");
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args)
{
	set_verbose(false);

	if (argc < 3) error("Not enough parameters");

	char *type = args[1];

	Mesh mesh;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[2], &mesh)) error("Loading mesh file '%s'\n", args[2]);

	if (strcmp(type, "sln") == 0) {
		// Testing on Exact solution which always gives the same value (values from Solution may differ by epsilon)
		ExactSolution ex_sln(&mesh, exact_solution);
		output.out(&ex_sln, "U");
	}
	else if (strcmp(type, "vec-sln") == 0) {
		// Testing on Exact solution which always gives the same value (values from Solution may differ by epsilon)
		ExactSolution ex_sln(&mesh, exact_vec_solution);
		output.out(&ex_sln, "U");
	}
	else if (strcmp(type, "3sln") == 0) {
		// Testing on Exact solution which always gives the same value (values from Solution may differ by epsilon)
		ExactSolution ex_sln0(&mesh, exact_solution0);
		ExactSolution ex_sln1(&mesh, exact_solution1);
		ExactSolution ex_sln2(&mesh, exact_solution2);
		output.out(&ex_sln0, &ex_sln1, &ex_sln2, "U");
	}
	else if (strcmp(type, "ord") == 0) {
		H1ShapesetLobattoHex shapeset;

		H1Space space(&mesh, &shapeset);
		space.set_bc_types(bc_types);
		space.set_essential_bc_values(essential_bc_values);

		order3_t order;
		if (mesh.elements[1]->get_mode() == MODE_HEXAHEDRON)
			order = order3_t(2, 3, 4);
		else if (mesh.elements[1]->get_mode() == MODE_TETRAHEDRON)
			order = order3_t(3);
		else
			error(H3D_ERR_NOT_IMPLEMENTED);
		space.set_uniform_order(order);

		output.out_orders(&space, "orders");
	}
	else if (strcmp(type, "bc") == 0) {
		output.out_bc(&mesh);
	}
	else if (strcmp(type, "mat") == 0) {
		StiffMatrix mat;
		test_mat(&mesh, mat);
		output.out(&mat);
	}
	else if (strcmp(type, "mm") == 0) {
		test_mm(&mesh);
	}

	return 0;
}
