// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// This file was written by:
// - David Andrs
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

//
// multi-mesh.cc
//
// -Laplace u1 + u2      = f1
//      -Laplace u2 + u3 = f2
//           -Laplace u3 = f3
//
// u1 = x^2 + y^2 + z^2
// u2 = (1-x^2)(1 - y^2)(1- z^2)
// u3 = x
//

#include "config.h"
#include <math.h>
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>
#include <common/utils.h>

// error should be smaller than this epsilon
#define EPS								10e-10F

#ifdef RHS2

double fnc(double x, double y, double z)
{
	return x*x + y*y + z*z;
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return fnc(x, y, z);
}

BCType bc_types(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return fnc(x, y, z);
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                    user_data_t<res_t> *data)
{
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (data->ext[0].fn[i] * v->fn[i]);
	return res;
}

#elif defined SYS

template<typename T>
T u1(T x, T y, T z)
{
	return (x*x + y*y + z*z);
//	return (1 - x*x) * (1 - y*y) * (1 - z*z);
}

template<typename T>
T u2(T x, T y, T z)
{
	return x;
//	return 3.0;
//	return (1 - x*x) * x*x * (1 - y*y) * y*y * (1 - z*z) * z*z;
}

// needed for calculation norms and used by visualizator
double exact_sln_fn_1(double x, double y, double z, double &dx, double &dy, double &dz)
{
//	dx = -2 * x * (1 - y*y) * (1 - z*z);
//	dy = -2 * (1 - x*x) * y * (1 - z*z);
//	dz = -2 * (1 - x*x) * (1 - y*y) * z;
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return u1(x, y, z);
}

double exact_sln_fn_2(double x, double y, double z, double &dx, double &dy, double &dz)
{
//	dx = 2 * x * (1 - x*x) * y*y * (1 - y*y) * z*z * (1 - z*z) - 2 * x*x*x * y*y * (1 - y*y) * z*z * (1 - z*z);
//	dy = 2 * x*x * (1 - x*x) * y * (1 - y*y) * z*z * (1 - z*z) - 2 * x*x * (1 - x*x) * y*y*y * z*z * (1 - z*z);
//	dz = 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * z * (1 - z*z) - 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * z*z*z;
	dx = 1;
	dy = 0;
	dz = 0;

	return u2(x, y, z);
}

//

BCType bc_types_1(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values_1(int ess_bdy_marker, double x, double y, double z) {
	return u1(x, y, z);
}

BCType bc_types_2(int marker)
{
	if (marker == 3) return BC_NATURAL;
	else return BC_ESSENTIAL;
}

scalar essential_bc_values_2(int ess_bdy_marker, double x, double y, double z) {
	return u2(x, y, z);
}

template<typename f_t, typename res_t>
res_t bilinear_form_1_1(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                        user_data_t<res_t> *data)
{
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t bilinear_form_1_2(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                        user_data_t<res_t> *data)
{
//	return 0.0;
	return int_u_v<f_t, res_t>(n, wt, u, v, e);
}

//template<typename T>
//T f1(T x, T y, T z)
//{
//	T ddxx = -2 * (1 - y*y) * (1 - z*z);
//	T ddyy = -2 * (1 - x*x) * (1 - z*z);
//	T ddzz = -2 * (1 - x*x) * (1 - y*y);
//
//	return -(ddxx + ddyy + ddzz) + u2(x, y, z);
//}

template<typename f_t, typename res_t>
res_t linear_form_1(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data)
{
//	return int_F_v<f_t, res_t>(n, wt, f1, u, e);
//	return -3.0 * int_u<f_t, res_t>(n, wt, u, e);
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * ((-6.0 + (e->x[i]) * u->fn[i]);
	return res;
//	int_u<f_t, res_t>(n, wt, u, e);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

//template<typename f_t, typename res_t>
//res_t bilinear_form_2_1(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
//                        user_data_t<res_t> *data)
//{
//	return int_u_v<f_t, res_t>(n, wt, u, v, e);
//}

template<typename f_t, typename res_t>
res_t bilinear_form_2_2(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                        user_data_t<res_t> *data)
{
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

//template<typename T>
//T f2(T x, T y, T z)
//{
//	T ddxx = 2 * (1 - x*x) * y*y * (1 - y*y) * z*z * (1 - z*z) - 10 * x*x * y*y * (1 - y*y) * z*z * (1 - z*z);
//	T ddyy = 2 * x*x * (1 - x*x) * (1 - y*y) * z*z * (1 - z*z) - 10 * x*x * (1 - x*x) * y*y * z*z * (1 - z*z);
//	T ddzz = 2 * x*x * (1 - x*x) * y*y * (1 - y*y) * (1 - z*z) - 10 * x*x * (1 - x*x) * y*y * (1 - y*y) * z*z;
//
//	return -(ddxx + ddyy + ddzz) + u1(x, y, z);
//}

template<typename f_t, typename res_t>
res_t linear_form_2(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data)
{
//	return int_F_v<f_t, res_t>(n, wt, f2, u, e);
//	return -6.0 * int_u<f_t, res_t>(n, wt, u, e);
	return 0;
}

template<typename f_t, typename res_t>
res_t linear_form_2_surf(int np, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data)
{
	res_t result = 0;
	for (int i = 0; i < np; i++) {
		res_t dx = 1;
		res_t dy = 0;
		res_t dz = 0;

		result += wt[i] * (u->fn[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i]));
	}
	return result;
}

#elif defined SYS3

template<typename T>
T u1(T x, T y, T z)
{
	return (x*x + y*y + z*z);
}

template<typename T>
T u2(T x, T y, T z)
{
	return (1 - x*x) * (1 - y*y) * (1 - z*z);
}

template<typename T>
T u3(T x, T y, T z)
{
	return x;
}

double exact_sln_fn_1(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;

	return u1(x, y, z);
}

double exact_sln_fn_2(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = -2 * x * (1 - y*y) * (1 - z*z);
	dy = -2 * (1 - x*x) * y * (1 - z*z);
	dz = -2 * (1 - x*x) * (1 - y*y) * z;

	return u2(x, y, z);
}

double exact_sln_fn_3(double x, double y, double z, double &dx, double &dy, double &dz)
{
	dx = 1;
	dy = 0;
	dz = 0;

	return u3(x, y, z);
}

//

BCType bc_types_1(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values_1(int ess_bdy_marker, double x, double y, double z) {
	return u1(x, y, z);
}

BCType bc_types_2(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values_2(int ess_bdy_marker, double x, double y, double z)
{
	return 0;
}

BCType bc_types_3(int marker)
{
	return BC_ESSENTIAL;
}

scalar essential_bc_values_3(int ess_bdy_marker, double x, double y, double z)
{
	return x;
}

// 1. eqn ------------------------------------------------------------------------------------------

template<typename f_t, typename res_t>
res_t biform_1_1(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                 user_data_t<res_t> *data)
{
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t biform_1_2(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                 user_data_t<res_t> *data)
{
	return int_u_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename T>
T rhs1(T x, T y, T z)
{
	return -6.0 + u2(x, y, z);
}

template<typename f_t, typename res_t>
res_t liform_1(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data)
{
	res_t res = 0.0;
	for (int i = 0; i < n; i++)
		res += wt[i] * (rhs1(e->x[i], e->y[i], e->z[i]) * v->fn[i]);
	return res;
}

// 2. eqn ------------------------------------------------------------------------------------------

template<typename f_t, typename res_t>
res_t biform_2_2(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                 user_data_t<res_t> *data)
{
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t biform_2_3(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                 user_data_t<res_t> *data)
{
	return int_u_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename T>
T rhs2(T x, T y, T z)
{
	T laplace = 2 * (1 - y*y) * (1 - z*z) + 2 * (1 - x*x) * (1 - z*z) + 2 * (1 - x*x) * (1 - y*y);
	return laplace + u3(x, y, z);
}

template<typename f_t, typename res_t>
res_t liform_2(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data)
{
	return int_F_v<f_t, res_t>(n, wt, rhs2, v, e);
}

// 3. eqn ------------------------------------------------------------------------------------------

template<typename f_t, typename res_t>
res_t biform_3_3(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                        user_data_t<res_t> *data)
{
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}


#endif

void out_fn(Solution *sln, const char *name)
{
#ifdef OUTPUT_DIR
	char of_name[512];
	sprintf(of_name, "%s/%s.vtk", OUTPUT_DIR, name);
	FILE *ofile = fopen(of_name, "w");
	if (ofile != NULL) {
		VtkOutputEngine output(ofile);
		output.out(sln, "Uh", FN_VAL_0);
		fclose(ofile);
	}
	else {
		warning("Can not open '%s' for writing.", of_name);
	}
#endif
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **args)
{
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &args, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	if (argc < 2) error("Not enough parameters");

	printf("* Loading mesh '%s'\n", args[1]);
	Mesh mesh1;
	Mesh3DReader mesh_loader;
	if (!mesh_loader.load(args[1], &mesh1)) error("Loading mesh file '%s'\n", args[1]);

	H1ShapesetLobattoHex shapeset;
#if defined RHS2
	printf("* Setting the space up\n");
	H1Space space(&mesh1, &shapeset);
	space.set_bc_types(bc_types);
	space.set_essential_bc_values(essential_bc_values);

	order3_t order(2, 2, 2);
	printf("  - Setting uniform order to (%d, %d, %d)\n", order.x, order.y, order.z);
	space.set_uniform_order(order);

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	printf("* Calculating a solution\n");

	// duplicate the mesh
	Mesh mesh2;
	mesh2.copy(mesh1);
	// do some changes
	mesh2.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);
	mesh2.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

	Solution fsln(&mesh2);
	fsln.set_const(-6.0);
#else
	// duplicate the mesh
	Mesh mesh2;
	mesh2.copy(mesh1);

	Mesh mesh3;
	mesh3.copy(mesh1);

	// change meshes
	mesh1.refine_all_elements(H3D_REFT_HEX_X);
	mesh2.refine_all_elements(H3D_REFT_HEX_Y);
	mesh3.refine_all_elements(H3D_REFT_HEX_Z);

	printf("* Setup spaces\n");
	H1Space space1(&mesh1, &shapeset);
	space1.set_bc_types(bc_types_1);
	space1.set_essential_bc_values(essential_bc_values_1);

	order3_t o1(2, 2, 2);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o1.x, o1.y, o1.z);
	space1.set_uniform_order(o1);

	H1Space space2(&mesh2, &shapeset);
	space2.set_bc_types(bc_types_2);
	space2.set_essential_bc_values(essential_bc_values_2);

	order3_t o2(2, 2, 2);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o2.x, o2.y, o2.z);
	space2.set_uniform_order(o2);

	H1Space space3(&mesh3, &shapeset);
	space3.set_bc_types(bc_types_3);
	space3.set_essential_bc_values(essential_bc_values_3);

	order3_t o3(1, 1, 1);
	printf("  - Setting uniform order to (%d, %d, %d)\n", o3.x, o3.y, o3.z);
	space3.set_uniform_order(o3);

	int ndofs = 0;
	ndofs += space1.assign_dofs();
	ndofs += space2.assign_dofs(ndofs);
	ndofs += space3.assign_dofs(ndofs);
	printf("  - Number of DOFs: %d\n", ndofs);
#endif

#if defined WITH_UMFPACK
	UMFPackMatrix mat;
	UMFPackVector rhs;
	UMFPackLinearSolver solver(&mat, &rhs);
#elif defined WITH_PARDISO
	PardisoMatrix mat;
	PardisoVector rhs;
	PardisoLinearSolver solver(&mat, &rhs);
#elif defined WITH_PETSC
	PetscMatrix mat;
	PetscVector rhs;
	PetscLinearSolver solver(&mat, &rhs);
#elif defined WITH_MUMPS
	MumpsMatrix mat;
	MumpsVector rhs;
	MumpsSolver solver(&mat, &rhs);
#endif

#ifdef RHS2
	WeakForm wf;
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>, ANY, &fsln);

	LinearProblem lp(&wf);
	lp.set_space(&space);
#elif defined SYS3
	WeakForm wf(3);
	wf.add_matrix_form(0, 0, biform_1_1<double, scalar>, biform_1_1<ord_t, ord_t>, SYM);
	wf.add_matrix_form(0, 1, biform_1_2<double, scalar>, biform_1_2<ord_t, ord_t>, UNSYM);
	wf.add_vector_form(0, liform_1<double, scalar>, liform_1<ord_t, ord_t>);

	wf.add_matrix_form(1, 1, biform_2_2<double, scalar>, biform_2_2<ord_t, ord_t>, SYM);
	wf.add_matrix_form(1, 2, biform_2_3<double, scalar>, biform_2_3<ord_t, ord_t>, UNSYM);
	wf.add_vector_form(1, liform_2<double, scalar>, liform_2<ord_t, ord_t>);

	wf.add_matrix_form(2, 2, biform_3_3<double, scalar>, biform_3_3<ord_t, ord_t>, SYM);

	LinearProblem lp(&wf);
	lp.set_spaces(Tuple<Space *>(&space1, &space2, &space3));
#endif

	// assemble stiffness matrix
	printf("  - assembling... "); fflush(stdout);
	Timer assemble_timer;
	assemble_timer.start();
	lp.assemble(&mat, &rhs);
	assemble_timer.stop();
	printf("%s (%lf secs)\n", assemble_timer.get_human_time(), assemble_timer.get_seconds());

	// solve the stiffness matrix
	printf("  - solving... "); fflush(stdout);
	Timer solve_timer;
	solve_timer.start();
	bool solved = solver.solve();
	solve_timer.stop();
	printf("%s (%lf secs)\n", solve_timer.get_human_time(), solve_timer.get_seconds());

	if (solved) {
#ifdef RHS2
		Timer sln_pre_tmr;
		Solution sln(&mesh1);
		sln_pre_tmr.start();
		sln.set_fe_solution(&space, solver.get_solution());
		sln_pre_tmr.stop();

		printf("* Solution:\n");

		ExactSolution ex_sln(&mesh1, exact_solution);
		// norm
		double h1_sln_norm = h1_norm(&sln);
		double h1_err_norm = h1_error(&sln, &ex_sln);
		printf("  - H1 solution norm:   % le\n", h1_sln_norm);
		printf("  - H1 error norm:      % le\n", h1_err_norm);

		double l2_sln_norm = l2_norm(&sln);
		double l2_err_norm = l2_error(&sln, &ex_sln);
		printf("  - L2 solution norm:   % le\n", l2_sln_norm);
		printf("  - L2 error norm:      % le\n", l2_err_norm);

		if (h1_err_norm > EPS || l2_err_norm > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}
#elif defined SYS3
		// solution 1
		Solution sln1(&mesh1);
		Solution sln2(&mesh2);
		Solution sln3(&mesh3);

		sln1.set_fe_solution(&space1, solver.get_solution());
		sln2.set_fe_solution(&space2, solver.get_solution());
		sln3.set_fe_solution(&space3, solver.get_solution());

		ExactSolution esln1(&mesh1, exact_sln_fn_1);
		ExactSolution esln2(&mesh2, exact_sln_fn_2);
		ExactSolution esln3(&mesh3, exact_sln_fn_3);

		// norm
		double h1_err_norm1 = h1_error(&sln1, &esln1);
		double h1_err_norm2 = h1_error(&sln2, &esln2);
		double h1_err_norm3 = h1_error(&sln3, &esln3);

		double l2_err_norm1 = l2_error(&sln1, &esln1);
		double l2_err_norm2 = l2_error(&sln2, &esln2);
		double l2_err_norm3 = l2_error(&sln3, &esln3);

		printf("  - H1 error norm:      % le\n", h1_err_norm1);
		printf("  - L2 error norm:      % le\n", l2_err_norm1);
		if (h1_err_norm1 > EPS || l2_err_norm1 > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}

		printf("  - H1 error norm:      % le\n", h1_err_norm2);
		printf("  - L2 error norm:      % le\n", l2_err_norm2);
		if (h1_err_norm2 > EPS || l2_err_norm2 > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}

		printf("  - H1 error norm:      % le\n", h1_err_norm3);
		printf("  - L2 error norm:      % le\n", l2_err_norm3);
		if (h1_err_norm3 > EPS || l2_err_norm3 > EPS) {
			// calculated solution is not enough precise
			res = ERR_FAILURE;
		}
#endif

#ifdef RHS2
		out_fn(&sln, "solution");
#elif defined SYS3
		out_fn(&sln1, "sln1");
		out_fn(&sln2, "sln2");
		out_fn(&sln3, "sln3");
#endif
	}
	else
		res = ERR_FAILURE;

#ifdef WITH_PETSC
	mat.free();
	rhs.free();

	PetscFinalize();
#endif

	return res;
}
