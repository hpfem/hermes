// This file is part of Hermes3D
//
// Copyright (c) 2010 hp-FEM group at the University of Nevada, Reno (UNR).
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

// Test to verify the hp-adaptivity
//

#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <hermes3d.h>
#include <common/trace.h>
#include <common/timer.h>
#include <common/error.h>

#undef REFERENCE_SOLUTION				// use ref. solution for guiding the adaptivity

const double TOLERANCE = 0.000001;		// error tolerance in percent
const double THRESHOLD = 0.3;			// error threshold for element refinement

#define ANISO_X						1
#define ANISO_Y						2
#define ANISO_Z						4

int init_q = 1;							// initial order
int init_p = 2;							// initial order
int aniso_type = 0;

// error should be smaller than this epsilon
#define EPS								10e-10F

double fnc(double x, double y, double z)
{
	switch (aniso_type) {
		case ANISO_X: return sin(x);
		case ANISO_Y: return sin(y);
		case ANISO_Z: return sin(z);

		case ANISO_X | ANISO_Y: return sin(x) * sin(y);
		case ANISO_X | ANISO_Z: return sin(x) * sin(z);
		case ANISO_Y | ANISO_Z: return sin(y) * sin(z);

		case ANISO_X | ANISO_Y | ANISO_Z: return sin(x) * sin(y) * sin(z);

		default:
			return 0;
	}
}

template<typename T>
T rhs(T x, T y, T z)
{
	T ddxx = 0;
	T ddyy = 0;
	T ddzz = 0;

	switch (aniso_type) {
		case ANISO_X: ddxx = -sin(x); break;
		case ANISO_Y: ddyy = -sin(y); break;
		case ANISO_Z: ddzz = -sin(z); break;

		case ANISO_X | ANISO_Y:
			ddxx = - sin(x) * sin(y);
			ddyy = - sin(x) * sin(y);
			break;

		case ANISO_X | ANISO_Z:
			ddxx = - sin(x) * sin(z);
			ddzz = - sin(x) * sin(z);
			break;

		case ANISO_Y | ANISO_Z:
			ddyy = - sin(y) * sin(z);
			ddzz = - sin(y) * sin(z);
			break;

		case ANISO_X | ANISO_Y | ANISO_Z:
			ddxx = - sin(x) * sin(y) * sin(z);
			ddyy = - sin(x) * sin(y) * sin(z);
			ddzz = - sin(x) * sin(y) * sin(z);
			break;
	}

	return ddxx + ddyy + ddzz;
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz)
{
	switch (aniso_type) {
		case ANISO_X:
			dx = cos(x);
			dy = 0;
			dz = 0;
			break;

		case ANISO_Y:
			dx = 0;
			dy = cos(y);
			dz = 0;
			break;

		case ANISO_Z:
			dx = 0;
			dy = 0;
			dz = cos(z);
			break;

		case ANISO_X | ANISO_Y:
			dx = cos(x) * sin(y);
			dy = sin(x) * cos(y);
			dz = 0;
			break;

		case ANISO_X | ANISO_Z:
			dx = cos(x) * sin(z);
			dy = 0;
			dz = sin(x) * cos(z);
			break;

		case ANISO_Y | ANISO_Z:
			dx = 0;
			dy = cos(y) * sin(z);
			dz = sin(y) * cos(z);
			break;

		case ANISO_X | ANISO_Y | ANISO_Z:
			dx = cos(x) * sin(y) * sin(z);
			dy = sin(x) * cos(y) * sin(z);
			dz = sin(x) * sin(y) * cos(z);
			break;
	}

	return fnc(x, y, z);
}

//

BCType bc_types(int marker)
{
	switch (aniso_type) {
		case ANISO_X:
			if (marker == 1 || marker == 2) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_Y:
			if (marker == 3 || marker == 4) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_Z:
			if (marker == 5 || marker == 6) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_X | ANISO_Y:
			if (marker == 1 || marker == 2 || marker == 3 || marker == 4) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_X | ANISO_Z:
			if (marker == 1 || marker == 2 || marker == 5 || marker == 6) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_Y | ANISO_Z:
			if (marker == 3 || marker == 4 || marker == 5 || marker == 6) return BC_ESSENTIAL;
			else return BC_NATURAL;

		case ANISO_X | ANISO_Y | ANISO_Z:
			return BC_ESSENTIAL;
	}

	return BC_ESSENTIAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return 0;
}

template<typename f_t, typename res_t>
res_t bilinear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e, user_data_t<res_t> *data)
{
	return int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data)
{
	return -int_F_v<f_t, res_t>(n, wt, rhs, u, e);
}

void out_fn(MeshFunction *x, const char *name, int i) {
#ifdef OUTPUT_DIR
	char of_name[1024];
	FILE *ofile;
	// mesh out
	sprintf(of_name, "%s/%s-%d.vtk", OUTPUT_DIR, name, i);
	ofile = fopen(of_name, "w");
	if (ofile != NULL) {
		VtkOutputEngine output(ofile);
		output.out(x, name);
		fclose(ofile);
	}
	else {
		warning("Can not open '%s' for writing.", of_name);
	}
#endif
}

void parse_aniso_type(char *str)
{
	int type = 0;
	if (strchr(str, 'x') != NULL) type |= ANISO_X;
	if (strchr(str, 'y') != NULL) type |= ANISO_Y;
	if (strchr(str, 'z') != NULL) type |= ANISO_Z;
	aniso_type = type;
}

bool check_order(const order3_t &spord)
{
	switch (aniso_type)
	{
		case ANISO_X:
			if (spord.y == init_q && spord.z == init_q) return true;
			break;

		case ANISO_Y:
			if (spord.x == init_q && spord.z == init_q) return true;
			break;

		case ANISO_Z:
			if (spord.x == init_q && spord.y == init_q) return true;
			break;

		case ANISO_X | ANISO_Y:
			if (spord.z == init_q) return true;
			break;

		case ANISO_X | ANISO_Z:
			if (spord.y == init_q) return true;
			break;

		case ANISO_Y | ANISO_Z:
			if (spord.x == init_q) return true;
			break;

		case ANISO_X | ANISO_Y | ANISO_Z:
			return true;
	}

	return false;
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
	int res = ERR_SUCCESS;

#ifdef WITH_PETSC
	PetscInitialize(&argc, &argv, (char *) PETSC_NULL, PETSC_NULL);
#endif
	set_verbose(false);

	if (argc < 3) error("Not enough parameters.");

	parse_aniso_type(argv[2]);

	H1ShapesetLobattoHex shapeset;

	printf("* Loading mesh '%s'\n", argv[1]);
	Mesh mesh;
	Mesh3DReader mloader;
	if (!mloader.load(argv[1], &mesh)) error("Loading mesh file '%s'\n", argv[1]);

	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_essential_bc_values(essential_bc_values);

	order3_t order;
	switch (aniso_type)
	{
		case ANISO_X:
			order = order3_t(init_p, init_q, init_q);
			break;

		case ANISO_Y:
			order = order3_t(init_q, init_p, init_q);
			break;

		case ANISO_Z:
			order = order3_t(init_q, init_q, init_p);
			break;

		case ANISO_X | ANISO_Y:
			order = order3_t(init_p, init_p, init_q);
			break;

		case ANISO_X | ANISO_Z:
			order = order3_t(init_p, init_q, init_p);
			break;

		case ANISO_Y | ANISO_Z:
			order = order3_t(init_q, init_p, init_p);
			break;

		case ANISO_X | ANISO_Y | ANISO_Z:
			order = order3_t(init_p, init_p, init_p);
			break;
	}
	printf("  - Setting uniform order to (%d, %d, %d)\n", order.x, order.y, order.z);
	space.set_uniform_order(order);

	WeakForm wf;
	wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<ord_t, ord_t>, SYM, ANY);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>, ANY);

	LinearProblem lp(&wf);
	lp.set_space(&space);

	bool done = false;
	int iter = 0;
	do {
		Timer assemble_timer("Assembling stiffness matrix");
		Timer solve_timer("Solving stiffness matrix");

		printf("\n=== Iter #%d ================================================================\n", iter);

		// check the we are doing all right
		FOR_ALL_ACTIVE_ELEMENTS(eid, &mesh) {
			order3_t spord = space.get_element_order(eid);
			printf("#%ld: order = (%d, %d, %d)\n", eid, spord.x, spord.y, spord.z);
		}

		if (mesh.get_num_elements() != 1) {
			res = ERR_FAILURE;
			printf("failed\n");
			break;
		}

		order3_t spord = space.get_element_order(1);
		if (!check_order(spord)) {
			res = ERR_FAILURE;
			printf("failed\n");
			break;
		}

		// we're good -> go ahead
		printf("\nSolution\n");

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

		int ndofs = space.assign_dofs();
		printf("  - Number of DOFs: %d\n", ndofs);

		// assemble stiffness matrix
		printf("  - Assembling... "); fflush(stdout);
		assemble_timer.reset();
		assemble_timer.start();
		lp.assemble(&mat, &rhs);
		assemble_timer.stop();
		printf("done in %s (%lf secs)\n", assemble_timer.get_human_time(), assemble_timer.get_seconds());

		// solve the stiffness matrix
		printf("  - Solving... "); fflush(stdout);
		solve_timer.reset();
		solve_timer.start();
		bool solved = solver.solve();
		solve_timer.stop();
		if (solved)
			printf("done in %s (%lf secs)\n", solve_timer.get_human_time(), solve_timer.get_seconds());
		else {
			res = ERR_FAILURE;
			printf("failed\n");
			break;
		}

		Solution sln(&mesh);
		sln.set_fe_solution(&space, solver.get_solution());

		printf("Reference solution\n");

		Mesh rmesh;
		rmesh.copy(mesh);
		rmesh.refine_all_elements(H3D_H3D_H3D_REFT_HEX_XYZ);

#ifdef REFERENCE_SOLUTION

#if defined WITH_UMFPACK
		UMFPackLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_PARDISO
		PardisoLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_PETSC
		PetscLinearSolver rsolver(&mat, &rhs);
#elif defined WITH_MUMPS
		MumpsSolver rsolver(&mat, &rhs);
#endif

		Space *rspace = space.dup(&rmesh);
		rspace->copy_orders(space, 1);

		LinearProblem rlp(&wf);
		rlp.set_space(rspace);

		int rndofs = rspace->assign_dofs();
		printf("  - Number of DOFs: %d\n", rndofs);

		printf("  - Assembling... "); fflush(stdout);
		assemble_timer.reset();
		assemble_timer.start();
		rlp.assemble(&mat, &rhs);
		assemble_timer.stop();
		printf("done in %s (%lf secs)\n", assemble_timer.get_human_time(), assemble_timer.get_seconds());

		printf("  - Solving... "); fflush(stdout);
		solve_timer.reset();
		solve_timer.start();
		bool rsolved = rsolver.solve();
		solve_timer.stop();
		if (rsolved)
			printf("done in %s (%lf secs)\n", solve_timer.get_human_time(), solve_timer.get_seconds());
		else {
			res = ERR_FAILURE;
			printf("failed\n");
			break;
		}
		Solution rsln(&rmesh);
		rsln.set_fe_solution(rspace, rsolver.get_solution());
#else
		ExactSolution rsln(&rmesh, exact_solution);
#endif

		{
			double h1_err_norm = h1_error(&sln, &rsln) * 100;
			printf("  - H1 error norm:      % le\n", h1_err_norm);
		}

		printf("Adaptivity:\n");
		H1Adapt hp(&space);
		hp.set_aniso(true);
		double tol = hp.calc_error(&sln, &rsln) * 100;
		printf("  - tolerance: "); fflush(stdout);
		printf("% lf\n", tol);
		if (tol < TOLERANCE || iter == 8) {
			ExactSolution ex_sln(&mesh, exact_solution);
			printf("\nDone\n");
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

			break;
		}

		Timer t("");
		printf("  - adapting... "); fflush(stdout);
		t.start();
		hp.adapt(THRESHOLD);
		t.stop();
		printf("done in %lf secs (refined %d element(s))\n", t.get_seconds(), hp.get_num_refined_elements());

		iter++;
	} while (!done && iter < 5);

#ifdef WITH_PETSC
	PetscFinalize();
#endif

	return res;
}
