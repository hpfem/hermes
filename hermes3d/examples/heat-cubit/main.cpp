#include "config.h"
#ifdef WITH_PETSC
#include <petsc.h>
#endif
#include <getopt.h>
#include <hermes3d.h>

// Solving a simple heat equation to demonstrate how to use CUBIT with Hermes3D
//
// Use mesh file from 'meshes/exodusII/cylynder2.e. Material IDs corresponds
// to elements markers, sideset IDs correspond to face (BC) markers
//

// commnad line arguments
bool do_output = true;				// generate output files (if true)
char *mesh_file_name = NULL;		// the name of the mesh file

// usage info

void usage() {
	printf("Usage\n");
	printf("\n");
	printf("  heat-cubit [options] <mesh-file>\n");
	printf("\n");
	printf("Options:\n");
	printf("  --no-output         - do not generate output files\n");
	printf("\n");
}

bool process_cmd_line(int argc, char **argv)
{
	static struct option long_options[] = {
		{ "no-output", no_argument, (int *) &do_output, false },
		{ 0, 0, 0, 0 }
	};

	// getopt_long stores the option index here.
	int option_index = 0;
	int c;
	while ((c = getopt_long(argc, argv, "", long_options, &option_index)) != -1) {
		switch (c) {
            case 0:
				break;

			case '?':
				// getopt_long already printed an error message
				break;

			default:
				return false;
		}
	}

	if (optind < argc) {
		mesh_file_name = argv[optind++];
		return true;
	}
	else
		return false;
}

BCType bc_types(int marker)
{
	if (marker == 1) return BC_ESSENTIAL;
	else return BC_NATURAL;
}

scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
	return 10;
}

template<typename f_t, typename res_t>
res_t bilinear_form1(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                    user_data_t<res_t> *data)
{
	return 10 * int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename f_t, typename res_t>
res_t bilinear_form2(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, fn_t<f_t> *v, geom_t<f_t> *e,
                    user_data_t<res_t> *data)
{
	return 0.5 * int_grad_u_grad_v<f_t, res_t>(n, wt, u, v, e);
}

template<typename T>
T rhs(T x, T y, T z)
{
	return 4.0;
}

template<typename f_t, typename res_t>
res_t linear_form(int n, double *wt, fn_t<res_t> *u_ext[], fn_t<f_t> *u, geom_t<f_t> *e, user_data_t<res_t> *data)
{
	return -int_F_v<f_t, res_t>(n, wt, rhs, u, e);
}

//

void out_fn(MeshFunction *x, const char *name)
{
	char of_name[1024];
	FILE *ofile;
	// mesh out
	sprintf(of_name, "%s.vtk", name);
	ofile = fopen(of_name, "w");
	if (ofile != NULL) {
		VtkOutputEngine output(ofile);
		output.out(x, name);
		fclose(ofile);
	}
	else {
		warning("Can not open '%s' for writing.", of_name);
	}
}

void out_bc(Mesh *mesh, const char *name)
{
	char of_name[1024];
	FILE *ofile;
	// mesh out
	sprintf(of_name, "%s.vtk", name);
	ofile = fopen(of_name, "w");
	if (ofile != NULL) {
		VtkOutputEngine output(ofile);
		output.out_bc(mesh, name);
		fclose(ofile);
	}
	else {
		warning("Can not open '%s' for writing.", of_name);
	}
}

// main ///////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
#ifdef WITH_PETSC
	PetscInitialize(&argc, &argv, (char *) PETSC_NULL, PETSC_NULL);
#endif

	if (!process_cmd_line(argc, argv)) {
		usage();
		return 0;
	}

	printf("* Loading mesh '%s'\n", mesh_file_name);
	Mesh mesh;
	ExodusIIReader mesh_loader;
	if (!mesh_loader.load(mesh_file_name, &mesh))
		error("Loading mesh file '%s'\n", mesh_file_name);

	H1ShapesetLobattoHex shapeset;

	printf("* Setting the space up\n");
	H1Space space(&mesh, &shapeset);
	space.set_bc_types(bc_types);
	space.set_essential_bc_values(essential_bc_values);
	space.set_uniform_order(order3_t(1, 1, 1));

	int ndofs = space.assign_dofs();
	printf("  - Number of DOFs: %d\n", ndofs);

	WeakForm wf;
	wf.add_matrix_form(bilinear_form1<double, scalar>, bilinear_form1<ord_t, ord_t>, SYM, 1);
	wf.add_matrix_form(bilinear_form2<double, scalar>, bilinear_form2<ord_t, ord_t>, SYM, 2);
	wf.add_vector_form(linear_form<double, scalar>, linear_form<ord_t, ord_t>, ANY);

	LinearProblem lp(&wf);
	lp.set_space(&space);

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

	// assemble stiffness matrix
	printf("  - assembling... "); fflush(stdout);
	Timer tmr_assemble;
	tmr_assemble.start();
	bool assembled = lp.assemble(&mat, &rhs);
	tmr_assemble.stop();
	if (assembled)
		printf("done in %s (%lf secs)\n", tmr_assemble.get_human_time(), tmr_assemble.get_seconds());
	else
		error("failed!");

	// solve the stiffness matrix
	printf("  - solving... "); fflush(stdout);
	Timer tmr_solve;
	tmr_solve.start();
	bool solved = solver.solve();
	tmr_solve.stop();

	if (solved) {
		printf("done in %s (%lf secs)\n", tmr_solve.get_human_time(), tmr_solve.get_seconds());
		double *s = solver.get_solution();

		Solution sln(&mesh);
		sln.set_fe_solution(&space, s);

		if (do_output) {
			printf("  - output... "); fflush(stdout);
			out_bc(&mesh, "bc");
			out_fn(&sln, "temp");
			printf("done\n");
		}
	}
	else {
		printf("failed\n");
	}

#ifdef WITH_PETSC
	mat.free();
	rhs.free();
	PetscFinalize();
#endif

	return 0;
}
