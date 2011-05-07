#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include "config.h"
//#include <getopt.h>
#include <hermes3d.h>

// Usage: $0 <mesh file> <element id> <refinement id> [<element id> <refinement id>...].

//#define DIRICHLET
#define NEWTON

// The following parameters can be changed:
const char* iterative_method = "bicgstab";        // Name of the iterative method employed by AztecOO (ignored
                                                  // by the other solvers). 
                                                  // Possibilities: gmres, cg, cgs, tfqmr, bicgstab.
const char* preconditioner = "jacobi";            // Name of the preconditioner employed by AztecOO (ignored by
                                                  // the other solvers). 
                                                  // Possibilities: none, jacobi, neumann, least-squares, or a
                                                  // preconditioner from IFPACK (see solver/aztecoo.h).
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

//#define X2_Y2_Z2
//#define XM_YN_ZO
#define XM_YN_ZO_2


// Everything is tested on the following geometry
//
//
//      7             6             11
//        +-----------+-----------+
//       /|          /|          /|
//      / |         / |         / |
//     /  |      5 /  |     10 /  |
//  4 +-----------+-----------+   |
//    |   |       |   |       |   |
//    |   +-------|---+-------|---+
//    |  / 3      |  / 2      |  / 9
//    | /         | /         | /
//    |/          |/          |/
//    +-----------+-----------+
//   0            1            8
//
//  ^ z
//  |
//  | / y
//  |/
//  +---> x

// vertices
Point3D vtcs[] = {
	{ -2, -1, -1 },
	{  0, -1, -1 },
	{  0,  1, -1 },
	{ -2,  1, -1 },
	{ -2, -1,  1 },
	{  0, -1,  1 },
	{  0,  1,  1 },
	{ -2,  1,  1 },
	{  2, -1, -1 },
	{  2,  1, -1 },
	{  2, -1,  1 },
	{  2,  1,  1 },
};

// mesh
unsigned int hexs[2][48][8] = {
	// hexs 1
	{
		{  0,  1,  2,  3,  4,  5,  6,  7 },
		{  1,  2,  3,  0,  5,  6,  7,  4 },
		{  2,  3,  0,  1,  6,  7,  4,  5 },
		{  3,  0,  1,  2,  7,  4,  5,  6 },

		{  4,  7,  6,  5,  0,  3,  2,  1 },
		{  5,  4,  7,  6,  1,  0,  3,  2 },
		{  6,  5,  4,  7,  2,  1,  0,  3 },
		{  7,  6,  5,  4,  3,  2,  1,  0 },

		{  1,  5,  6,  2,  0,  4,  7,  3 },
		{  5,  6,  2,  1,  4,  7,  3,  0 },
		{  6,  2,  1,  5,  7,  3,  0,  4 },
		{  2,  1,  5,  6,  3,  0,  4,  7 },

		{  0,  3,  7,  4,  1,  2,  6,  5 },
		{  4,  0,  3,  7,  5,  1,  2,  6 },
		{  7,  4,  0,  3,  6,  5,  1,  2 },
		{  3,  7,  4,  0,  2,  6,  5,  1 },

		{  1,  0,  4,  5,  2,  3,  7,  6 },
		{  0,  4,  5,  1,  3,  7,  6,  2 },
		{  4,  5,  1,  0,  7,  6,  2,  3 },
		{  5,  1,  0,  4,  6,  2,  3,  7 },

		{  2,  6,  7,  3,  1,  5,  4,  0 },
		{  3,  2,  6,  7,  0,  1,  5,  4 },
		{  7,  3,  2,  6,  4,  0,  1,  5 },
		{  6,  7,  3,  2,  5,  4,  0,  1 },

		{  0,  3,  7,  4,  1,  2,  6,  5 },
		{  4,  0,  3,  7,  5,  1,  2,  6 },
		{  7,  4,  0,  3,  6,  5,  1,  2 },
		{  3,  7,  4,  0,  2,  6,  5,  1 },

		{  1,  5,  6,  2,  0,  4,  7,  3 },
		{  5,  6,  2,  1,  4,  7,  3,  0 },
		{  6,  2,  1,  5,  7,  3,  0,  4 },
		{  2,  1,  5,  6,  3,  0,  4,  7 },

		{  7,  6,  5,  4,  3,  2,  1,  0 },
		{  4,  7,  6,  5,  0,  3,  2,  1 },
		{  5,  4,  7,  6,  1,  0,  3,  2 },
		{  6,  5,  4,  7,  2,  1,  0,  3 },

		{  3,  0,  1,  2,  7,  4,  5,  6 },
		{  0,  1,  2,  3,  4,  5,  6,  7 },
		{  1,  2,  3,  0,  5,  6,  7,  4 },
		{  2,  3,  0,  1,  6,  7,  4,  5 },

		{  2,  6,  7,  3,  1,  5,  4,  0 },
		{  6,  7,  3,  2,  5,  4,  0,  1 },
		{  7,  3,  2,  6,  4,  0,  1,  5 },
		{  3,  2,  6,  7,  0,  1,  5,  4 },

		{  1,  0,  4,  5,  2,  3,  7,  6 },
		{  5,  1,  0,  4,  6,  2,  3,  7 },
		{  4,  5,  1,  0,  7,  6,  2,  3 },
		{  0,  4,  5,  1,  3,  7,  6,  2 }
	},
	// hexs 2
	{
		{  1,  8,  9,  2,  5, 10, 11,  6 },
		{  8,  9,  2,  1, 10, 11,  6,  5 },
		{  9,  2,  1,  8, 11,  6,  5, 10 },
		{  2,  1,  8,  9,  6,  5, 10, 11 },

		{  5,  6, 11, 10,  1,  2,  9,  8 },
		{ 10,  5,  6, 11,  8,  1,  2,  9 },
		{ 11, 10,  5,  6,  9,  8,  1,  2 },
		{  6, 11, 10,  5,  2,  9,  8,  1 },

		{  8, 10, 11,  9,  1,  5,  6,  2 },
		{ 10, 11,  9,  8,  5,  6,  2,  1 },
		{ 11,  9,  8, 10,  6,  2,  1,  5 },
		{  9,  8, 10, 11,  2,  1,  5,  6 },

		{  1,  2,  6,  5,  8,  9, 11, 10 },
		{  5,  1,  2,  6, 10,  8,  9, 11 },
		{  6,  5,  1,  2, 11, 10,  8,  9 },
		{  2,  6,  5,  1,  9, 11, 10,  8 },

		{  8,  1,  5, 10,  9,  2,  6, 11 },
		{  1,  5, 10,  8,  2,  6, 11,  9 },
		{  5, 10,  8,  1,  6, 11,  9,  2 },
		{ 10,  8,  1,  5, 11,  9,  2,  6 },

		{  9, 11,  6,  2,  8, 10,  5,  1 },
		{  2,  9, 11,  6,  1,  8, 10,  5 },
		{  6,  2,  9, 11,  5,  1,  8, 10 },
		{ 11,  6,  2,  9, 10,  5,  1,  8 },

		{  1,  2,  6,  5,  8,  9, 11, 10 },
		{  5,  1,  2,  6, 10,  8,  9, 11 },
		{  6,  5,  1,  2, 11, 10,  8,  9 },
		{  2,  6,  5,  1,  9, 11, 10,  8 },

		{  8, 10, 11,  9,  1,  5,  6,  2 },
		{ 10, 11,  9,  8,  5,  6,  2,  1 },
		{ 11,  9,  8, 10,  6,  2,  1,  5 },
		{  9,  8, 10, 11,  2,  1,  5,  6 },

		{  6, 11, 10,  5,  2,  9,  8,  1 },
		{  5,  6, 11, 10,  1,  2,  9,  8 },
		{ 10,  5,  6, 11,  8,  1,  2,  9 },
		{ 11, 10,  5,  6,  9,  8,  1,  2 },

		{  2,  1,  8,  9,  6,  5, 10, 11 },
		{  1,  8,  9,  2,  5, 10, 11,  6 },
		{  8,  9,  2,  1, 10, 11,  6,  5 },
		{  9,  2,  1,  8, 11,  6,  5,  10 },

		{  9, 11,  6,  2,  8, 10,  5,  1 },
		{ 11,  6,  2,  9, 10,  5,  1,  8 },
		{  6,  2,  9, 11,  5,  1,  8, 10 },
		{  2,  9, 11,  6,  1,  8, 10,  5 },

		{  8,  1,  5, 10,  9,  2,  6, 11 },
		{ 10,  8,  1,  5, 11,  9,  2,  6 },
		{  5, 10,  8,  1,  6, 11,  9,  2 },
		{  1,  5, 10,  8,  2,  6, 11,  9 }
	}
};

unsigned int bnd[10][5] = {
	{ 0,  3,  7,  4, 1 },
	{ 8,  9, 11, 10, 2 },
	{ 0,  1,  5,  4, 3 },
	{ 1,  8, 10,  5, 3 },
	{ 3,  2,  6,  7, 4 },
	{ 2,  9, 11,  6, 4 },
	{ 0,  1,  2,  3, 5 },
	{ 1,  8,  9,  2, 5 },
	{ 4,  5,  6,  7, 6 },
	{ 5, 10, 11,  6, 6 }
};


int m = 2, n = 2, o = 2;

template<typename T>
T fnc(T x, T y, T z) {
#ifdef XM_YN_ZO
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 3) * z + pow(z, 4);
#elif defined XM_YN_ZO_2
	return pow(x, m) * pow(y, n) * pow(z, o) + pow(x, 2) * pow(y, 3) - pow(x, 2) * z + pow(z, 4);
#elif defined X2_Y2_Z2
	return x*x + y*y + z*z;
#endif
}

template<typename T>
T dfnc(T x, T y, T z) {
#ifdef XM_YN_ZO
	T ddxx = m * (m - 1) * pow(x, m - 2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 6 * x * z;
	T ddyy = n * (n - 1) * pow(x, m) * pow(y, n - 2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = o * (o - 1) * pow(x, m) * pow(y, n) * pow(z, o - 2) + 12 * pow(z, 2);
	return -(ddxx + ddyy + ddzz);

#elif defined XM_YN_ZO_2
	T ddxx = m*(m-1) * pow(x, m-2) * pow(y, n) * pow(z, o) + 2 * pow(y, 3) - 2 * z;
	T ddyy = n*(n-1) * pow(x, m) * pow(y, n-2) * pow(z, o) + 6 * pow(x, 2) * y;
	T ddzz = o*(o-1) * pow(x, m) * pow(y, n) * pow(z, o-2) + 12 * pow(z, 2);
	return -(ddxx + ddyy + ddzz);
#elif defined X2_Y2_Z2
	return -6.0;
#endif
}

// needed for calculation norms and used by visualizator
double exact_solution(double x, double y, double z, double &dx, double &dy, double &dz) {
#ifdef XM_YN_ZO
	dx = m * pow(x, m - 1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 3 * pow(x, 2) * z;
	dy = n * pow(x, m) * pow(y, n - 1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o - 1) - pow(x, 3) + 4 * pow(z, 3);
#elif defined XM_YN_ZO_2
	dx = m * pow(x, m-1) * pow(y, n) * pow(z, o) + 2 * x * pow(y, 3) - 2 * x * z;
	dy = n * pow(x, m) * pow(y, n-1) * pow(z, o) + 3 * pow(x, 2) * pow(y, 2);
	dz = o * pow(x, m) * pow(y, n) * pow(z, o-1) - pow(x, 2) + 4 * pow(z, 3);
#elif defined X2_Y2_Z2
	dx = 2 * x;
	dy = 2 * y;
	dz = 2 * z;
#endif

	return fnc(x, y, z);
}

// Boundary condition types.
BCType bc_types(int marker) {
#ifdef DIRICHLET
	return H3D_BC_ESSENTIAL;
#elif defined NEWTON
	return H3D_BC_NATURAL;
#endif
}

// Dirichlet boundary conditions.
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z) {
#ifdef DIRICHLET
	return fnc(x, y, z);
#else
	return 0;
#endif
}

template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) {
	return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar bilinear_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *data) {
	return int_u_v<Real, Scalar>(n, wt, u, v, e);
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *data) {
	return int_F_v<Real, Scalar>(n, wt, dfnc, u, e);
}

template<typename Real, typename Scalar>
Scalar linear_form_surf(int np, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *data) {
	Scalar result = 0;
	for (int i = 0; i < np; i++) {
#ifdef XM_YN_ZO
		Scalar dx = m * pow(e->x[i], m - 1) * pow(e->y[i], n) * pow(e->z[i], o) + 2 * e->x[i] * pow(e->y[i], 3) - 3 * pow(e->x[i], 2) * e->z[i];
		Scalar dy = n * pow(e->x[i], m) * pow(e->y[i], n - 1) * pow(e->z[i], o) + 3 * pow(e->x[i], 2) * pow(e->y[i], 2);
		Scalar dz = o * pow(e->x[i], m) * pow(e->y[i], n) * pow(e->z[i], o - 1) - pow(e->x[i], 3) + 4 * pow(e->z[i], 3);
#elif defined XM_YN_ZO_2
		Scalar dx = m * pow(e->x[i], m-1) * pow(e->y[i], n) * pow(e->z[i], o) + 2 * e->x[i] * pow(e->y[i], 3) - 2 * e->x[i] * e->z[i];
		Scalar dy = n * pow(e->x[i], m) * pow(e->y[i], n-1) * pow(e->z[i], o) + 3 * pow(e->x[i], 2) * pow(e->y[i], 2);
		Scalar dz = o * pow(e->x[i], m) * pow(e->y[i], n) * pow(e->z[i], o-1) - pow(e->x[i], 2) + 4 * pow(e->z[i], 3);
#elif defined X2_Y2_Z2
		Scalar dx = 2 * e->x[i];
		Scalar dy = 2 * e->y[i];
		Scalar dz = 2 * e->z[i];
#endif
		result += wt[i] * (u->val[i] * (dx * e->nx[i] + dy * e->ny[i] + dz * e->nz[i] + fnc(e->x[i], e->y[i], e->z[i])));
	}
	return result;
}

// helpers ////////////////////////////////////////////////////////////////////////////////////////

int parse_reft(char *str) {
	if (strcasecmp(str, "x") == 0) return H3D_REFT_HEX_X;
	else if (strcasecmp(str, "y") == 0) return H3D_REFT_HEX_Y;
	else if (strcasecmp(str, "z") == 0) return H3D_REFT_HEX_Z;
	else if (strcasecmp(str, "xy") == 0 || strcasecmp(str, "yx") == 0) return H3D_H3D_REFT_HEX_XY;
	else if (strcasecmp(str, "xz") == 0 || strcasecmp(str, "zx") == 0) return H3D_H3D_REFT_HEX_XZ;
	else if (strcasecmp(str, "yz") == 0 || strcasecmp(str, "zy") == 0) return H3D_H3D_REFT_HEX_YZ;
	else if (strcasecmp(str, "xyz") == 0) return H3D_H3D_H3D_REFT_HEX_XYZ;
	else return H3D_REFT_HEX_NONE;
}

const double EPS = 10e-12;

int main(int argc, char **args) 
{
  // Test variable.
  int success_test = 1;

	for (int i = 0; i < 48; i++) {
		for (int j = 0; j < 48; j++) {
			info("Config: %d, %d ", i, j);

			Mesh mesh;

			for (unsigned int k = 0; k < countof(vtcs); k++)
				mesh.add_vertex(vtcs[k].x, vtcs[k].y, vtcs[k].z);
			unsigned int h1[] = {
					hexs[0][i][0] + 1, hexs[0][i][1] + 1, hexs[0][i][2] + 1, hexs[0][i][3] + 1,
					hexs[0][i][4] + 1, hexs[0][i][5] + 1, hexs[0][i][6] + 1, hexs[0][i][7] + 1 };
			mesh.add_hex(h1);
			unsigned int h2[] = {
					hexs[1][j][0] + 1, hexs[1][j][1] + 1, hexs[1][j][2] + 1, hexs[1][j][3] + 1,
					hexs[1][j][4] + 1, hexs[1][j][5] + 1, hexs[1][j][6] + 1, hexs[1][j][7] + 1 };
			mesh.add_hex(h2);
			// bc
			for (unsigned int k = 0; k < countof(bnd); k++) {
				unsigned int facet_idxs[Quad::NUM_VERTICES] = { bnd[k][0] + 1, bnd[k][1] + 1, bnd[k][2] + 1, bnd[k][3] + 1 };
				mesh.add_quad_boundary(facet_idxs, bnd[k][4]);
			}

			mesh.ugh();

      // Initialize the space.
			H1Space space(&mesh, bc_types, essential_bc_values);
			
#ifdef XM_YN_ZO
			Ord3 ord(4, 4, 4);
#elif defined XM_YN_ZO_2
			Ord3 ord(4, 4, 4);
#elif defined X2_Y2_Z2
			Ord3 ord(2, 2, 2);
#endif
			space.set_uniform_order(ord);

      // Initialize the weak formulation.
      WeakForm wf;
#ifdef DIRICHLET
      wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
      wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>);
#elif defined NEWTON
      wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
      wf.add_matrix_form_surf(bilinear_form_surf<double, scalar>, bilinear_form_surf<Ord, Ord>);
      wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>);
      wf.add_vector_form_surf(linear_form_surf<double, scalar>, linear_form_surf<Ord, Ord>);
#endif

      // Initialize the FE problem.
      bool is_linear = true;
      DiscreteProblem dp(&wf, &space, is_linear);

      // Set up the solver, matrix, and rhs according to the solver selection.
      SparseMatrix* matrix = create_matrix(matrix_solver);
      Vector* rhs = create_vector(matrix_solver);
      Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);
      
      // Initialize the preconditioner in the case of SOLVER_AZTECOO.
      if (matrix_solver == SOLVER_AZTECOO) 
      {
        ((AztecOOSolver*) solver)->set_solver(iterative_method);
        ((AztecOOSolver*) solver)->set_precond(preconditioner);
        // Using default iteration parameters (see solver/aztecoo.h).
      }

      // Assemble the linear problem.
      info("Assembling (ndof: %d).", Space::get_num_dofs(&space));
      dp.assemble(matrix, rhs);
        
      // Solve the linear system. If successful, obtain the solution.
      info("Solving.");
      Solution sln(space.get_mesh());
      if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), &space, &sln);
      else error ("Matrix solver failed.\n");


      ExactSolution ex_sln(&mesh, exact_solution);

      // Calculate exact error.
      info("Calculating exact error.");
      Adapt *adaptivity = new Adapt(&space, HERMES_H1_NORM);
      bool solutions_for_adapt = false;
      double err_exact = adaptivity->calc_err_exact(&sln, &ex_sln, solutions_for_adapt, HERMES_TOTAL_ERROR_ABS);

      if (err_exact > EPS)
      {
        // Calculated solution is not precise enough.
	      success_test = 0;
        info("failed, error:%g", err_exact);
      }
      else
        info("passed");

      // Clean up.
      delete matrix;
      delete rhs;
      delete solver;
      delete adaptivity;
		}
	}
  
  if (success_test) {
    info("Success!");
    return ERR_SUCCESS;
  }
  else {
    info("Failure!");
    return ERR_FAILURE;
  }
}

