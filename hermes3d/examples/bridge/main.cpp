#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#include <hermes3d.h>

// This example shows the Bridge model
// a hexahedral mesh in CTU format.
//
// PDE: Laplace equation -Laplace u = f, where f = CONST_F.
//

// Right-hand side.
const double CONST_F = -1.0;  

// Boundary markers.
int bdy_supported = 1;
const int P_INIT_X = 2,
          P_INIT_Y = 2,
          P_INIT_Z = 2;                           // Initial polynomial degree of all mesh elements.

MatrixSolverType matrix_solver = SOLVER_UMFPACK;  // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                  // SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.                                               

// Boundary condition types. 
BCType bc_types(int marker) 
{
  return (marker == bdy_supported) ? BC_ESSENTIAL : BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values. 
scalar essential_bc_values(int ess_bdy_marker, double x, double y, double z)
{
  return 0;
}

#include "forms.cpp"

int main(int argc, char **args) 
{
	// Time measurement.
	TimePeriod cpu_time;
	cpu_time.tick();

    // Load the mesh. 
    Mesh mesh;

    //CTUReader mloader;
    H3DReader mloader;

    std::cout << "Loading mesh ...\n" ;
    //mloader.load("./data/most-sup.top", &mesh);
    mloader.load("./data/most.mesh3d", &mesh);

/*
    // Create H1 space with default shapeset.
    H1Space space(&mesh, bc_types, essential_bc_values, Ord3(P_INIT_X, P_INIT_Y, P_INIT_Z));

    // Construct initial solution and set it to zero.
    Solution sln_prev(&mesh);
    sln_prev.set_zero();

    // Initialize weak formulation. 
    WeakForm wf;
    wf.add_matrix_form(bilinear_form<double, scalar>, bilinear_form<Ord, Ord>, HERMES_SYM);
    wf.add_vector_form(linear_form<double, scalar>, linear_form<Ord, Ord>, HERMES_ANY, &sln_prev);

    // Initialize discrete problem.
    bool is_linear = true;
    DiscreteProblem dp(&wf, &space, is_linear);
*/

    // Set up the solver, matrix, and rhs according to the solver selection.
    SparseMatrix* matrix = create_matrix(matrix_solver);
    Vector* rhs = create_vector(matrix_solver);
    Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

    // Time measurement.
    cpu_time.tick();

    // Print timing information.
    info("Solutions and mesh with polynomial orders saved. Total running time: %g s", cpu_time.accumulated());

    // Clean up.
    delete matrix;
    delete rhs;
    delete solver;

    return 0;
}
