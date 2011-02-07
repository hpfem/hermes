#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

// This example shows how to model harmonic steady state in parallel plate waveguide.
// The Helmholtz equation is solved and there are demonstrated two typical boundary
// conditions used in high frequency domain.
//
// PDE: Helmholtz equation for electric field
//
//    Delta E  + (omega^2*mu*epsilon - j*omega*sigma*mu)*E = 0
//
// BC:              Gamma_1
//             ----------------------------
//  Gamma_3    |                           |  Gamma_4
//             ----------------------------
//                  Gamma_2
//
//     1) Dirichlet boundary condition Ex = 0 (perfect eletric conductor) on Gamma_1 and Gamma_2.
//     2) Essential (Dirichlet) boundary condition on Gamma_3
//          Ex(y) = E_0 * cos(y*M_PI/h), where h is height of the waveguide
//     3) Newton boundary condition (impedance matching) on Gamma_4
//          dE/dn = j*beta*E
//
// The following parameters can be changed:

const int P_INIT = 6;                                  // Initial polynomial degree of all elements.
const int INIT_REF_NUM = 3;                            // Number of initial mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;       // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PARDISO, SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const int BDY_PERFECT = 1, BDY_LEFT = 2, BDY_IMPEDANCE = 3;

// Problem parameters.
const double epsr = 1.0;                    // Relative permittivity
const double eps0 = 8.85418782e-12;         // Permittivity of vacuum F/m
const double mur = 1.0;                     // Relative permeablity
const double mu0 = 4*M_PI*1e-7;             // Permeability of vacuum H/m
const double frequency = 3e9;               // Frequency MHz
const double omega = 2*M_PI * frequency;    // Angular velocity
const double sigma = 0;                     // Conductivity Ohm/m
const double beta = 54;                     // Propagation constant
const double E0 = 100;                      // Input electric intensity
const double h = 0.1;                       // Height of waveguide


// Distribution of electric field for TE1 mode
scalar essential_bc_values(double x, double y)
{    
    return E0*cos(y*M_PI/h); //
}

// Weak forms.
#include "forms.cpp"

int main(int argc, char* argv[])
{
    // Load the mesh.
    Mesh mesh;
    H2DReader mloader;
    mloader.load("domain.mesh", &mesh);

    // Perform uniform mesh refinement.
    for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(2); // 2 is for vertical split.

    // Enter boundary markers.
    BCTypes bc_types;    
    bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_LEFT, BDY_PERFECT));
    bc_types.add_bc_newton(Hermes::vector<int>(BDY_IMPEDANCE));

    // Enter Dirichlet boundary values;
    BCValues bc_values_r;
    bc_values_r.add_const(BDY_PERFECT, 0.0);
    bc_values_r.add_function(BDY_LEFT, essential_bc_values);

    BCValues bc_values_i;
    bc_values_i.add_const(BDY_LEFT, 0.0);
    bc_values_i.add_const(BDY_PERFECT, 0.0);

    // Create H1 shapeset.
    H1Space e_r_space(&mesh, &bc_types, &bc_values_r, P_INIT);
    H1Space e_i_space(&mesh, &bc_types, &bc_values_i, P_INIT);
    info("ndof = %d.", Space::get_num_dofs(Hermes::vector<Space *>(&e_r_space, &e_i_space)));

    // Initialize the weak formulation
    // Weak forms for real and imaginary parts

    WeakForm wf(2);
    wf.add_matrix_form(0, 0, callback(matrix_form_real_real));
    wf.add_matrix_form(0, 1, callback(matrix_form_real_imag));
    wf.add_matrix_form(1, 1, callback(matrix_form_imag_imag));
    wf.add_matrix_form(1, 0, callback(matrix_form_imag_real));

    // Impedance matching - Newton boundary condition
    wf.add_matrix_form_surf(0, 1, callback(matrix_form_surface_imag_real), BDY_IMPEDANCE);
    wf.add_matrix_form_surf(1, 0, callback(matrix_form_surface_real_imag), BDY_IMPEDANCE);

    int success = 1;
    for (int p_init = 2; p_init <= 10; p_init++) {

        printf("********* p_init = %d *********\n", p_init);
        e_r_space.set_uniform_order(p_init);
        e_i_space.set_uniform_order(p_init);

        // Initialize the FE problem.
        bool is_linear = true;
        DiscreteProblem dp(&wf, Hermes::vector<Space *>(&e_r_space, &e_i_space), is_linear);

        // Set up the solver, matrix, and rhs according to the solver selection.
        SparseMatrix* matrix = create_matrix(matrix_solver);
        Vector* rhs = create_vector(matrix_solver);
        Solver* solver = create_linear_solver(matrix_solver, matrix, rhs);

        // Initialize the solutions.
        Solution e_r_sln, e_i_sln;

        // Assemble the stiffness matrix and right-hand side vector.
        info("Assembling the stiffness matrix and right-hand side vector.");
        bool rhsonly = false;
        dp.assemble(matrix, rhs, rhsonly);

        // Solve the linear system and if successful, obtain the solutions.
        info("Solving the matrix problem.");
        if(solver->solve())
            Solution::vector_to_solutions(solver->get_solution(),
                                          Hermes::vector<Space *>(&e_r_space, &e_i_space),
                                          Hermes::vector<Solution *>(&e_r_sln, &e_i_sln));
        else
            error ("Matrix solver failed.\n");

        int ndof = Space::get_num_dofs(&e_r_space);
        printf("ndof = %d\n", ndof);
        double sum = 0;
        for (int i=0; i < ndof; i++) sum += solver->get_solution()[i];
        printf("coefficient sum = %g\n", sum);

        // Actual test. The values of 'sum' depend on the
        // current shapeset. If you change the shapes  et,
        // you need to correct these numbers.
        if (p_init == 2 && fabs(sum - 29.6682) > 1e-3) success = 0;
        if (p_init == 3 && fabs(sum - 153.823) > 1e-3) success = 0;
        if (p_init == 4 && fabs(sum - 112.541) > 1e-3) success = 0;
        if (p_init == 5 && fabs(sum - 117.415) > 1e-3) success = 0;
        if (p_init == 6 && fabs(sum - 118.358) > 1e-3) success = 0;
        if (p_init == 7 && fabs(sum - 118.213) > 1e-3) success = 0;
        if (p_init == 8 && fabs(sum - 118.198) > 1e-3) success = 0;
        if (p_init == 9 && fabs(sum - 118.2) > 1e-3) success = 0;
        if (p_init == 10 && fabs(sum - 118.2) > 1e-3) success = 0;
    }

    if (success == 1) {
        printf("Success!\n");
        return ERR_SUCCESS;
    }
    else {
        printf("Failure!\n");
        return ERR_FAILURE;
    }
}
