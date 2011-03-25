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
//          Ex(y) = E_0 * cos(y*M_PI/h), where h is height of the waveguide ()
//     3) Newton boundary condition (impedance matching) on Gamma_4
//          dE/dn = j*beta*E
//
// The following parameters can be changed:

const int P_INIT = 6;                                  // Initial polynomial degree of all elements.
const int INIT_REF_NUM = 3;                            // Number of initial mesh refinements.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;       // Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
                                                       // SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.

// Boundary markers.
const std::string BDY_PERFECT = "1", BDY_LEFT = "2", BDY_IMPEDANCE = "3";

// Problem parameters.
const double epsr = 1.0;                    // Relative permittivity
const double eps0 = 8.85418782e-12;         // Permittivity of vacuum F/m
const double eps = epsr * eps0;
const double mur = 1.0;                     // Relative permeablity
const double mu0 = 4*M_PI*1e-7;             // Permeability of vacuum H/m
const double mu = mur * mu0;
const double frequency = 3e9;               // Frequency MHz
const double omega = 2*M_PI * frequency;    // Angular velocity
const double sigma = 0;                     // Conductivity Ohm/m
const double beta = 54;                     // Propagation constant
const double E0 = 100;                        // Input electric intensity
const double h = 0.1;                       // Height of waveguide

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
    // Load the mesh.
    Mesh mesh;
    H2DReader mloader;
    mloader.load("../domain.mesh", &mesh);

    // Perform uniform mesh refinement.
    for(int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements(2); // 2 is for vertical split.

    // Initialize boundary conditions
    // Initialize boundary conditions
    DefaultEssentialBCConst bc1(BDY_PERFECT, 0.0);
    EssentialBCNonConst bc2(BDY_LEFT);
    //NaturalBoundaryCondition bc3(BDY_IMPEDANCE);
    EssentialBCs bcs(Hermes::vector<EssentialBC *>(&bc1, &bc2));

    // Create an H1 space with default shapeset.
    H1Space e_r_space(&mesh, &bcs, P_INIT);
    H1Space e_i_space(&mesh, &bcs, P_INIT);
    int ndof = Space::get_num_dofs(&e_r_space);
    info("ndof = %d", ndof);

    // Initialize the weak formulation
    // Weak forms for real and imaginary parts

    // Initialize the weak formulation.
    WeakFormHelmholtz wf(eps, mu, omega, sigma, beta, E0, h);

    // Initialize the solutions.
    Solution e_r_sln, e_i_sln;

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
        if (p_init == 2 && fabs(sum - 65.582) > 1e-3) success = 0;
        if (p_init == 3 && fabs(sum + 48.9119) > 1e-3) success = 0;
        if (p_init == 4 && fabs(sum + 46.299) > 1e-3) success = 0;
        if (p_init == 5 && fabs(sum + 39.8476) > 1e-3) success = 0;
        if (p_init == 6 && fabs(sum + 40.5802) > 1e-3) success = 0;
        if (p_init == 7 && fabs(sum + 40.8093) > 1e-3) success = 0;
        if (p_init == 8 && fabs(sum + 40.7958) > 1e-3) success = 0;
        if (p_init == 9 && fabs(sum + 40.7928) > 1e-3) success = 0;
        if (p_init == 10 && fabs(sum + 40.7929) > 1e-3) success = 0;
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
