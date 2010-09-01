#include "hermes1d.h"

// ********************************************************************
// This example shows solution of harmonic steady state on the homoegenous transmision line.
// Wave propagation is described by two linear differencial equation with complex coefficients.
// dU(x)/dx = -(R+j\omega L)I(x)
// dI(x)/dx = -(G+j\omega C)U(x)
// These equations are rewrited into four equations with real coefficients
// Ur' - R*Ir + \omega*L*Ii = 0
// Ui' - R*Ii - \omega*L*Ir = 0
// Ir' - G*Ur + \omega*C*Ui = 0
// Ii' - G*Ui - \omega*C*Ur = 0

// in an interval (0, 10) equipped with Dirichlet bdy conditions
// x1(0) = 1, x2(0) = 0, x3(0) = 0, x4(0) = 0

double L=25e-9;               // induktance [H/m]
double C=10e-12;              // capacitance [F/m]
double G=1e-9;                // conductance [S/m]
double R=1e-3;                // resistance [Ohm/m]
double l=10 ;                  // length of the line [m]
double omega=2*M_PI*3e8;        //
double Zl=60;                  // load impedance[Ohm]

//double Val_newton_alpha_U_Re=-R/Zl;
//double Val_newton_alpha_U_Im=-omega*L/Zl;
//double Val_newton_alpha_I_Re=-G*Zl;
//double Val_newton_alpha_I_Im=-omega*C*Zl;

// General input:
static int N_eq = 4;
int N_elem = 1000;          // number of elements
double A = 0, B = l;        // domain end points
int P_init = 2;             // initial polynomal degree

// Matrix solver
const int MATRIX_SOLVER = 1;            // 0... default (LU decomposition)
                                        // 1... UMFPACK
                                        // 2... CG (no preconditioning)
                                        // Only relevant for iterative matrix solvers:
const double MATRIX_SOLVER_TOL = 1e-7;  // Tolerance for residual in L2 norm
const int MATRIX_SOLVER_MAXITER = 150;  // Max. number of iterations

// Newton's method
double NEWTON_TOL = 1e-2;
int NEWTON_MAXITER = 150;

// Boundary conditions
double Val_dir_left_1 = 1;  // real part of the voltage at the beginnig of the line
double Val_dir_left_2 = 0;  // imaginary part of the voltage at the beginnig of the line
double Val_dir_left_3 = 0;  // real part of the voltage at the beginnig of the line
double Val_dir_left_4 = 0;  // imaginary part of the voltage at the beginnig of the line

//At the end of the line is an indirect boundary condition U(l) = I(l)*Zl see below

// Weak forms for Jacobi matrix and residual
#include "forms.cpp"

/******************************************************************************/
int main() {
    // create mesh
    Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
    mesh->set_bc_left_dirichlet(0, Val_dir_left_1);
    mesh->set_bc_left_dirichlet(1, Val_dir_left_2);
    mesh->set_bc_left_dirichlet(2, Val_dir_left_3);
    mesh->set_bc_left_dirichlet(3, Val_dir_left_4);
    printf("N_dof = %d\n", mesh->assign_dofs());

    // register weak forms
    DiscreteProblem *dp = new DiscreteProblem();
    dp->add_matrix_form(0, 0, jacobian_1_1);
    dp->add_matrix_form(0, 2, jacobian_1_3);
    dp->add_matrix_form(0, 3, jacobian_1_4);
    dp->add_matrix_form(1, 1, jacobian_2_2);
    dp->add_matrix_form(1, 2, jacobian_2_3);
    dp->add_matrix_form(1, 3, jacobian_2_4);
    dp->add_matrix_form(2, 0, jacobian_3_1);
    dp->add_matrix_form(2, 1, jacobian_3_2);
    dp->add_matrix_form(2, 2, jacobian_3_3);
    dp->add_matrix_form(3, 0, jacobian_4_1);
    dp->add_matrix_form(3, 1, jacobian_4_2);
    dp->add_matrix_form(3, 3, jacobian_4_4);
    dp->add_vector_form(0, residual_1);
    dp->add_vector_form(1, residual_2);
    dp->add_vector_form(2, residual_3);
    dp->add_vector_form(3, residual_4);
    dp->add_matrix_form_surf(0, 0, jacobian_surf_right_U_Re, BOUNDARY_RIGHT);
    dp->add_matrix_form_surf(0, 2, jacobian_surf_right_U_Im, BOUNDARY_RIGHT);
    dp->add_matrix_form_surf(1, 1, jacobian_surf_right_I_Re, BOUNDARY_RIGHT);
    dp->add_matrix_form_surf(1, 3, jacobian_surf_right_I_Im, BOUNDARY_RIGHT);

    // Newton's loop
    newton(dp, mesh, MATRIX_SOLVER, MATRIX_SOLVER_TOL, MATRIX_SOLVER_MAXITER,
           NEWTON_TOL, NEWTON_MAXITER);

    // Plot the solution
    Linearizer l(mesh);
    l.plot_solution("solution.gp");

    printf("Done.\n");
    return 1;
}
