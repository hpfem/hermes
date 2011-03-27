#include "forms.h"

/* We solve A*x = E*B*x, where A, B are matrices composed of 2x2 blocks */

static double kappa = 1;
static double Z = 1.0;

// speed of light (in atomic units)
static double c = 137.03599911;
// Hydrogen atom:
#define _V(r) (-Z/(r))
// Harmonic oscillator
//#define _V(r) ((r)*(r))

double A_00(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        double r = x[i];
        // specify V:
        double V = _V(r);
        val += u[i]*V*v[i]*weights[i];
    }
    return val;
}

double A_11(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        double r = x[i];
        // specify V:
        double V = _V(r);
        val += u[i]*(V-2*c*c)*v[i]*weights[i];
    }
    return val;
}

double A_01(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        double r = x[i];
        val += c*(-u[i]*dvdx[i] + u[i]*(kappa/r)*v[i])*weights[i];
    }
    return val;
}

double A_10(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        double r = x[i];
        val += c*(u[i]*dvdx[i] + u[i]*(kappa/r)*v[i])*weights[i];
    }
    return val;
}

double B_00(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for(int i = 0; i < num; i++) {
        double r = x[i];
        val += u[i]*v[i]*weights[i];
    }
    return val;
}

void assemble_dirac(Mesh *mesh, Matrix *A, Matrix *B,
        int _kappa, int _Z)
{
    kappa = _kappa;
    Z = _Z;
    mesh->set_bc_left_dirichlet(0, 0);
    mesh->set_bc_right_dirichlet(0, 0);

    // variable for the total number of DOF
    int N_dof = mesh->assign_dofs();

    // register weak forms
    DiscreteProblem *dp1 = new DiscreteProblem();
    dp1->add_matrix_form(0, 0, A_00);
    dp1->add_matrix_form(0, 1, A_01);
    dp1->add_matrix_form(1, 0, A_10);
    dp1->add_matrix_form(1, 1, A_11);
    DiscreteProblem *dp2 = new DiscreteProblem();
    dp2->add_matrix_form(0, 0, B_00);
    dp2->add_matrix_form(1, 1, B_00);

    printf("Assembling A, B. ndofs: %d\n", N_dof);
    dp1->assemble_matrix(mesh, A);
    dp2->assemble_matrix(mesh, B);
}
