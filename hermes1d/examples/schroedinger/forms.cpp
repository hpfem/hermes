#include "forms.h"
#undef _XOPEN_SOURCE
#include "_forms_api.h"

static double l = 0;
static double Z = 1;

/*
double Z_eff(double x)
{
    return 47.0-x/50;
}
*/

double lhs_R(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        double coeff;
        double r = x[i];
        // specify r^2*V:
        double rrV = r*r*potential_V(r); // r^2 * (-1/r) = -r
        // Hydrogen:
        //double rrV = -Z_eff(r)*r; // r^2 * (-1/r) = -r
        // Harmonic oscillator:
        //double rrV = r*r*r*r; // r^2 * (r^2) = r^4
        coeff = 0.5*r*r*dudx[i]*dvdx[i] + (rrV + 0.5 * (l + 1)*l) *u[i]*v[i];
        val += coeff*weights[i];
    }
    return val;
}

double rhs_R(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for(int i = 0; i < num; i++) {
        double r = x[i];
        val += u[i]*v[i]*r*r*weights[i];
    }
    return val;
}

/* ------------------------------------------------- */

double lhs_rR(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        double coeff;
        double r = x[i];
        // specify V:
        // Hydrogen:
        if (fabs(r) < 1e-12)
            throw std::runtime_error("Internal error: r is too small");
        double V = -Z/r;
        // Harmonic oscillator:
        //double V = r*r;
        coeff = 0.5*dudx[i]*dvdx[i] + (V + (l + 1)*l/(2*r*r)) *u[i]*v[i];
        val += coeff*weights[i];
    }
    return val;
}

double rhs_rR(int num, double *x, double *weights,
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

static DiscreteProblem *dp1 = NULL;
static DiscreteProblem *dp2 = NULL;
static int _equation_type = -1;

/*
   Assembles the discrete problem A*x = E*B*x
*/
void assemble_schroedinger(Mesh *mesh, Matrix *A, Matrix *B, int _l, int _Z,
        int equation_type)
{
    import__forms();
    l = _l;
    Z = _Z;
    if (_equation_type == -1) {
        // Initialize the problem
        dp1 = new DiscreteProblem();
        dp2 = new DiscreteProblem();
        if (equation_type == eqn_type_R) {
            dp1->add_matrix_form(0, 0, lhs_R);
            dp2->add_matrix_form(0, 0, rhs_R);
        } else if (equation_type == eqn_type_rR) {
            dp1->add_matrix_form(0, 0, lhs_rR);
            dp2->add_matrix_form(0, 0, rhs_rR);
        } else
            throw std::runtime_error("Unknown equation type");
        _equation_type = equation_type;
    } else {
        // Make sanity checks
        if (_equation_type != equation_type)
            throw std::runtime_error("Equation types mismatch");
    }

    int N_dof = mesh->assign_dofs();
    printf("Assembling A, B. ndofs: %d\n", N_dof);
    dp1->assemble_matrix(mesh, A);
    printf(" A is done.\n");
    dp2->assemble_matrix(mesh, B);
    printf(" B is done.\n");
    printf("  Done assembling.\n");
}
