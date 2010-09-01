#include "hermes1d.h"

#include "python_api.h"

#include "h1d_wrapper_api.h"

int N_elem = 2;                         // number of elements
static int N_eq = 2;
double A = 0, B = M_PI;                     // domain end points
int P_init = 7;                           // initial polynomal degree

/*
This example solves u'' + k^2 u = 0 as an eigenvalue problem, by reducing it to
two first order equations using:

u = u1
u' = u1' = k*u2
u'' = k*u2'

so:

     - u2' = k*u1
 u1'       = k*u2

in an interval (0, pi) equipped with Dirichlet bdy conditions
u(0) = 0, u(pi) = 0
The exact solutions are u(x) = sin(k*x) for k = 0, +-1, +-2, +-3, ...

The weak formulation is:

\int          - u2' * v1 = k * \int u1 * v1
\int u1' * v2            = k * \int u2 * v2

We solve A*x = E*B*x, where A, B are matrices composed of 2x2 blocks.

*/


double kappa = 1;

double A_01(int num, double *x, double *weights,
                double *u, double *dudx, double *v, double *dvdx,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                void *user_data)
{
    double val = 0;
    for(int i = 0; i<num; i++) {
        double r = x[i];
        val += (-u[i]*dvdx[i])*weights[i];
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
        val += (u[i]*dvdx[i])*weights[i];
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

/******************************************************************************/
int main(int argc, char* argv[]) {
  // create mesh
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  // you can set the zero dirichlet at the right hand side
  mesh->set_bc_left_dirichlet(0, 0);
  mesh->set_bc_right_dirichlet(0, 0);

  // variable for the total number of DOF
  int N_dof = mesh->assign_dofs();
  printf("ndofs: %d\n", N_dof);

  // register weak forms
  DiscreteProblem *dp1 = new DiscreteProblem();
  dp1->add_matrix_form(0, 1, A_01);
  dp1->add_matrix_form(1, 0, A_10);
  DiscreteProblem *dp2 = new DiscreteProblem();
  dp2->add_matrix_form(0, 0, B_00);
  dp2->add_matrix_form(1, 1, B_00);

  CooMatrix *mat1 = new CooMatrix(N_dof);
  CooMatrix *mat2 = new CooMatrix(N_dof);
  double *y_prev = new double[N_dof];

  dp1->assemble_matrix(mesh, mat1);
  dp2->assemble_matrix(mesh, mat2);

  Python p;

  p.exec("print 'Python initialized'");
  p.push("A", c2py_CooMatrix(mat1));
  p.push("B", c2py_CooMatrix(mat2));
  p.exec("from hermes1d.solvers.eigen import solve_eig_numpy, solve_eig_pysparse");
  p.exec("eigs = solve_eig_numpy(A.to_scipy_coo(), B.to_scipy_coo())");
  //p.exec("eigs = solve_eig_pysparse(A.to_scipy_coo(), B.to_scipy_coo())");
  p.exec("from utils import show_eigs");
  p.exec("show_eigs(eigs)");

  if (import_hermes1d__h1d_wrapper__h1d_wrapper())
      throw std::runtime_error("Can't import hermes1d");
  p.push("mesh",  c2py_Mesh(mesh));
  printf("2\n");
  p.exec("from plot import plot_eigs, plot_file");
  p.exec("plot_eigs(mesh, eigs)");
  printf("Done.\n");
  return 0;
}
