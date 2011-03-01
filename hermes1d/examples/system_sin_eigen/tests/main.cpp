#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_VERBOSE
#define HERMES_REPORT_FILE "application.log"
#include "hermes1d.h"
#include "h1d_wrapper_api.h"

// This test makes sure that example "system_sin_eigen" works correctly.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int NELEM = 2;                            // Number of elements.
static int NEQ = 2;
double A = 0, B = M_PI;                   // Domain end points.
int P_init = 7;                           // Initial polynomal degree.

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
  // Create mesh, set Dirichlet BC, enumerate basis functions.
  Mesh *mesh = new Mesh(A, B, NELEM, P_init, NEQ);
  mesh->set_bc_left_dirichlet(0, 0);
  mesh->set_bc_right_dirichlet(0, 0);

  // Obtain the number of degrees of freedom.
  int N_dof = mesh->assign_dofs();
  info("ndofs: %d", N_dof);

  // Initialize the FE problem.
  DiscreteProblem *dp1 = new DiscreteProblem();
  dp1->add_matrix_form(0, 1, A_01);
  dp1->add_matrix_form(1, 0, A_10);
  DiscreteProblem *dp2 = new DiscreteProblem();
  dp2->add_matrix_form(0, 0, B_00);
  dp2->add_matrix_form(1, 1, B_00);

  // Allocate matrices and vector for previous solution.
  CooMatrix *mat1 = new CooMatrix(N_dof);
  CooMatrix *mat2 = new CooMatrix(N_dof);

  // Assemble the Jacobian matrices and residual vectors.
  dp1->assemble_matrix(mesh, mat1);
  dp2->assemble_matrix(mesh, mat2);

  Python p;

  p.exec("print 'Python initialized'");
  p.push("A", c2py_CooMatrix(mat1));
  p.push("B", c2py_CooMatrix(mat2));
  p.exec("from hermes1d.solvers.eigen import solve_eig_numpy, solve_eig_pysparse");
  p.exec("eigs = solve_eig_numpy(A.to_scipy_coo(), B.to_scipy_coo())");
  p.exec("from utils import show_eigs");
  p.exec("show_eigs(eigs)");

  if (import_hermes1d__h1d_wrapper__h1d_wrapper())
      throw std::runtime_error("Can't import hermes1d");
  p.push("mesh",  c2py_Mesh(mesh));
  printf("2\n");
  p.exec("from plot import plot_eigs, plot_file");
  p.exec("plot_eigs(mesh, eigs)");

  return ERROR_FAILURE;
}
