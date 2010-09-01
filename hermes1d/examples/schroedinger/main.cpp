#include "hermes1d.h"

#include "python_api.h"

#include "h1d_wrapper_api.h"

static int N_eq = 1;
int N_elem = 40;                         // number of elements
double A = 0, B = 20;                     // domain end points
int P_init = 2;                           // initial polynomal degree

double l = 0;

double lhs(int num, double *x, double *weights,
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
        // Hydrogen:
        double rrV = -r; // r^2 * (-1/r) = -r
        // Harmonic oscillator:
        //double rrV = r*r*r*r; // r^2 * (r^2) = r^4
        coeff = 0.5*r*r*dudx[i]*dvdx[i] + (rrV + 0.5 * (l + 1)*l) *u[i]*v[i];
        val += coeff*weights[i];
    }
    return val;
}

double rhs(int num, double *x, double *weights,
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

double E;

double residual(int num, double *x, double *weights,
                double u_prev[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double du_prevdx[MAX_SLN_NUM][MAX_EQN_NUM][MAX_QUAD_PTS_NUM],
                double *v, double *dvdx, void *user_data)
{
    double *u = &u_prev[0][0][0];
    double *dudx = &du_prevdx[0][0][0];
    double val = 0;
    for(int i = 0; i<num; i++) {
        double coeff;
        coeff = 0.5*x[i]*x[i]*dudx[i]*dvdx[i] -u[i]*v[i]*x[i]
            + 0.5 * (l + 1)*l *u[i]*v[i];
        coeff -= E*u[i]*v[i]*x[i]*x[i];
        val += coeff*weights[i];
    }
    return val;
}

/******************************************************************************/
int main(int argc, char* argv[]) {
  // create mesh
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  // you can set the zero dirichlet at the right hand side
  //mesh.set_bc_right_dirichlet(0, 0);

  // variable for the total number of DOF
  int N_dof = mesh->assign_dofs();
  printf("ndofs: %d\n", N_dof);

  // register weak forms
  DiscreteProblem *dp1 = new DiscreteProblem();
  dp1->add_matrix_form(0, 0, lhs);
  DiscreteProblem *dp2 = new DiscreteProblem();
  dp2->add_matrix_form(0, 0, rhs);

  DiscreteProblem *dp3 = new DiscreteProblem();
  dp3->add_vector_form(0, residual);

  // allocate Jacobi matrix and residual
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
  p.exec("E, v = eigs[0]");
  int n;

  double *res = new double[N_dof];
  E = py2c_double(p.pull("E"));
  printf("E=%.10f\n", E);
  E = -0.5;
  dp3->assemble_vector(mesh, res);
  // calculate L2 norm of residual vector
  double res_norm = 0;
  for(int i=0; i<N_dof; i++) res_norm += res[i]*res[i];
  res_norm = sqrt(res_norm);
  printf("L2 norm of the residual: %f\n", res_norm);


  Linearizer l(mesh);
  const char *out_filename = "solution.gp";
  l.plot_solution(out_filename);

  printf("still ok\n");
  if (import_hermes1d__h1d_wrapper__h1d_wrapper())
      throw std::runtime_error("Can't import hermes1d");
  p.push("mesh",  c2py_Mesh(mesh));
  printf("2\n");
  p.exec("from plot import plot_eigs, plot_file");
  p.exec("plot_eigs(mesh, eigs)");
  printf("Done.\n");
  return 0;
}
