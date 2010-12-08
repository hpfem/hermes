#include "hermes1d.h"

#include "python_api.h"
#include "h1d_wrapper_api.h"
#include "forms.h"

int N_elem = 40;                         // number of elements
static int N_eq = 2;
double A = 0, B = 50;                     // domain end points
int P_init = 5;                           // initial polynomal degree
int kappa = 1;
int Z = 1;


/******************************************************************************/
int main(int argc, char* argv[]) {
  // create mesh
  Mesh *mesh = new Mesh(A, B, N_elem, P_init, N_eq);
  // We need to do this here in order to determine N_dof, otherwise it is not
  // needed, as assemble_dirac() will do it again
  mesh->set_bc_left_dirichlet(0, 0);
  mesh->set_bc_right_dirichlet(0, 0);
  int N_dof = mesh->assign_dofs();

  CooMatrix *mat1 = new CooMatrix(N_dof);
  CooMatrix *mat2 = new CooMatrix(N_dof);

  assemble_dirac(mesh, mat1, mat2, kappa, Z);

  Python p;

  p.exec("print 'Python initialized'");
  p.push("A", c2py_CooMatrix(mat1));
  p.push("B", c2py_CooMatrix(mat2));
  p.exec("from utils import solve_eig_numpy, solve_eig_pysparse");
  p.exec("eigs = solve_eig_numpy(A.to_scipy_coo(), B.to_scipy_coo())");
  //p.exec("eigs = solve_eig_pysparse(A.to_scipy_coo(), B.to_scipy_coo())");
  p.exec("from utils import show_eigs");
  p.exec("show_eigs(eigs)");

  if (import_hermes1d__h1d_wrapper__h1d_wrapper())
      throw std::runtime_error("Can't import hermes1d");
  p.push("mesh",  c2py_Mesh(mesh));
  p.exec("from plot import plot_eigs, plot_file");
  p.exec("plot_eigs(mesh, eigs)");
  printf("Done.\n");
  return 0;
}
