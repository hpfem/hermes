#include "schroedinger.h"

namespace {

using Teuchos::RCP;
using Teuchos::Ptr;
using Teuchos::rcp;
using Teuchos::null;

using Schroedinger::Potential;

RCP<Potential> global_potential=null;

double bilinear_form_left(int n, double *wt, Func<double> *u_ext[],
 Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<double> *ext)
{
  double result = 0;
  for (int i = 0; i < n; i++) {
    double x = e->x[i];
    double y = e->y[i];
    double V = global_potential->get_value(x, y);
    result += wt[i] * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]
                       + V * u->val[i] * v->val[i]);
  }
  return result;
}

Ord bilinear_form_left_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                           Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(20);
/*
  return wt[i] * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i]
             + u->val[i] * v->val[i])
*/
}

template<typename Real, typename Scalar>
Scalar bilinear_form_right(int n, double *wt, Func<Scalar> *u_ext[],
     Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_u_v<Real, Scalar>(n, wt, u, v);
}

} // anonymous namespace


namespace Schroedinger {


void ModuleSchroedinger::assemble(const Ptr<SparseMatrix> &A,
        const Ptr<SparseMatrix> &B)
{

    const int BDY_BOTTOM = 1;
    const int BDY_RIGHT = 2;
    const int BDY_TOP = 3;
    const int BDY_LEFT = 4;
    int P_INIT = 2;

    Mesh mesh;
    H2DReader mloader;
    mloader.load("domain.mesh", &mesh);
    mesh.refine_all_elements();
    mesh.refine_all_elements();
    mesh.refine_all_elements();

    BCTypes bc_types;
    bc_types.add_bc_dirichlet(Hermes::vector<int>(BDY_BOTTOM, BDY_RIGHT,
                BDY_TOP, BDY_LEFT));
    BCValues bc_values;
    bc_values.add_zero(Hermes::vector<int>(BDY_BOTTOM, BDY_RIGHT,
                BDY_TOP, BDY_LEFT));

    H1Space space(&mesh, &bc_types, &bc_values, P_INIT);
    int ndof = Space::get_num_dofs(&space);
    printf("ndof: %d\n", ndof);

    global_potential = this->potential;
    WeakForm wf_left, wf_right;
    wf_left.add_matrix_form(bilinear_form_left, bilinear_form_left_ord);
    wf_right.add_matrix_form(callback(bilinear_form_right));
    DiscreteProblem dp_left(&wf_left, &space, true);
    DiscreteProblem dp_right(&wf_right, &space, true);
    dp_left.assemble(A.getRawPtr());
    dp_right.assemble(B.getRawPtr());
}

} // namespace Schroedinger
