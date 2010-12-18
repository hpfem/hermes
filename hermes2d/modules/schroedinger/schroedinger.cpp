#include "schroedinger.h"

namespace {

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


void ModuleSchroedinger::assemble(const Ptr<Matrix> &A, const Ptr<Matrix> &B)
{
    global_potential = this->potential;
}
