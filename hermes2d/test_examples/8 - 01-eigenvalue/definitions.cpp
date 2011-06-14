#include "definitions.h"

WeakFormEigenLeft::WeakFormEigenLeft() : WeakForm(1) 
{
  add_matrix_form(new WeakFormsH1::DefaultJacobianDiffusion(0, 0));
  add_matrix_form(new MatrixFormPotential(0, 0));
}

template<typename Real, typename Scalar>
Scalar WeakFormEigenLeft::MatrixFormPotential::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                                                           Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
{
  Scalar result = 0;
  for (int i = 0; i < n; i++) 
  {
    Real x = e->x[i];
    Real y = e->y[i];
    result += wt[i] * (x*x + y*y) * u->val[i] * v->val[i];
  }
  return result;
}

scalar WeakFormEigenLeft::MatrixFormPotential::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                                                     Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
{
  return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
}

Ord WeakFormEigenLeft::MatrixFormPotential::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                                                Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
{
  return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
}


WeakFormEigenRight::WeakFormEigenRight() : WeakForm(1) 
{
  add_matrix_form(new WeakFormsH1::DefaultMatrixFormVol(0, 0));
}

