#include "hermes2d.h"

using namespace WeakFormsH1;

/* Weak forms */

class WeakFormEigenLeft : public WeakForm
{
public:
  WeakFormEigenLeft() : WeakForm(1) {
    add_matrix_form(new DefaultJacobianDiffusion(0, 0));
    add_matrix_form(new MatrixFormPotential(0, 0));
  };

private:
  class MatrixFormPotential : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormPotential(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      for (int i=0;  i < n; i++) {
        Real x = e->x[i];
        Real y = e->y[i];
	result += wt[i] * (x*x + y*y) * u->val[i] * v->val[i];
      }
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };
};

class WeakFormEigenRight : public WeakForm
{
public:
  WeakFormEigenRight() : WeakForm(1) {
    add_matrix_form(new DefaultMatrixFormVol(0, 0));
  };
};
