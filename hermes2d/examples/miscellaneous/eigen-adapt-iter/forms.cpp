#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormS : public WeakForm
{
public:
  WeakFormS() : WeakForm(1) {
    add_matrix_form(new MatrixFormVolS(0, 0));
  };

private:
  class MatrixFormVolS : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolS(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      scalar result = 0;
      for (int i = 0; i < n; i++) {
        double x = e->x[i];
        double y = e->y[i];
        result += wt[i] * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i] 
                           + V(x, y) * u->val[i] * v->val[i]);
      }
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(30);
    }

    double V(double x, double y) {
      return 0;
      //double r = sqrt(x*x + y*y);
      //return -1./(0.001 + r*r);
    }
  };
};

class WeakFormM : public WeakForm
{
public:
  WeakFormM() : WeakForm(1) {
    add_matrix_form(new MatrixFormVolM(0, 0));
  };

private:
  class MatrixFormVolM : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolM(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };
};