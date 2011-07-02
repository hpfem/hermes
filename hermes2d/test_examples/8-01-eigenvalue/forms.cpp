#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormEigenLeft : public WeakForm
{
public:
  WeakFormEigenLeft() : WeakForm(1) {
    add_matrix_form(new MatrixFormVolLeft(0, 0));
  };

private:
  class MatrixFormVolLeft : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolLeft(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) {
      return matrix_form<Scalar, Scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };
};

class WeakFormEigenRight : public WeakForm
{
public:
  WeakFormEigenRight() : WeakForm(1) {
    add_matrix_form(new MatrixFormVolRight(0, 0));
  };

private:
  class MatrixFormVolRight : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolRight(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_u_v<Real, Scalar>(n, wt, u, v);
    }

    Scalar value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) {
      return matrix_form<Scalar, Scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };
};