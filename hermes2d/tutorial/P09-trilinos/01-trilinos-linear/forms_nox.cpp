#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormPoissonNox : public WeakForm
{
public:
  WeakFormPoissonNox(bool is_matfree) : WeakForm(1) {
    this->is_matfree = is_matfree;
    add_matrix_form(new MatrixFormVolPoissonNox(0, 0));
    add_vector_form(new VectorFormVolPoissonNox(0));
  };

private:
  class MatrixFormVolPoissonNox : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolPoissonNox(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class VectorFormVolPoissonNox : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolPoissonNox(int i) : WeakForm::VectorFormVol(i) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Func<Scalar>* u = u_ext[0];
      Scalar result = 0;
      for (int i = 0; i < n; i++) {
        result += wt[i] * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
        result -= wt[i] * (rhs<Real>(e->x[i], e->y[i]) * v->val[i]);
      }
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    template<typename Real>
    Real rhs(Real x, Real y) {
      return -4.0;
    }
  };
};