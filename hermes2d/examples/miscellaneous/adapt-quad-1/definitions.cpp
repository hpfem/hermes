#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

class WeakFormPoisson : public WeakForm
{
public:
  // Problem parameters.
  double const_f;

  WeakFormPoisson(double const_f, bool adapt_eval, int adapt_order_increase, double adapt_rel_error_tol)
    : WeakForm(1), const_f(const_f) {
    MatrixFormVolPoisson* matrix_form = new MatrixFormVolPoisson(0, 0);
    VectorFormVolPoisson* vector_form = new VectorFormVolPoisson(0);
    if(adapt_eval) {
      matrix_form->adapt_order_increase = adapt_order_increase;
      vector_form->adapt_order_increase = adapt_order_increase;
      matrix_form->adapt_rel_error_tol = adapt_rel_error_tol;
      vector_form->adapt_rel_error_tol = adapt_rel_error_tol;
    }
    add_matrix_form(matrix_form);
    add_vector_form(vector_form);
  };

private:
  class MatrixFormVolPoisson : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolPoisson(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class VectorFormVolPoisson : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolPoisson(int i) : WeakForm::VectorFormVol(i) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      return static_cast<WeakFormPoisson *>(wf)->const_f * int_v<Real, Scalar>(n, wt, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
  };
};
