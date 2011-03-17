#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

// Exact solution.
#include "exact_functions.cpp"

class WeakFormPerturbedPoisson : public WeakForm
{
public:
  WeakFormPerturbedPoisson(MyRightHandSide* my_rhs) : WeakForm(1) {
    add_matrix_form(new MatrixFormVolPerturbedPoisson(0, 0, my_rhs->K));
    add_vector_form(new VectorFormVolPerturbedPoisson(0, my_rhs));
  };

private:
  class MatrixFormVolPerturbedPoisson : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolPerturbedPoisson(int i, int j, double K) : WeakForm::MatrixFormVol(i, j), K(K) {
      sym = HERMES_SYM;
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) + K*K * int_u_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Member.
    double K;
  };

  class VectorFormVolPerturbedPoisson : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolPerturbedPoisson(int i, MyRightHandSide* my_rhs) : WeakForm::VectorFormVol(i), 
          my_rhs(my_rhs) { }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (my_rhs->rhs(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(24);
    }

    // Member.
    MyRightHandSide* my_rhs;
  };
};
