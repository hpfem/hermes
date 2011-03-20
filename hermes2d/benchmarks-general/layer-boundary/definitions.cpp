#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

// Exact solution to the 1D problem -u'' + K*K*u = K*K in (-1,1) with zero Dirichlet BC.
class MyExactFunction
{
public:
  MyExactFunction(double K) : K(K) {};

  double uhat(double x) {
    return 1. - (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
  }
  double duhat_dx(double x) {
    return -K * (exp(K*x) - exp(-K*x)) / (exp(K) + exp(-K));
  }
  double dduhat_dxx(double x) {
    return -K*K * (exp(K*x) + exp(-K*x)) / (exp(K) + exp(-K));
  }

  // Member.
  double K;
};

// Right-hand side for the 2D equation -Laplace u + K*K*u = 0 with zero Dirichlet BC.
class MyRightHandSide : public MyExactFunction
{
public:
  MyRightHandSide(double K) : MyExactFunction(K) {};

  double value(double x, double y) {
    return -(dduhat_dxx(x)*uhat(y) + uhat(x)*dduhat_dxx(y)) + K*K*uhat(x)*uhat(y);
  }
};

// Exact solution to the 2D equation.
class MyExactSolution : public ExactSolutionScalar, public MyExactFunction
{
public:
  MyExactSolution(Mesh* mesh, double K) : ExactSolutionScalar(mesh), MyExactFunction(K) {};

  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    dx = duhat_dx(x) * uhat(y);
    dy = uhat(x) * duhat_dx(y);
    return uhat(x) * uhat(y);
  };
};

// Weak forms for the 2D equation with Dirichlet boundary conditions.
class MyWeakFormPerturbedPoisson : public WeakForm
{
public:
  MyWeakFormPerturbedPoisson(MyRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new MyMatrixFormVolPerturbedPoisson(0, 0, rhs->K));
    add_vector_form(new MyVectorFormVolPerturbedPoisson(0, rhs));
  };

private:
  class MyMatrixFormVolPerturbedPoisson : public WeakForm::MatrixFormVol
  {
  public:
    MyMatrixFormVolPerturbedPoisson(int i, int j, double K) : WeakForm::MatrixFormVol(i, j), K(K) {
      sym = HERMES_SYM;
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) + K*K * int_u_v<Real, 
             Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Member.
    double K;
  };

  class MyVectorFormVolPerturbedPoisson : public WeakForm::VectorFormVol
  {
  public:
    MyVectorFormVolPerturbedPoisson(int i, MyRightHandSide* rhs) : WeakForm::VectorFormVol(i), 
          rhs(rhs) { }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (rhs->value(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) {
      return Ord(24);
    }

    // Member.
    MyRightHandSide* rhs;
  };
};
