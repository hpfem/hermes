#include "weakform/weakform.h"
#include "weakform_library/laplace.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

/* Exact solution */

// Exact solution to the 1D problem -u'' + K*K*u = K*K in (-1,1) with zero Dirichlet BC.
class CustomExactFunction
{
public:
  CustomExactFunction(double K) : K(K) {};

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

// Exact solution to the 2D equation.
class CustomExactSolution : public ExactSolutionScalar, public CustomExactFunction
{
public:
  CustomExactSolution(Mesh* mesh, double K) : ExactSolutionScalar(mesh), CustomExactFunction(K) {};

  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    dx = duhat_dx(x) * uhat(y);
    dy = uhat(x) * duhat_dx(y);
    return uhat(x) * uhat(y);
  };
};

/* Right-hand side */

// Right-hand side for the 2D equation -Laplace u + K*K*u = 0 with zero Dirichlet BC.
class CustomRightHandSide : public CustomExactFunction
{
public:
  CustomRightHandSide(double K) : CustomExactFunction(K) {};

  double value(double x, double y) {
    return -(dduhat_dxx(x)*uhat(y) + uhat(x)*dduhat_dxx(y)) + K*K*uhat(x)*uhat(y);
  }
};

/* Weak forms */

// Weak forms for the 2D equation with Dirichlet boundary conditions.
class CustomWeakFormPerturbedPoisson : public WeakForm
{
public:
  CustomWeakFormPerturbedPoisson(CustomRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new DefaultMatrixFormVolConst(0, 0));
    add_matrix_form(new DefaultMatrixFormVolMassConst(0, 0, rhs->K*rhs->K));
    add_vector_form(new CustomVectorFormVolPerturbedPoisson(0, rhs));
  };

private:
  class CustomVectorFormVolPerturbedPoisson : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVolPerturbedPoisson(int i, CustomRightHandSide* rhs) 
          : WeakForm::VectorFormVol(i), rhs(rhs) { }

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
    CustomRightHandSide* rhs;
  };
};
