#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

// Right-hand side for the 2D equation -Laplace u = f with Dirichlet BC.
class CustomRightHandSide
{
public:
  CustomRightHandSide(double K) : K(K) {};

  template<typename Real>
  Real value(Real x, Real y) {
    Real u = atan(this->K * x);
    Real dudx = 1. / (1 + (this->K * x) * (this->K * x)) * this->K;
    Real ddudxx = - this->K / (1 + (this->K * x) * (this->K * x)) / (1 + (this->K * x) * (this->K * x)) * 2. * this->K * this->K * x;
    return - ddudxx + u;
  }

  // Member.
  double K;
};

// Exact solution (needed in the Dirichlet condition).
class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double K) : ExactSolutionScalar(mesh), K(K) {};

  virtual scalar value (double x, double y) const {
    return atan(this->K * x);
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 1./(1 + (this->K * x)*(this->K * x)) * this->K;
    dy = 0;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(20);
  }

  // Members.
  double K;
};

// Weak forms.
class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(CustomRightHandSide* rhs, std::string bdy_left_right, double K) : WeakForm(1) {
    add_matrix_form(new CustomMatrixFormVol(0, 0));
    add_vector_form(new CustomVectorFormVol(0, rhs));
    add_vector_form_surf(new CustomVectorFormSurfRight(0, K, bdy_left_right));
    add_vector_form_surf(new CustomVectorFormSurfLeft(0, K, bdy_left_right));
  };

private:
  class CustomMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVol(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v) + 
             int_u_v<Real, Scalar>(n, wt, u, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class CustomVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVol(int i, CustomRightHandSide* rhs) : WeakForm::VectorFormVol(i), rhs(rhs) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * rhs->value(e->x[i], e->y[i]) * v->val[i];
      return result;
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Members.
    CustomRightHandSide* rhs;
  };

  class CustomVectorFormSurfRight : public WeakForm::VectorFormSurf
  {
  public:
    CustomVectorFormSurfRight(int i, double K, std::string area) : WeakForm::VectorFormSurf(i, area), K(K) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar dfdx_at_1 = 1. / (1 + (this->K * 1.) * (this->K * 1.)) * this->K;
      return - dfdx_at_1 * int_v<Real, Scalar>(n, wt, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Members.
    double K;
  };

  class CustomVectorFormSurfLeft : public WeakForm::VectorFormSurf
  {
  public:
    CustomVectorFormSurfLeft(int i, double K, std::string area) : WeakForm::VectorFormSurf(i, area), K(K) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
      Scalar dfdx_at_minus_1 = -1. / (1 + (-this->K * 1.) * (-this->K * 1.)) * this->K;
      return - dfdx_at_minus_1 * int_v<Real, Scalar>(n, wt, v);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Members.
    double K;
  };
};