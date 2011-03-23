#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

// Right-hand side for the 2D equation -Laplace u = f with Dirichlet BC.
class CustomRightHandSide
{
public:
  CustomRightHandSide(double slope) : slope(slope) {};

  template<typename Real>
  Real value(Real x, Real y) {
    Real t2 = sqr(y + 0.25) + sqr(x - 1.25);
    Real t = sqrt(t2);
    Real u = (sqr(M_PI - 3.0*t)*sqr(slope) + 9.0);
    return 27.0/2.0 * sqr(2.0*y + 0.5) * (M_PI - 3.0*t) * pow(slope, 3.0) / (sqr(u) * t2) +
           27.0/2.0 * sqr(2.0*x - 2.5) * (M_PI - 3.0*t) * pow(slope, 3.0) / (sqr(u) * t2) -
           9.0/4.0 * sqr(2.0*y + 0.5) * slope / (u * pow(t,3.0)) -
           9.0/4.0 * sqr(2.0*x - 2.5) * slope / (u * pow(t,3.0)) +
           18.0 * slope / (u * t);
  }

  // Member.
  double slope;
};

// Exact solution (needed in the Dirichlet condition).
class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double slope) : ExactSolutionScalar(mesh), slope(slope) {};

  // Exact solution.
  double value(double x, double y) {
    return atan(slope * (sqrt(sqr(x-1.25) + sqr(y+0.25)) - M_PI/3));
  }

  // Exact solution with derivatives.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
      double t = sqrt(sqr(x-1.25) + sqr(y+0.25));
      double u = t * (sqr(slope) * sqr(t - M_PI/3) + 1);
      dx = slope * (x-1.25) / u;
      dy = slope * (y+0.25) / u;
      return value(x, y);
  };

  // Member.
  double slope;
};

class EssentialBCNonConst : public EssentialBC {
public:
  EssentialBCNonConst(std::string marker, CustomExactSolution* exact_solution) : EssentialBC(Hermes::vector<std::string>()),
    exact_solution(exact_solution) {
    markers.push_back(marker);
  }

  ~EssentialBCNonConst() { }

  inline EssentialBCValueType get_value_type() const { return EssentialBC::BC_FUNCTION; }

  scalar function(double x, double y) const {
    return exact_solution->value(x, y);
  }

  // Member.
  CustomExactSolution* exact_solution;
};

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(CustomRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new CustomMatrixFormVolPoisson(0, 0));
    add_vector_form(new CustomVectorFormVolPoisson(0, rhs));
  };

private:
  class CustomMatrixFormVolPoisson : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVolPoisson(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_SYM) {}

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class CustomVectorFormVolPoisson : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVolPoisson(int i, CustomRightHandSide* rhs) : WeakForm::VectorFormVol(i), rhs(rhs) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result -= wt[i] * (rhs->value<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Member.
    CustomRightHandSide* rhs;
  };
};

