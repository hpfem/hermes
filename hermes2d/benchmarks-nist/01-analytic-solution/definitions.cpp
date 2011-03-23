#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform/sample_weak_forms.h"

// Right-hand side for the 2D equation -Laplace u = f with Dirichlet BC.
class MyRightHandSide
{
public:
  MyRightHandSide(double poly_deg) : poly_deg(poly_deg) {};

  template<typename Real>
  Real value(Real x, Real y) {
    Real a = pow(2.0, 4.0*poly_deg);
    Real b = pow(x-1.0, 8.0);
    Real c = (38.0*pow(x, 2.0) - 38.0*x + 9.0);
    Real d = pow(y-1.0, poly_deg);
    Real e = pow(y-1.0, 8.0);
    Real f = (38.0*pow(y, 2.0) - 38.0*y + 9.0);
    Real g = pow(x-1.0, poly_deg);

    return poly_deg*a*pow(x, 8.0)*b*c*pow(y, poly_deg)*d 
           + poly_deg*a*pow(y, 8.0)*e*f*pow(x, poly_deg)*g;
  }

  // Member.
  double poly_deg;
};

// Exact solution (needed for the Dirichlet condition).
class MyExactSolution : public ExactSolutionScalar
{
public:
  MyExactSolution(Mesh* mesh, double poly_deg) 
            : ExactSolutionScalar(mesh), poly_deg(poly_deg) {};

  // Exact solution.
  double value(double x, double y) {
    return pow(2, 4 * poly_deg) * pow(x, poly_deg) * pow(1 - x, poly_deg) 
           * pow(y, poly_deg) * pow(1 - y, poly_deg);
  }

  // Exact solution with derivatives.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    double A = pow((1.0-y), poly_deg);
    double B = pow((1.0-x), poly_deg);
    double C = pow(y, poly_deg);
    double D = pow(x, poly_deg);

    dx = ((poly_deg * pow(16.0, poly_deg)*A*C) / (x-1.0) 
         + (poly_deg*pow(16.0, poly_deg)*A*C)/x)*B*D;
    dy = ((poly_deg*pow(16.0, poly_deg)*B*D)/(y-1.0)
         + (poly_deg*pow(16.0, poly_deg)*B*D)/y)*A*C;

    return value(x, y);
  };

  // Member.
  double poly_deg;
};

class EssentialBCNonConst : public EssentialBC {
public:
  EssentialBCNonConst(std::string marker, MyExactSolution* exact_solution) 
              : EssentialBC(Hermes::vector<std::string>()),
    exact_solution(exact_solution) {
    markers.push_back(marker);
  }

  ~EssentialBCNonConst() { }

  inline EssentialBCValueType get_value_type() const { return EssentialBC::BC_FUNCTION; }

  scalar function(double x, double y) const {
    return exact_solution->value(x, y);
  }

  // Member.
  MyExactSolution* exact_solution;
};

class CustomWeakFormPoisson : public WeakFormLaplace
{
public:
  CustomWeakFormPoisson(MyRightHandSide* rhs) : WeakFormLaplace() {
    add_vector_form(new MyVectorFormVolPoisson(0, rhs));
  };

private:
  class MyVectorFormVolPoisson : public WeakForm::VectorFormVol
  {
  public:
    MyVectorFormVolPoisson(int i, MyRightHandSide* rhs) : WeakForm::VectorFormVol(i), 
          rhs(rhs) { }

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

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    // Member.
    MyRightHandSide* rhs;
  };
};

