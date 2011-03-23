#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform/sample_weak_forms.h"

// Right-hand side.
class MyRightHandSide
{
public:
  MyRightHandSide(double alpha, double x_loc, double y_loc, double r_zero) 
    : alpha(alpha), x_loc(x_loc), y_loc(y_loc), r_zero(r_zero) { };

  template<typename Real>
  Real value(Real x, Real y) {  
    Real a = pow(x - x_loc, 2);
    Real b = pow(y - y_loc, 2);
    Real c = sqrt(a + b);
    Real d = ((alpha*x - (alpha * x_loc)) * (2*x - (2 * x_loc)));
    Real e = ((alpha*y - (alpha * y_loc)) * (2*y - (2 * y_loc)));
    Real f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);
    Real g = (alpha * c - (alpha * r_zero));

    return ((alpha/(c * f)) - (d/(2 * pow(a + b, 1.5) * f)) - ((alpha * d * g)/((a + b) * pow(f, 2))) + 
           (alpha/(c * f)) - (e/(2 * pow(a + b, 1.5) * f)) - ((alpha * e * g)/((a + b) * pow(f, 2))));
  }

  double alpha, x_loc, y_loc, r_zero;
};


// Exact solution.
class MyExactSolution : public ExactSolutionScalar
{
public:
  MyExactSolution(Mesh* mesh, double alpha, double x_loc, double y_loc, double r_zero) 
             : ExactSolutionScalar(mesh), alpha(alpha), x_loc(x_loc), y_loc(y_loc), r_zero(r_zero) { }

  // Exact solution value.
  double value(double x, double y) {
    return atan(alpha * (sqrt(pow(x - x_loc, 2) + pow(y - y_loc, 2)) - r_zero));
  };

  // Exact solution value and derivative.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    double a = pow(x - x_loc, 2);
    double b = pow(y - y_loc, 2);
    double c = sqrt(a + b);
    double d = (alpha*x - (alpha * x_loc));
    double e = (alpha*y - (alpha * y_loc));
    double f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);

    dx = (d/(c * f));
    dy = (e/(c * f));

    return value(x, y);
  };

  // Members.
  double alpha, x_loc, y_loc, r_zero;
};

class CustomWeakFormPoisson : public WeakFormLaplace
{
public:
  CustomWeakFormPoisson(MyRightHandSide* rhs) : WeakFormLaplace()
  {
    add_vector_form(new MyVectorFormVolPoisson(0, rhs));
  };

private:
  class MyVectorFormVolPoisson : public WeakForm::VectorFormVol
  {
  public:
    MyVectorFormVolPoisson(int i, MyRightHandSide* rhs) : WeakForm::VectorFormVol(i), rhs(rhs) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result -= wt[i] * (rhs->value<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }
    
    // Members.
    MyRightHandSide* rhs;
  };
};

// Essential boundary conditions.
class EssentialBCNonConstant : public EssentialBC
{
public:
  EssentialBCNonConstant(std::string marker, MyExactSolution* exact_solution) : 
        EssentialBC(Hermes::vector<std::string>()), exact_solution(exact_solution) 
  {
    markers.push_back(marker);
  };
  
  ~EssentialBCNonConstant() {};

  virtual EssentialBCValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar function(double x, double y) const {
    return exact_solution->value(x, y);
  };

  MyExactSolution* exact_solution;
};
