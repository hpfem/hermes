#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform/sample_weak_forms.h"

// Exact solution.
class MyExactFunction
{
public:
  MyExactFunction(double k, double alpha) : k(k), alpha(alpha) {};

  double fn(double x, double y) {
    if (x <= 0) 
      return cos(k * y);
    else 
      return cos(k * y) + pow(x, alpha);
  }

  // Members.
  double k, alpha;
};

// Right-hand side.
class MyRightHandSide : public MyExactFunction
{
public:
  MyRightHandSide(double k, double alpha) : MyExactFunction(k, alpha) {};

  double value(double x, double y) {
    if (x < 0)
      return fn(x, y) * k * k;
    else 
      return fn(x, y) * k * k - alpha *(alpha - 1) * pow(x, alpha - 2.) - k * k * pow(x, alpha);
  }
};

// Exact solution.
class MyExactSolution : public ExactSolutionScalar, public MyExactFunction
{
public:
  MyExactSolution(Mesh* mesh, double k, double alpha) : ExactSolutionScalar(mesh), MyExactFunction(k, alpha) {};

  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    if (x <= 0) 
      dx = 0;
    else 
      dx = alpha * pow(x, alpha - 1);
    dy = -sin(k * y) * k;

    return fn(x, y);
  };
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
        result += wt[i] * rhs->value(e->x[i], e->y[i]) * v->val[i];
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(30);
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
    return exact_solution->fn(x, y);
  };

  MyExactSolution* exact_solution;
};
