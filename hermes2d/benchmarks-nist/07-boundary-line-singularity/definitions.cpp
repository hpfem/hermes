#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform/sample_weak_forms.h"

// Exact solution.
class MyExactSolution : public ExactSolutionScalar
{
public:
  MyExactSolution(Mesh* mesh, double alpha) : ExactSolutionScalar(mesh), alpha(alpha) {};

  // Exact solution.
  double value(double x, double y) {
    return pow(x, alpha);
  };

  // Exact solution with derivatives.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    dx = (alpha/(pow(x, 0.4)));
    dy = 0;
    return value(x, y);
  };

  // Members.
  double alpha;
};

class CustomWeakFormPoisson : public WeakFormLaplace
{
public:
  CustomWeakFormPoisson() : WeakFormLaplace() {
    add_vector_form(new MyVectorFormVol(0));
  };

private:
  class MyVectorFormVol : public WeakForm::VectorFormVol
  {
  public:
    MyVectorFormVol(int i) : WeakForm::VectorFormVol(i) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result -= wt[i] * (rhs<Real>(e->x[i], e->y[i]) * v->val[i]);
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
      return (-0.24/(pow(x, 1.4)));
    }
  };
};

// Essential boundary conditions.
class EssentialBCNonConst : public EssentialBC
{
public:
  EssentialBCNonConst(std::string marker, MyExactSolution* exact_solution) : 
        EssentialBC(Hermes::vector<std::string>()), exact_solution(exact_solution) 
  {
    markers.push_back(marker);
  };
  
  ~EssentialBCNonConst() {};

  virtual EssentialBCValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar function(double x, double y) const {
    return exact_solution->value(x, y);
  };

  MyExactSolution* exact_solution;
};
