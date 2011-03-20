#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform/sample_weak_forms.h"

// Exact solution.
class MyExactSolution : public ExactSolutionScalar
{
public:
  MyExactSolution(Mesh* mesh, double ALPHA_P, double X_LOC, double Y_LOC) 
         : ExactSolutionScalar(mesh), ALPHA_P(ALPHA_P), X_LOC(X_LOC), Y_LOC(Y_LOC) {};

  double fn(double x, double y) {
    return exp(-ALPHA_P * (pow((x - X_LOC), 2) + pow((y - Y_LOC), 2)));
  };

  // Function representing an exact scalar-valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    double a = -ALPHA_P * ( (x - X_LOC) * (x - X_LOC) + (y - Y_LOC) * (y - Y_LOC));

    dx = -exp(a) * (2 * ALPHA_P * (x - X_LOC));
    dy = -exp(a) * (2 * ALPHA_P * (y - Y_LOC));

    return fn(x, y);
  };

  // Members.
  double ALPHA_P;
  double X_LOC;
  double Y_LOC;
};

// Weak forms.
class MyWeakFormPoisson : public WeakFormLaplace
{
public:
  MyWeakFormPoisson(double ALPHA_P, double X_LOC, double Y_LOC) : WeakFormLaplace()
  {
    MyVectorFormVolPoisson* wfp = new MyVectorFormVolPoisson(0, ALPHA_P, X_LOC, Y_LOC);
    add_vector_form(wfp);
  };

private:
  class MyVectorFormVolPoisson : public WeakForm::VectorFormVol
  {
  public:
    MyVectorFormVolPoisson(int i, double ALPHA_P, double X_LOC, double Y_LOC) 
            : WeakForm::VectorFormVol(i), 
    ALPHA_P(ALPHA_P), X_LOC(X_LOC), Y_LOC(Y_LOC) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                       Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result -= wt[i] * (rhs<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                 Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, 
            ExtData<Ord> *ext) {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    template<typename Real>
    Real rhs(Real x, Real y) {
      Real a_P = (-ALPHA_P * pow((x - X_LOC), 2) - ALPHA_P * pow((y - Y_LOC), 2));
  
      return 4 * exp(a_P) * ALPHA_P * (ALPHA_P * (x - X_LOC) * (x - X_LOC) + ALPHA_P * (y - Y_LOC) * (y - Y_LOC) - 1);
    }
    
    // Members.
    double ALPHA_P;
    double X_LOC;
    double Y_LOC;
  };
};

// Boundary conditions (use values provided by exact solution).
class EssentialBCNonConstantExact : public EssentialBC
{
public:
  EssentialBCNonConstantExact(std::string marker, MyExactSolution* exact_solution) : 
        EssentialBC(Hermes::vector<std::string>()), exact_solution(exact_solution) 
  {
    markers.push_back(marker);
  };
  
  ~EssentialBCNonConstantExact() {};

  virtual EssentialBCValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar function(double x, double y) const {
    return exact_solution->fn(x, y);
  };

  MyExactSolution* exact_solution;
};