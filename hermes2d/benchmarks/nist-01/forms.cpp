#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

// Exact solution.
#include "exact_solution.cpp"

class WeakFormPoisson : public WeakForm
{
public:
  WeakFormPoisson(double EXACT_SOL_P) : WeakForm(1)
  {
    add_matrix_form(new MatrixFormVolPoisson(0, 0));
    
    VectorFormVolPoisson* wfp= new VectorFormVolPoisson(0);
    wfp->EXACT_SOL_P = EXACT_SOL_P;
    add_vector_form(wfp);
  };

private:
  class MatrixFormVolPoisson : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolPoisson(int i, int j) : WeakForm::MatrixFormVol(i, j) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };

  class VectorFormVolPoisson : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolPoisson(int i) : WeakForm::VectorFormVol(i) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (rhs<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
    }

    template<typename Real>
    Real rhs(Real x, Real y)
    {
      Real a = pow(2.0, 4.0*EXACT_SOL_P);
      Real b = pow(x-1.0, 8.0);
      Real c = (38.0*pow(x, 2.0) - 38.0*x + 9.0);
      Real d = pow(y-1.0, EXACT_SOL_P);
      Real e = pow(y-1.0, 8.0);
      Real f = (38.0*pow(y, 2.0) - 38.0*y + 9.0);
      Real g = pow(x-1.0, EXACT_SOL_P);

      return EXACT_SOL_P*a*pow(x, 8.0)*b*c*pow(y, EXACT_SOL_P)*d + EXACT_SOL_P*a*pow(y, 8.0)*e*f*pow(x,EXACT_SOL_P)*g;
    };
    
    // Members.
    double EXACT_SOL_P;
  };
};

class DirichletFunctionBoundaryCondition : public DirichletBoundaryCondition
{
public:
  DirichletFunctionBoundaryCondition(std::string marker) : 
        DirichletBoundaryCondition(Hermes::vector<std::string>()){
    markers.push_back(marker);
  };
  
  ~DirichletFunctionBoundaryCondition() {};

  virtual BoundaryConditionValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar function(double x, double y) const {
    return fn(x, y);
  };
};