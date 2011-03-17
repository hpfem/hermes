#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

// Exact solution.
#include "exact_solution.cpp"

class WeakFormPoisson : public WeakForm
{
public:
  WeakFormPoisson() : WeakForm(1) {
    add_matrix_form(new MatrixFormVolPoisson(0, 0));
  };

private:
  class MatrixFormVolPoisson : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolPoisson(int i, int j) : WeakForm::MatrixFormVol(i, j) {
      sym = HERMES_SYM;
    }

template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }
  };
};

class DirichletFunctionBoundaryConditionExact : public DirichletBoundaryCondition
{
public:
  DirichletFunctionBoundaryConditionExact(std::string marker, ExactSolutionNIST02* exact_solution) : 
        DirichletBoundaryCondition(Hermes::vector<std::string>()), exact_solution(exact_solution) 
  {
    markers.push_back(marker);
  };
  
  ~DirichletFunctionBoundaryConditionExact() {};

  virtual BoundaryConditionValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar function(double x, double y) const {
    return exact_solution->fn(x, y);
  };

  ExactSolutionNIST02* exact_solution;
};