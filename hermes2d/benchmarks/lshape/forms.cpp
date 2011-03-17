#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

// Exact solution.
#include "exact_solution.cpp"

class WeakFormLinearLShape : public WeakForm
{
public:
  // Problem parameters.
  double const_f;

  WeakFormLinearLShape() : WeakForm(1) {
    add_matrix_form(new MatrixFormVolLShape(0, 0));
  };

private:
  class MatrixFormVolLShape : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolLShape(int i, int j) : WeakForm::MatrixFormVol(i, j) {
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

class DirichletFunctionBoundaryCondition : public DirichletBoundaryCondition {
public:
  DirichletFunctionBoundaryCondition(std::string marker, ExactSolutionLShape* exact_solution)
        : DirichletBoundaryCondition(Hermes::vector<std::string>()), exact_solution(exact_solution) {
    markers.push_back(marker);
  }

  ~DirichletFunctionBoundaryCondition() { }

  inline BoundaryConditionValueType get_value_type() const { return BoundaryCondition::BC_FUNCTION; }

  scalar function(double x, double y) const {
    return exact_solution->fn(x, y);
  }

  // Member.
  ExactSolutionLShape* exact_solution;
};




















// Bilinear form corresponding to the Laplace equation.
template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}
