#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

class WeakFormPoisson : public WeakForm
{
public:
  WeakFormPoisson(double slope) : WeakForm(1) {
    add_matrix_form(new MatrixFormVolPoisson(0, 0));
    add_vector_form(new VectorFormVolPoisson(0, slope));
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

  class VectorFormVolPoisson : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolPoisson(int i, double slope) : WeakForm::VectorFormVol(i), 
          slope(slope) { }

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
};

// Exact solution (needed in the Dirichlet condition).
#include "exact_solution.cpp"

class DirichletFunctionBoundaryCondition : public DirichletBoundaryCondition {
public:
  DirichletFunctionBoundaryCondition(std::string marker, ExactSolutionPoisson* exact_solution) : DirichletBoundaryCondition(Hermes::vector<std::string>()),
    exact_solution(exact_solution) {
      markers.push_back(marker);
  }

  ~DirichletFunctionBoundaryCondition() { }

  inline BoundaryConditionValueType get_value_type() const { return BoundaryCondition::BC_FUNCTION; }

  scalar function(double x, double y) const {
    return exact_solution->fn(x, y);
  }

  // Member.
  ExactSolutionPoisson* exact_solution;
};
