#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

// Exact solution.
#include "exact_solution.cpp"

class WeakFormNIST06 : public WeakForm
{
public:
  WeakFormNIST06(ExactSolutionNIST06* exact_solution) : WeakForm(1) {
    add_matrix_form(new MatrixFormVolNIST06(0, 0, exact_solution->epsilon));
    
    add_vector_form(new VectorFormVolNIST06(0, exact_solution));
  };

private:
  class MatrixFormVolNIST06 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolNIST06(int i, int j, double epsilon) : WeakForm::MatrixFormVol(i, j), epsilon(epsilon) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar val = 0;
      for (int i=0; i < n; i++) {
        val = val + wt[i] * epsilon * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
        val = val + wt[i] * (2*u->dx[i] + u->dy[i]) * v->val[i];
      }

      return val;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    double epsilon;
  };

  class VectorFormVolNIST06 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolNIST06(int i, ExactSolutionNIST06* exact_solution) : WeakForm::VectorFormVol(i), exact_solution(exact_solution) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (exact_solution->rhs<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(30);
    }

    // Members.
    ExactSolutionNIST06* exact_solution;
  };
};

class DirichletFunctionBoundaryConditionExact : public DirichletBoundaryCondition
{
public:
  DirichletFunctionBoundaryConditionExact(std::string marker, ExactSolutionNIST06* exact_solution) : 
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

  ExactSolutionNIST06* exact_solution;
};