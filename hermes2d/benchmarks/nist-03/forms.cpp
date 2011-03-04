#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

// Exact solution.
#include "exact_solution.cpp"

class WeakFormNIST03 : public WeakForm
{
public:
  WeakFormNIST03(ExactSolutionNIST03U* exact_solution_u, ExactSolutionNIST03V* exact_solution_v) : WeakForm(2)
  {
    add_matrix_form(new MatrixFormVolNIST03_0_0(exact_solution_u->A, exact_solution_u->B));
    add_matrix_form(new MatrixFormVolNIST03_0_1(exact_solution_u->C));
    add_matrix_form(new MatrixFormVolNIST03_1_1(exact_solution_u->A, exact_solution_u->B));

    add_vector_form(new VectorFormVolNIST03_0(exact_solution_u));
    add_vector_form(new VectorFormVolNIST03_1(exact_solution_v));
  }

private:
  class MatrixFormVolNIST03_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolNIST03_0_0(double A, double B) : WeakForm::MatrixFormVol(0, 0), A(A), B(B) { }

template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (A * u->dx[i] * v->dx[i] + B * u->dy[i] * v->dy[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
       return Ord(30);
    }

    // Members.
    double A, B;
  };

  class MatrixFormVolNIST03_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolNIST03_0_1(double C) : WeakForm::MatrixFormVol(0, 1), C(C) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (C * u->dx[i] * v->dy[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
       return Ord(30);
    }

    // Member.
    double C;
  };

  class MatrixFormVolNIST03_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolNIST03_1_1(double A, double B) : WeakForm::MatrixFormVol(1, 1), A(A), B(B) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
    {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (B * u->dx[i] * v->dx[i] + A * u->dy[i] * v->dy[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext)
    {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
    {
       return Ord(30);
    }
  
    // Members.
    double A, B;
  };

  class VectorFormVolNIST03_0 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolNIST03_0(ExactSolutionNIST03U* exact_solution) : WeakForm::VectorFormVol(0), 
    exact_solution(exact_solution) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (rhs_0<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(30);
    }

    template<typename Real>
    Real rhs_0(Real x, Real y)
    {
      return exact_solution->A * exact_solution->dudxdudx(x, y) + exact_solution->B * exact_solution->dudydudy(x, y) + exact_solution->C * exact_solution->dvdxdvdy(x, y);
    }

    // Member.
    ExactSolutionNIST03U* exact_solution;
  };

  class VectorFormVolNIST03_1 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolNIST03_1(ExactSolutionNIST03V* exact_solution) : WeakForm::VectorFormVol(1), 
    exact_solution(exact_solution) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (rhs_1<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(30);
    }

    template<typename Real>
    Real rhs_1(Real x, Real y)
    {
      return exact_solution->B * exact_solution->dvdxdvdx(x, y) + exact_solution->A * exact_solution->dvdydvdy(x, y) + exact_solution->C * exact_solution->dudxdudy(x, y);
    }

    // Member.
    ExactSolutionNIST03V* exact_solution;
  };
};

class DirichletFunctionBoundaryConditionExactU : public DirichletBoundaryCondition
{
public:
  DirichletFunctionBoundaryConditionExactU(std::string marker, ExactSolutionNIST03U* exact_solution) : 
        DirichletBoundaryCondition(Hermes::vector<std::string>()), exact_solution(exact_solution) 
  {
    markers.push_back(marker);
  }
  
  ~DirichletFunctionBoundaryConditionExactU() {}

  virtual BoundaryConditionValueType get_value_type() const { 
    return BC_FUNCTION; 
  }

  virtual scalar function(double x, double y) const {
    return exact_solution->u_fn(x, y);
  }

  ExactSolutionNIST03U* exact_solution;
};

class DirichletFunctionBoundaryConditionExactV : public DirichletBoundaryCondition
{
public:
  DirichletFunctionBoundaryConditionExactV(std::string marker, ExactSolutionNIST03V* exact_solution) : 
        DirichletBoundaryCondition(Hermes::vector<std::string>()), exact_solution(exact_solution) 
  {
    markers.push_back(marker);
  }
  
  ~DirichletFunctionBoundaryConditionExactV() {}

  virtual BoundaryConditionValueType get_value_type() const { 
    return BC_FUNCTION; 
  }

  virtual scalar function(double x, double y) const {
    return exact_solution->v_fn(x, y);
  }

  ExactSolutionNIST03V* exact_solution;
};