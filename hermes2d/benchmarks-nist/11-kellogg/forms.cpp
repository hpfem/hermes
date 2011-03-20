#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"

// Exact solution.
#include "exact_solution.cpp"

class WeakFormNIST11 : public WeakForm
{
public:
  WeakFormNIST11(double r) : WeakForm(1)
  {
    add_matrix_form(new MatrixFormVol_I_III(0, 0, r));
    
    add_matrix_form(new MatrixFormVol_II_IV(0, 0));
  };

private:
  class MatrixFormVol_I_III : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVol_I_III(int i, int j, double r) : WeakForm::MatrixFormVol(i, j), r(r) { 
      sym = HERMES_SYM;
      area = "0";
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      return r * int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Member.
    double r;
  };

  class MatrixFormVol_II_IV : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVol_II_IV(int i, int j) : WeakForm::MatrixFormVol(i, j) {
      sym = HERMES_SYM;
      area = "1";
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

class EssentialBCNonConstantExact : public EssentialBC
{
public:
  EssentialBCNonConstantExact(std::string marker, ExactSolutionNIST11* exact_solution) : 
        EssentialBC(Hermes::vector<std::string>()), exact_solution(exact_solution) 
  {
    markers.push_back(marker);
  };
  
  ~EssentialBCNonConstantExact() {};

  virtual BoundaryConditionValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar function(double x, double y) const {
    return exact_solution->fn(x, y);
  };

  ExactSolutionNIST11* exact_solution;
};