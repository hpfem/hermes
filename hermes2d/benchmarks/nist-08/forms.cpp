#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/boundaryconditions.h"

// Exact solution.
#include "exact_solution.cpp"

class WeakFormNIST08 : public WeakForm
{
public:
  WeakFormNIST08(double alpha) : WeakForm(1) {
    add_matrix_form(new MatrixFormVolNIST08(0, 0, alpha));
    
    add_vector_form(new VectorFormVolNIST08(0, alpha));
  };

private:
  class MatrixFormVolNIST08 : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolNIST08(int i, int j, double alpha) : WeakForm::MatrixFormVol(i, j), alpha(alpha) {
      sym = HERMES_SYM;
    }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar val = 0;
      for (int i=0; i < n; i++) {
        Scalar x = e->x[i];
        Scalar y = e->y[i];
        Scalar r = sqrt(x*x + y*y);
        Scalar h = 1/(alpha + r);
        Scalar grad_u_grad_v = u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i];
        val += wt[i] * (grad_u_grad_v - pow(h, 4) * u->val[i] * v->val[i]);
      }

      return val;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    // Member.
    double alpha;
  };

  class VectorFormVolNIST08 : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolNIST08(int i, double alpha) : WeakForm::VectorFormVol(i), alpha(alpha) { }

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (rhs<Real>(e->x[i], e->y[i]) * v->val[i]);
      return result;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return vector_form<scalar, scalar>(n, wt, u_ext, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) {
      return Ord(30);
    }

    template<typename Real>
    Real rhs(Real x, Real y) {
      return -sin(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/pow((alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4) 
             + 2*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(1.0/2.0))) 
             + pow(x,2)*sin(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4)*(pow(x,2) + pow(y,2))) 
             + pow(y,2)*sin(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),4)*(pow(x,2) + pow(y,2))) 
             - pow(x,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(3.0/2.0))) 
             - pow(y,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),2)*pow((pow(x,2) + pow(y,2)),(3.0/2.0))) 
             - 2*pow(x,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),3)*(pow(x,2) + pow(y,2))) 
             - 2*pow(y,2)*cos(1.0/(alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))))/(pow((alpha + pow((pow(x,2) + pow(y,2)),(1.0/2.0))),3)*(pow(x,2) + pow(y,2)));
    }
    
    // Member.
    double alpha;
  };
};

class DirichletNonConstantExact : public DirichletBoundaryCondition
{
public:
  DirichletNonConstantExact(std::string marker, ExactSolutionNIST08* exact_solution) : 
        DirichletBoundaryCondition(Hermes::vector<std::string>()), exact_solution(exact_solution) 
  {
    markers.push_back(marker);
  };
  
  ~DirichletNonConstantExact() {};

  virtual BoundaryConditionValueType get_value_type() const { 
    return BC_FUNCTION; 
  };

  virtual scalar function(double x, double y) const {
    return exact_solution->fn(x, y);
  };

  ExactSolutionNIST08* exact_solution;
};