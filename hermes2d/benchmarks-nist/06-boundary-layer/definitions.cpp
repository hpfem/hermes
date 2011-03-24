#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"

/* Right-hand side */

class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide(double epsilon) : DefaultNonConstRightHandSide(), epsilon(epsilon) {};

  virtual double value(double x, double y) const{
    return -epsilon*(-2*pow(M_PI,2)*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*cos(M_PI*(x + y)) 
           + 2*M_PI*(1 - exp(-(1 - x)/epsilon))*exp(-(1 - y)/epsilon)*sin(M_PI*(x + y))/epsilon 
           + 2*M_PI*(1 - exp(-(1 - y)/epsilon))*exp(-(1 - x)/epsilon)*sin(M_PI*(x + y))/epsilon 
           - (1 - exp(-(1 - y)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - x)/epsilon)/pow(epsilon,2) 
           - (1 - exp(-(1 - x)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - y)/epsilon)/pow(epsilon,2)) 
           - 3*M_PI*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*sin(M_PI*(x + y)) 
           - 2*(1 - exp(-(1 - y)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - x)/epsilon)/epsilon 
           - (1 - exp(-(1 - x)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - y)/epsilon)/epsilon;
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(8);
  }

  double epsilon;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double epsilon) 
        : ExactSolutionScalar(mesh), epsilon(epsilon) {};

  virtual scalar value(double x, double y) const {
    return (1 - exp(-(1-x)/epsilon)) * (1 - exp(-(1-y)/epsilon)) * cos(M_PI * (x + y)); 
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = -M_PI*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*sin(M_PI*(x + y)) 
         - (1 - exp(-(1 - y)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - x)/epsilon)/epsilon;
    dy = -M_PI*(1 - exp(-(1 - x)/epsilon))*(1 - exp(-(1 - y)/epsilon))*sin(M_PI*(x + y)) 
         - (1 - exp(-(1 - x)/epsilon))*cos(M_PI*(x + y))*exp(-(1 - y)/epsilon)/epsilon;
  };

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(8);
  }

  scalar epsilon;
};

/* Weak forms */

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(CustomRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new CustomMatrixFormVol(0, 0, rhs->epsilon));
    add_vector_form(new DefaultVectorFormVolNonConst(0, rhs));
  };

private:
  class CustomMatrixFormVol : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVol(int i, int j, double epsilon) 
          : WeakForm::MatrixFormVol(i, j), epsilon(epsilon) { }

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) {
      Scalar val = 0;
      for (int i=0; i < n; i++) {
        val = val + wt[i] * epsilon * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
        val = val + wt[i] * (2*u->dx[i] + u->dy[i]) * v->val[i];
      }

      return val;
    }

    scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) {
      return matrix_form<scalar, scalar>(n, wt, u_ext, u, v, e, ext);
    }

    Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
    }

    double epsilon;
  };
};

