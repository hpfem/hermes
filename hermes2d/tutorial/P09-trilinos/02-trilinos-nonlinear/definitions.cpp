#include "weakform/weakform.h"
#include "integrals/h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/h1.h"

/* Right-hand side */

class CustomRightHandSide: public DefaultFunction
{
public:
  CustomRightHandSide()
    : DefaultFunction() {};

  virtual double value(double x, double y) const {
    return - kx(x,y) * dudx(x,y) - ky(x,y) * dudy(x,y) - k(x,y) * (dudxx(x,y) + dudyy(x,y));
  }

  virtual Ord ord(Ord x, Ord y) const {
    return - kx(x,y) * dudx(x,y) - ky(x,y) * dudy(x,y) - k(x,y) * (dudxx(x,y) + dudyy(x,y));
  }

  template<typename Real>
  Real u(Real x, Real y) const {  return (x - x*x) * (y - y*y);  }
  template<typename Real>
  Real dudx(Real x, Real y) const {  return (1- 2*x) * y * (1 - y);  }
  template<typename Real>
  Real dudy(Real x, Real y) const {  return (1- 2*y) * x * (1 - x);  }

  template<typename Real>
  Real dudxx(Real x, Real y) const {  return -2.0 * (y-y*y);  }
  template<typename Real>
  Real dudyy(Real x, Real y) const {  return -2.0 * (x-x*x);  }
  template<typename Real>
  Real dudxy(Real x, Real y) const {  return (1- 2*y) * (1 - 2*x);  }

  template<typename Real>
  Real k(Real x, Real y) const {  return 1.0 / sqrt(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)));  }
  template<typename Real>
  Real kx(Real x, Real y) const {  return -0.5 * pow(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)), -1.5) *
                   (2.0 * dudx(x,y) * dudxx(x,y) + 2.0 * dudy(x,y) * dudxy(x,y));  }
  template<typename Real>
  Real ky(Real x, Real y) const {  return -0.5 * pow(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)), -1.5) *
                   (2.0 * dudx(x,y) * dudxy(x,y) + 2.0 * dudy(x,y) * dudyy(x,y));  }
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh)
            : ExactSolutionScalar(mesh) {};

  double value(double x, double y) const {
    return  x * y * (1-x) * (1-y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = (1- 2*x) * y * (1 - y);
	  dy = (1- 2*y) * x * (1 - x);
	  
  };

  virtual Ord ord(Ord x, Ord y) const {
    return (1- 2*x) * y * (1 - y);
  }
};

/* Initial solution */

class CustomInitialSolution : public ExactSolutionScalar
{
public:
  CustomInitialSolution(Mesh* mesh)
            : ExactSolutionScalar(mesh) {};

  double value(double x, double y) const {
    return  0;
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = 0;
	  dy = 0;
	};

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(1);
  }
};

/* Weak form */
class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(bool JFNK = false, bool precondition_jacobian = false, bool precondition_jacobian_approx = false) : WeakForm(1, JFNK)
  {
    // Jacobian forms - volumetric.
    if(!JFNK || precondition_jacobian)
      add_matrix_form(new JacobianFormVol(0, 0));

    // Residual forms - volumetric.
    add_vector_form(new ResidualFormVol(0, new CustomRightHandSide()));

    // Preconditioning form.
    if(JFNK && precondition_jacobian_approx)
      add_matrix_form(new PrecondFormVol(0, 0));
  }

  ~CustomWeakForm() {}

private:
  class JacobianFormVol : public WeakForm::MatrixFormVol
  {
  public:
    JacobianFormVol(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * ( -0.5 * pow(1.0 + sqr(u_ext[0]->dx[i]) + sqr(u_ext[0]->dy[i]), -1.5) * 
                           (2.0 * u_ext[0]->dx[i] * u->dx[i] + 2.0 * u_ext[0]->dx[i] * u->dx[i])
                           * (u_ext[0]->dx[i] * v->dx[i] + u_ext[0]->dy[i] * v->dy[i]) +
                           (pow(1.0 + sqr(u_ext[0]->dx[i]) + sqr(u_ext[0]->dy[i]), -0.5))
                           * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]) );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }
  };

  class ResidualFormVol : public WeakForm::VectorFormVol
  {
  public:
    ResidualFormVol(int i, CustomRightHandSide* rhs) : WeakForm::VectorFormVol(i), rhs(rhs) {}

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const {
      Func<scalar>* u = u_ext[0];
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * ((pow(1.0 + sqr(u->dx[i]) + sqr(u->dy[i]), -0.5)) * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i])
                           - rhs->value(e->x[i], e->y[i]) * v->val[i] );
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the test function and solution plus two.
      return Ord(10);
    }

  private:
    // Problem parameters.
    CustomRightHandSide* rhs;
  };

  class PrecondFormVol : public WeakForm::MatrixFormVol
  {
  public:
    PrecondFormVol(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM) { }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const {
      // Returning the sum of the degrees of the basis and test function plus two.
      return Ord(10);
    }
  };
};