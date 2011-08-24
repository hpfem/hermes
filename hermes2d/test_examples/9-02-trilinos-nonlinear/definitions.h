#include "hermes2d.h"

/* Right-hand side */

class CustomRightHandSide: public HermesFunction
{
public:
  CustomRightHandSide() : HermesFunction() {};

  virtual double value(double x, double y) const;

  virtual Ord ord(Ord x, Ord y) const;

  template<typename Real>
  Real u(Real x, Real y) const;
  template<typename Real>
  Real dudx(Real x, Real y) const;
  template<typename Real>
  Real dudy(Real x, Real y) const;

  template<typename Real>
  Real dudxx(Real x, Real y) const;
  template<typename Real>
  Real dudyy(Real x, Real y) const;
  template<typename Real>
  Real dudxy(Real x, Real y) const;

  template<typename Real>
  Real k(Real x, Real y) const;
  template<typename Real>
  Real kx(Real x, Real y) const;
  template<typename Real>
  Real ky(Real x, Real y) const;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh)
             : ExactSolutionScalar(mesh) {};

  double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
};

/* Initial solution */

class CustomInitialSolution : public ExactSolutionScalar
{
public:
  CustomInitialSolution(Mesh* mesh)
               : ExactSolutionScalar(mesh) {};

  double value(double x, double y) const;

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const;

  virtual Ord ord(Ord x, Ord y) const;
};

/* Weak form */
class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(bool JFNK = false, bool precondition_jacobian = false, bool precondition_jacobian_approx = false);

  ~CustomWeakForm() {}

private:
  class JacobianFormVol : public WeakForm::MatrixFormVol
  {
  public:
    JacobianFormVol(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const;
  };

  class ResidualFormVol : public WeakForm::VectorFormVol
  {
  public:
    ResidualFormVol(int i, CustomRightHandSide* rhs) : WeakForm::VectorFormVol(i), rhs(rhs) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    // Problem parameters.
    CustomRightHandSide* rhs;
  };

  class PrecondFormVol : public WeakForm::MatrixFormVol
  {
  public:
    PrecondFormVol(int i, int j) : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM) {};

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
            Geom<Ord> *e, ExtData<Ord> *ext) const;
  };
};

