#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Preconditioners;
using namespace Hermes::Hermes2D;

class CustomRightHandSide : public Hermes2DFunction<double>
{
public:
  CustomRightHandSide() : Hermes2DFunction<double>() {};

  virtual double value(double x, double y) const;

  virtual Hermes::Ord ord(Ord x, Ord y) const;

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

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};

  double value(double x, double y) const;

  virtual void derivatives (double x, double y, double& dx, double& dy) const;

  virtual Ord ord(double x, double y) const;

  MeshFunction<double>* clone() const { return new CustomExactSolution(mesh); }
};

/* Weak form */
class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(bool JFNK = false, bool precondition_jacobian = false, bool precondition_jacobian_approx = false);

  ~CustomWeakForm() {}

private:
  class JacobianFormVol : public MatrixFormVol<double>
  {
  public:
    JacobianFormVol(int i, int j) : MatrixFormVol<double>(i, j) 
    {
      this->setSymFlag(HERMES_SYM);
    };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;

    virtual MatrixFormVol<double>* clone() const { return new JacobianFormVol(i, j); }
  };

  class ResidualFormVol : public VectorFormVol<double>
  {
  public:
    ResidualFormVol(int i, CustomRightHandSide* rhs) : VectorFormVol<double>(i), rhs(rhs) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
      GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;

    virtual VectorFormVol<double>* clone() const { return new ResidualFormVol(i, rhs); }
  private:
    // Problem parameters.
    CustomRightHandSide* rhs;
  };

  class PrecondFormVol : public MatrixFormVol<double>
  {
  public:
    PrecondFormVol(int i, int j) : MatrixFormVol<double>(i, j) 
    {
      this->setSymFlag(HERMES_SYM);
    };

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
      Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
      GeomVol<Ord> *e, Func<Ord> **ext) const;

    virtual MatrixFormVol<double>* clone() const { return new PrecondFormVol(i, j); }
  };
};

