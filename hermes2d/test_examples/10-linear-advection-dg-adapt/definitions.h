#include "hermes2d.h"

/* Namespaces used */

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh), t(0.0) {};
  ~CustomExactSolution();

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual double value(double x, double y) const;

  virtual Ord ord(double x, double y) const;

  MeshFunction<double>* clone() const;

  void setTime(double time);

private:
  double t;
};

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(MeshFunctionSharedPtr<double> previousSln, MeshSharedPtr mesh);

  WeakForm<double>* clone() const;

private:
  MeshFunctionSharedPtr<double> previousSln;
  class CustomMatrixFormVol : public MatrixFormVol<double>
  {
  public:
    CustomMatrixFormVol(int i, int j) : MatrixFormVol<double>(i, j) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, GeomVol<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, GeomVol<Ord> *e, Func<Ord> **ext) const;

    MatrixFormVol<double>* clone() const;
  };

  class CustomVectorFormVol : public VectorFormVol<double>
  {
  public:
    CustomVectorFormVol(int i) : VectorFormVol<double>(i) {};

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, GeomVol<Real> *e, Func<Scalar> **ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomVol<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomVol<Ord> *e, Func<Ord> **ext) const;

    VectorFormVol<double>* clone() const;
  };

  class CustomVectorFormInterface : public VectorFormDG<double>
  {
  public:
    CustomVectorFormInterface(int i) : VectorFormDG<double>(i)
    {
    };

    template<typename Real, typename Scalar>
    Scalar vector_form(int n, double *wt, DiscontinuousFunc<Scalar> **u_ext, Func<Real> *v, InterfaceGeom<Real> *e, DiscontinuousFunc<Scalar> **ext) const;

    virtual double value(int n, double *wt, DiscontinuousFunc<double> **u_ext, Func<double> *v, InterfaceGeom<double> *e, DiscontinuousFunc<double> **ext) const;

    virtual Ord ord(int n, double *wt, DiscontinuousFunc<Ord> **u_ext, Func<Ord> *v, InterfaceGeom<Ord> *e, DiscontinuousFunc<Ord> **ext) const;

    VectorFormDG<double>* clone() const;

  };

  class CustomVectorFormSurface : public VectorFormSurf<double>
  {
  public:
    CustomVectorFormSurface(int i) : VectorFormSurf<double>(i)  {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v, GeomSurf<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, GeomSurf<Ord> *e, Func<Ord> **ext) const;

    VectorFormSurf<double>* clone() const;

    template<typename Real>
    Real F(Real x, Real y) const;
  };
  
  double laxF_flux(double u_cent, double u_neib, double a_dot_n) const;

  MeshSharedPtr mesh;
};
