#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;
using namespace Hermes::Hermes2D::RefinementSelectors;

/* Custom function that is used in the exact solution and in right-hand side */

class CustomExactFunction1
{
public:
  CustomExactFunction1() {};

  double val(double x);
  
  double dx(double x);
  
  double ddxx(double x);
};

class CustomExactFunction2
{
public:
  CustomExactFunction2(double K) : K(K) {};

  double val(double x);
  
  double dx(double x);
  
  double ddxx(double x);

  double K;
};

/* Right-hand side */

class CustomRightHandSide1: public Hermes2DFunction<double>
{
public:
  CustomRightHandSide1(double K, double d_u, double sigma);

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  ~CustomRightHandSide1();

  CustomExactFunction1* cef1;
  CustomExactFunction2* cef2;
  double d_u, sigma;
};

class CustomRightHandSide2: public Hermes2DFunction<double>
{
public:
  CustomRightHandSide2(double K, double d_v);

  virtual double value(double x, double y) const;

  virtual Ord value(Ord x, Ord y) const;

  ~CustomRightHandSide2();

  CustomExactFunction1* cef1;
  CustomExactFunction2* cef2;
  double d_v;
};

/* Exact solution */

class ExactSolutionFitzHughNagumo1 : public ExactSolutionScalar<double>
{
public:
  ExactSolutionFitzHughNagumo1(MeshSharedPtr mesh);

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord(double x, double y) const;

  ~ExactSolutionFitzHughNagumo1();
	
	virtual MeshFunction<double>* clone() const;

  CustomExactFunction1* cef1;
};

class ExactSolutionFitzHughNagumo2 : public ExactSolutionScalar<double>
{
public:
  ExactSolutionFitzHughNagumo2(MeshSharedPtr mesh, double K);

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord(double x, double y) const;

	virtual MeshFunction<double>* clone() const;

  ~ExactSolutionFitzHughNagumo2();

  CustomExactFunction2* cef2;
	double K;
};

/* Weak forms */

class CustomResidual1 : public VectorFormVol<double>
{
public:
  CustomResidual1(double d_u, double sigma, CustomRightHandSide1* g1)
    : VectorFormVol<double>(0), d_u(d_u), sigma(sigma), g1(g1) {};

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                       Geom<double> *e, Func<double> **ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                  Geom<Ord> *e, Func<Ord> **ext) const;

  virtual VectorFormVol<double>* clone() const;

private:
  double d_u;
  double sigma;
  CustomRightHandSide1* g1;
};

class CustomResidual2 : public VectorFormVol<double>
{
public:
  CustomResidual2(double d_v, CustomRightHandSide2* g2)
    : VectorFormVol<double>(1), d_v(d_v), g2(g2) {};

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *v,
                       Geom<double> *e, Func<double> **ext) const;
  
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                  Geom<Ord> *e, Func<Ord> **ext) const;
  
  virtual VectorFormVol<double>* clone() const;
  
private:
  double d_v;
  CustomRightHandSide2* g2;
};

class CustomWeakForm : public WeakForm<double>
{
public:
  CustomWeakForm(CustomRightHandSide1* g1, CustomRightHandSide2* g2);
};
