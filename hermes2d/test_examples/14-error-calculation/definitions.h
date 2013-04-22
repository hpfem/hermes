#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

class CustomNormFormVol : public NormFormVol<double>
{
public:
  CustomNormFormVol(int i, int j);

  virtual double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const;
};

class CustomNormFormSurf : public NormFormSurf<double>
{
public:
  CustomNormFormSurf(int i, int j);

  virtual double value(int n, double *wt, Func<double> *u, Func<double> *v, Geom<double> *e) const;
};

class CustomNormFormDG : public NormFormDG<double>
{
public:
  CustomNormFormDG(int i, int j);

  virtual double value(int n, double *wt, DiscontinuousFunc<double> *u, DiscontinuousFunc<double> *v, Geom<double> *e) const;
};

/*
  This represents such a function.

<- mesh_size ->
  ---------
  | 2 | 0 |
  ---------
  | 0 | 2 |
  ---------
  */

class CustomExactSolutionScalar : public ExactSolutionScalar<double>
{
public:
  CustomExactSolutionScalar(MeshSharedPtr mesh);

  virtual double value(double x, double y) const;

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual Ord ord(Ord x, Ord y) const;

  ~CustomExactSolutionScalar();
	
	virtual MeshFunction<double>* clone() const;

  virtual void set_active_element(Element* e);
private:
  // For jumps detection.
  int element_id;
};