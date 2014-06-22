#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

enum VelocityComponent
  {
    velX = 1,
    velY = 2
  };

class CouplingErrorFormVelocity : public MatrixFormVol<double>
{
public:
  CouplingErrorFormVelocity(VelocityComponent component, double c_p);

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *sln_i, Func<double> *sln_j, Geom<double> *e, Func<double>* *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *sln_i, Func<Ord> *sln_j, Geom<Ord> *e, Func<Ord>* *ext) const;

  MatrixFormVol<double>* clone() const;

  VelocityComponent component;
  double c_p;
};

class CouplingErrorFormTemperature : public MatrixFormVol<double>
{
public:
  CouplingErrorFormTemperature(VelocityComponent component, double c_p);

  virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *sln_i, Func<double> *sln_j, Geom<double> *e, Func<double>* *ext) const;

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *sln_i, Func<Ord> *sln_j, Geom<Ord> *e, Func<Ord>* *ext) const;

  MatrixFormVol<double>* clone() const;

  VelocityComponent component;
  double c_p;
};

