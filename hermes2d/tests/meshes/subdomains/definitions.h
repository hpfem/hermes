#include "hermes2d.h"

/* Weak forms */

class CustomVectorFormVol : public Hermes::Hermes2D::VectorFormVol<double>
{
public:
  CustomVectorFormVol(int i);

  virtual double value(int n, double *wt, Hermes::Hermes2D::Func<double> *u_ext[], Hermes::Hermes2D::Func<double> *v,
    Hermes::Hermes2D::Geom<double> *e, Hermes::Hermes2D::ExtData<double> *ext) const;

  virtual Hermes::Ord ord(int n, double *wt, Hermes::Hermes2D::Func<Hermes::Ord> *u_ext[], Hermes::Hermes2D::Func<Hermes::Ord> *v,
    Hermes::Hermes2D::Geom<Hermes::Ord> *e, Hermes::Hermes2D::ExtData<Hermes::Ord> *ext) const;
};


class CustomWeakFormPoisson : public Hermes::Hermes2D::WeakForm<double>
{
public:
  CustomWeakFormPoisson();
};

