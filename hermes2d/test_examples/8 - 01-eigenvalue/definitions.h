#include "hermes2d.h"

/* Weak forms */

class WeakFormEigenLeft : public WeakForm<double>
{
public:
  WeakFormEigenLeft();

private:
  class MatrixFormPotential : public MatrixFormVol<double>
  {
  public:
    MatrixFormPotential(int i, int j) : MatrixFormVol<double>(i, j) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<double> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
  };
};

class WeakFormEigenRight : public WeakForm<double>
{
public:
  WeakFormEigenRight();
};

