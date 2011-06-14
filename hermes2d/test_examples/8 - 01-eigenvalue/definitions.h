#include "hermes2d.h"

/* Weak forms */

class WeakFormEigenLeft : public WeakForm
{
public:
  WeakFormEigenLeft();

private:
  class MatrixFormPotential : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormPotential(int i, int j) : WeakForm::MatrixFormVol(i, j) {};

    template<typename Real, typename Scalar>
    Scalar matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, 
                       Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const;

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const;
  };
};

class WeakFormEigenRight : public WeakForm
{
public:
  WeakFormEigenRight();
};

