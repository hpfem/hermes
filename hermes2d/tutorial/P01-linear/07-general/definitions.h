#include "hermes2d.h"

/* Essential boundary conditions */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition {
public:
  CustomEssentialBCNonConst(std::string marker);

  inline EssentialBCValueType get_value_type() const;

  virtual scalar value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const;
};

/* Weak forms */

class CustomWeakFormGeneral : public WeakForm
{
public:
  CustomWeakFormGeneral(std::string bdy_vertical);

private:
  class MatrixFormVolGeneral : public WeakForm::MatrixFormVol
  {
  public:
    MatrixFormVolGeneral(int i, int j);

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;
  };

  class VectorFormVolGeneral : public WeakForm::VectorFormVol
  {
  public:
    VectorFormVolGeneral(int i);

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double rhs(double x, double y) const;
  };

  class VectorFormSurfGeneral : public WeakForm::VectorFormSurf
  {
  public:
    VectorFormSurfGeneral(int i, std::string area = HERMES_ANY);

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                         Geom<double> *e, ExtData<scalar> *ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const;

  private:
    double g_N(double x, double y) const;
  };
};
