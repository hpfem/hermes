#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::WeakFormsH1;
using namespace Hermes::Hermes2D::Views;

/* Nonlinearity lambda(u) = Hermes::pow(u, alpha) */

class CustomNonlinearity : public Hermes1DFunction<double>
{
public:
  CustomNonlinearity(double alpha);

  virtual double value(double u) const;

  virtual Ord value(Ord u) const;

  protected:
    double alpha;
};

/* Weak forms */

// NOTE: The linear problem in the Picard's method is 
//       solved using the Newton's method.

class CustomWeakFormPicard : public WeakForm<double>
{
public:
  CustomWeakFormPicard(MeshFunctionSharedPtr<double> prev_iter_sln, Hermes1DFunction<double>* lambda, Hermes2DFunction<double>* f);

private:
  class CustomJacobian : public MatrixFormVol<double>
  {
  public:
    CustomJacobian(int i, int j, Hermes1DFunction<double>* lambda) : MatrixFormVol<double>(i, j), lambda(lambda) {};

    virtual double value(int n, double *wt, Func<double> *u_ext[], Func<double> *u,
                         Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                    Geom<Ord> *e, Func<Ord> **ext) const;
    
    MatrixFormVol<double>* clone() const;

    protected:
      Hermes1DFunction<double>* lambda;
  };

  class CustomResidual : public VectorFormVol<double>
  {
  public:
    CustomResidual(int i, Hermes1DFunction<double>* lambda, Hermes2DFunction<double>* f) 
      : VectorFormVol<double>(i), lambda(lambda), f(f) 
    {
    }

    virtual double value(int n, double *wt, Func<double> *u_ext[],
                         Func<double> *v, Geom<double> *e, Func<double> **ext) const;

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                    Geom<Ord> *e, Func<Ord> **ext) const;

    VectorFormVol<double>* clone() const;

  private:
      Hermes1DFunction<double>* lambda;
      Hermes2DFunction<double>* f;
  };
};

/* Essential boundary conditions */

class CustomEssentialBCNonConst : public EssentialBoundaryCondition<double> 
{
public:
  CustomEssentialBCNonConst(std::string marker) 
           : EssentialBoundaryCondition<double>(Hermes::vector<std::string>()) 
  {
    this->markers.push_back(marker);
  };

  virtual EssentialBoundaryCondition<double>::EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y, double n_x, double n_y, 
                       double t_x, double t_y) const;
};
