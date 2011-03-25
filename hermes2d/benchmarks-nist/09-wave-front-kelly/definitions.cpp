#include "weakform/weakform.h"
#include "integrals/integrals_h1.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/laplace.h"
#include "adapt/kelly_type_adapt.h"

/* Right-hand side */

class CustomRightHandSide: public DefaultNonConstRightHandSide
{
public:
  CustomRightHandSide(double alpha, double x_loc, double y_loc, double r_zero) 
    : DefaultNonConstRightHandSide(), alpha(alpha), x_loc(x_loc), y_loc(y_loc), r_zero(r_zero) { };

  virtual double value(double x, double y) const {  
    double a = pow(x - x_loc, 2);
    double b = pow(y - y_loc, 2);
    double c = sqrt(a + b);
    double d = ((alpha*x - (alpha * x_loc)) * (2*x - (2 * x_loc)));
    double e = ((alpha*y - (alpha * y_loc)) * (2*y - (2 * y_loc)));
    double f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);
    double g = (alpha * c - (alpha * r_zero));

    return -(((alpha/(c * f)) - (d/(2 * pow(a + b, 1.5) * f)) 
           - ((alpha * d * g)/((a + b) * pow(f, 2))) + 
	      (alpha/(c * f)) - (e/(2 * pow(a + b, 1.5) * f)) 
           - ((alpha * e * g)/((a + b) * pow(f, 2)))));
  }

  virtual Ord ord (Ord x, Ord y) const {
    return Ord(8);  
  }
  double alpha, x_loc, y_loc, r_zero;
};

/* Exact solution */

class CustomExactSolution : public ExactSolutionScalar
{
public:
  CustomExactSolution(Mesh* mesh, double alpha, double x_loc, double 
                      y_loc, double r_zero) 
             : ExactSolutionScalar(mesh), alpha(alpha), x_loc(x_loc), 
                                   y_loc(y_loc), r_zero(r_zero) { }

  virtual scalar value(double x, double y) const {
    return atan(alpha * (sqrt(pow(x - x_loc, 2) + pow(y - y_loc, 2)) - r_zero));
  };

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    double a = pow(x - x_loc, 2);
    double b = pow(y - y_loc, 2);
    double c = sqrt(a + b);
    double d = (alpha*x - (alpha * x_loc));
    double e = (alpha*y - (alpha * y_loc));
    double f = (pow(alpha*c - (alpha * r_zero), 2) + 1.0);

    dx = (d/(c * f));
    dy = (e/(c * f));
  };

  Ord ord (Ord x, Ord y) const {
    return Ord(8);  
  }

  double alpha, x_loc, y_loc, r_zero;
};


/* Weak forms */

class CustomWeakFormPoisson : public WeakForm
{
public:
  CustomWeakFormPoisson(DefaultNonConstRightHandSide* rhs) : WeakForm(1) {
    add_matrix_form(new DefaultMatrixFormStiffness(0, 0));
    add_vector_form(new DefaultVectorFormVolNonConst(0, rhs));
  };
};

/* Bilinear form inducing the energy norm */

class EnergyErrorForm : public Adapt::MatrixFormVolError
{
public:
  EnergyErrorForm(WeakForm* problem_wf) : MatrixFormVolError(HERMES_UNSET_NORM)
  {
    this->form = problem_wf->mfvol[0];
  }
  
  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
                        Func<scalar> *u, Func<scalar> *v, Geom<double> *e,
                        ExtData<scalar> *ext) const 
  {
    return this->form->value(n, wt, u_ext, u, v, e, ext);
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                  Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e,
                  ExtData<Ord> *ext) const 
  {
    return this->form->ord(n, wt, u_ext, u, v, e, ext);
  }
    
private:
  WeakForm::MatrixFormVol* form;
};

/* Linear form for the residual error estimator */

class ResidualErrorForm : public KellyTypeAdapt::ErrorEstimatorForm
{
public:
  ResidualErrorForm(CustomRightHandSide* rhs) : ErrorEstimatorForm(0), rhs(rhs) {};
  
  scalar residual_estimator(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Geom<double> *e, ExtData<scalar> *ext) const
  {
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
    scalar result = 0.;
    
    for (int i = 0; i < n; i++)
      result += wt[i] * sqr(rhs->value(e->x[i], e->y[i]) - u->laplace[i] );
    
    return result * sqr(e->diam);
#else
    error("Define H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h if you want to use second derivatives of shape functions in weak forms.");
#endif
  }
  
  Ord residual_estimator(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Geom<Ord> *e, ExtData<Ord> *ext) const
  {
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
    Ord result = 0.;
    
    for (int i = 0; i < n; i++)
      result += wt[i] * sqr(rhs->ord(e->x[i], e->y[i]) - u->laplace[i] );
    
    return result * sqr(e->diam);
#else
    error("Define H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h if you want to use second derivatives of shape functions in weak forms.");
#endif
  }

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
              Func<scalar> *u, Geom<double> *e,
              ExtData<scalar> *ext) const 
  {
    return residual_estimator(n, wt, u_ext, u, e, ext);
  }
  
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                    Func<Ord> *u, Geom<Ord> *e,
                    ExtData<Ord> *ext) const 
  {
    return residual_estimator(n, wt, u_ext, u, e, ext);
  }
  
private:  
  CustomRightHandSide* rhs;
  
};

// Linear form for the interface error estimator.
class InterfaceErrorForm : public KellyTypeAdapt::ErrorEstimatorForm
{
public:
  InterfaceErrorForm() : ErrorEstimatorForm(0, H2D_DG_INNER_EDGE) {};
  
  template<typename Real, typename Scalar>
  Scalar interface_estimator(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Geom<Real> *e, ExtData<Scalar> *ext)
  {
    Scalar result = 0.;
    for (int i = 0; i < n; i++)
      result += wt[i] * sqr( e->nx[i] * (u->get_dx_central(i) - u->get_dx_neighbor(i)) +
                             e->ny[i] * (u->get_dy_central(i) - u->get_dy_neighbor(i))  );
    return result * e->diam / 24.;
  }

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],
              Func<scalar> *u, Geom<double> *e,
              ExtData<scalar> *ext)
  {
    return interface_estimator<double, scalar>(n, wt, u_ext, u, e, ext);
  }
  
  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[],
                    Func<Ord> *u, Geom<Ord> *e,
                    ExtData<Ord> *ext)
  {
    return interface_estimator<Ord, Ord>(n, wt, u_ext, u, e, ext);
  }  
};
