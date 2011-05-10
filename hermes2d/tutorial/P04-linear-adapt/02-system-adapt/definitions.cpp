#include "hermes2d.h"

using namespace WeakFormsH1;

/* Custom function that is used in the exact solution and in right-hand side */

class CustomExactFunction1
{
public:
  CustomExactFunction1() { };

  double val(double x) {
    return cos(M_PI*x/2);
  }
  double dx(double x) {
    return -sin(M_PI*x/2)*(M_PI/2.);
  }
  double ddxx(double x) {
    return -cos(M_PI*x/2)*(M_PI/2.)*(M_PI/2.);
  }
};

class CustomExactFunction2
{
public:
  CustomExactFunction2(double K) : K(K) {};

  double val(double x) {
    return 1. - (exp(K*x) + exp(-K*x))/(exp(K) + exp(-K));
  }
  double dx(double x) {
    return -K*(exp(K*x) - exp(-K*x))/(exp(K) + exp(-K));
  }
  double ddxx(double x) {
    return -K*K*(exp(K*x) + exp(-K*x))/(exp(K) + exp(-K));
  }

  double K;
};

/* Right-hand side */

class CustomRightHandSide1: public DefaultFunction
{
public:
  CustomRightHandSide1(double K, double d_u, double sigma)
    : DefaultFunction(), d_u(d_u), sigma(sigma) 
  {
    cef1 = new CustomExactFunction1();
    cef2 = new CustomExactFunction2(K);
  };

  virtual scalar value(double x, double y) const 
  {
    double Laplace_u = cef1->ddxx(x) * cef1->val(y)
                       + cef1->val(x) * cef1->ddxx(y);
    double u = cef1->val(x) * cef1->val(y);
    double v = cef2->val(x) * cef2->val(y);
    return -d_u * d_u * Laplace_u - u + sigma * v;
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(10);
  }

  ~CustomRightHandSide1() { delete cef1; delete cef2;}

  CustomExactFunction1* cef1;
  CustomExactFunction2* cef2;
  double d_u, sigma;
};

class CustomRightHandSide2: public DefaultFunction
{
public:
  CustomRightHandSide2(double K, double d_v)
    : DefaultFunction(), d_v(d_v) 
  {
    cef1 = new CustomExactFunction1();
    cef2 = new CustomExactFunction2(K);
  };

  virtual scalar value(double x, double y) const 
  {
    double Laplace_v = cef2->ddxx(x) * cef2->val(y)
                       + cef2->val(x) * cef2->ddxx(y);
    double u = cef1->val(x) * cef1->val(y);
    double v = cef2->val(x) * cef2->val(y);
    return -d_v*d_v * Laplace_v - u + v;
  }

  virtual Ord ord(Ord x, Ord y) const {
    return Ord(10);
  }

  ~CustomRightHandSide2() { delete cef1; delete cef2;}

  CustomExactFunction1* cef1;
  CustomExactFunction2* cef2;
  double d_v;
};

/* Exact solution */

class ExactSolutionFitzHughNagumo1 : public ExactSolutionScalar
{
public:
  ExactSolutionFitzHughNagumo1(Mesh* mesh)
       : ExactSolutionScalar(mesh) 
  {
    cef1 = new CustomExactFunction1();
  }

  virtual scalar value (double x, double y) const 
  {
    return cef1->val(x)*cef1->val(y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const 
  {
    dx = cef1->dx(x)*cef1->val(y);
    dy = cef1->val(x)*cef1->ddxx(y);
  }

  virtual Ord ord(Ord x, Ord y) const 
  {
    return Ord(10);
  }

  ~ExactSolutionFitzHughNagumo1() 
  {
    delete cef1;
  }

  CustomExactFunction1* cef1;
};

class ExactSolutionFitzHughNagumo2 : public ExactSolutionScalar
{
public:
  ExactSolutionFitzHughNagumo2(Mesh* mesh, double K)
       : ExactSolutionScalar(mesh) 
  {
    cef2 = new CustomExactFunction2(K);
  }
  virtual scalar value (double x, double y) const 
  {
    return cef2->val(x)*cef2->val(y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const 
  {
    dx = cef2->dx(x)*cef2->val(y);
    dy = cef2->val(x)*cef2->dx(y);
  }

  virtual Ord ord(Ord x, Ord y) const 
  {
    return Ord(10);
  }

  ~ExactSolutionFitzHughNagumo2() 
  {
    delete cef2;
  }

  CustomExactFunction2* cef2;
};

/* Weak forms */

class CustomResidual1 : public WeakForm::VectorFormVol
{
public:
  CustomResidual1(double d_u, double sigma, CustomRightHandSide1* g1)
    : WeakForm::VectorFormVol(0, HERMES_ANY), d_u(d_u), sigma(sigma), g1(g1) { };

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<scalar> *ext) const
  {
     scalar result = 0;
     for (int i = 0; i < n; i++) {
       result += wt[i] * (    sqr(d_u) * (u_ext[0]->dx[i]*v->dx[i] + u_ext[0]->dy[i]*v->dy[i]) 
                            - u_ext[0]->val[i]*v->val[i] 
                            + sigma*u_ext[1]->val[i]*v->val[i]
                            - g1->value(e->x[i], e->y[i])*v->val[i]
     	                 );
     }
 
     return result;
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
     Ord result = 0;
     for (int i = 0; i < n; i++) {
       result += wt[i] * (    sqr(d_u) * (u_ext[0]->dx[i]*v->dx[i] + u_ext[0]->dy[i]*v->dy[i]) 
                            - u_ext[0]->val[i]*v->val[i] 
                            + sigma*u_ext[1]->val[i]*v->val[i]
                            - g1->ord(e->x[i], e->y[i])*v->val[i]
			 );
     }

     return result;
  }

  // This is to make the form usable in rk_time_step().
  virtual WeakForm::VectorFormVol* clone() {
    return new CustomResidual1(*this);
  }

  private:
    double d_u;
    double sigma;
    CustomRightHandSide1* g1;
};

class CustomResidual2 : public WeakForm::VectorFormVol
{
public:
  CustomResidual2(double d_v, CustomRightHandSide2* g2)
    : WeakForm::VectorFormVol(1, HERMES_ANY), d_v(d_v), g2(g2) { };

  virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                       Geom<double> *e, ExtData<scalar> *ext) const
  {
     scalar result = 0;
     for (int i = 0; i < n; i++) {
       result += wt[i] * (    sqr(d_v) * (u_ext[1]->dx[i]*v->dx[i] + u_ext[1]->dy[i]*v->dy[i]) 
                            - u_ext[0]->val[i]*v->val[i] 
                            + u_ext[1]->val[i]*v->val[i]
                            - g2->value(e->x[i], e->y[i])*v->val[i]
     	                 );
     }
 
     return result;
  }

  virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
     Ord result = 0;
     for (int i = 0; i < n; i++) {
       result += wt[i] * (    sqr(d_v) * (u_ext[1]->dx[i]*v->dx[i] + u_ext[1]->dy[i]*v->dy[i]) 
                            - u_ext[0]->val[i]*v->val[i] 
                            + u_ext[1]->val[i]*v->val[i]
                            - g2->ord(e->x[i], e->y[i])*v->val[i]
     	                 );
     }

     return result;
  }

  // This is to make the form usable in rk_time_step().
  virtual WeakForm::VectorFormVol* clone() {
    return new CustomResidual2(*this);
  }

  private:
    double d_v;
    CustomRightHandSide2* g2;
};

class WeakFormFitzHughNagumo : public WeakForm
{
public:
  WeakFormFitzHughNagumo(CustomRightHandSide1* g1, CustomRightHandSide2* g2) : WeakForm(2) 
  {
    // Jacobian.
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, HERMES_ANY, D_u * D_u));
    add_matrix_form(new DefaultMatrixFormVol(0, 0, HERMES_ANY, -1.0));
    add_matrix_form(new DefaultMatrixFormVol(0, 1, HERMES_ANY, g1->sigma, HERMES_DEFAULT_FUNCTION, HERMES_NONSYM));
    add_matrix_form(new DefaultMatrixFormVol(1, 0, HERMES_ANY, -1.0, HERMES_DEFAULT_FUNCTION, HERMES_NONSYM));
    add_matrix_form(new DefaultJacobianDiffusion(1, 1, HERMES_ANY, D_v * D_v));
    add_matrix_form(new DefaultMatrixFormVol(1, 1, HERMES_ANY, 1.0));

    // Residual.
    add_vector_form(new CustomResidual1(D_u, g1->sigma, g1));
    add_vector_form(new CustomResidual2(D_v, g2));
  }
};
