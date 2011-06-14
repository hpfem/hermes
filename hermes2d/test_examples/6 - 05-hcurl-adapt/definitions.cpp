#include "weakform/weakform.h"
#include "integrals/hcurl.h"
#include "boundaryconditions/essential_bcs.h"
#include "weakform_library/weakforms_hcurl.h"

/* Exact solution */

#include "bessel.cpp"

static void exact_sol_val(double x, double y, scalar& e0, scalar& e1)
{
  double t1 = x*x;
  double t2 = y*y;
  double t4 = sqrt(t1+t2);
  double t5 = jv(-1.0/3.0,t4);
  double t6 = 1/t4;
  double t7 = jv(2.0/3.0,t4);
  double t11 = (t5-2.0/3.0*t6*t7)*t6;
  double t12 = atan2(y,x);
  if (t12 < 0) t12 += 2.0*M_PI;
  double t13 = 2.0/3.0*t12;
  double t14 = cos(t13);
  double t17 = sin(t13);
  double t18 = t7*t17;
  double t20 = 1/t1;
  double t23 = 1/(1.0+t2*t20);
  e0 = t11*y*t14-2.0/3.0*t18/x*t23;
  e1 = -t11*x*t14-2.0/3.0*t18*y*t20*t23;
}

static void exact_sol_der(double x, double y, scalar& e1dx, scalar& e0dy)
{
  double t1 = x*x;
  double t2 = y*y;
  double t3 = t1+t2;
  double t4 = sqrt(t3);
  double t5 = jv(2.0/3.0,t4);
  double t6 = 1/t4;
  double t7 = jv(-1.0/3.0,t4);
  double t11 = (-t5-t6*t7/3.0)*t6;
  double t14 = 1/t4/t3;
  double t15 = t14*t5;
  double t21 = t7-2.0/3.0*t6*t5;
  double t22 = 1/t3*t21;
  double t27 = atan2(y,x);
  if (t27 < 0) t27 += 2.0*M_PI;
  double t28 = 2.0/3.0*t27;
  double t29 = cos(t28);
  double t32 = t21*t14;
  double t35 = t21*t6;
  double t36 = t35*t29;
  double t39 = sin(t28);
  double t41 = 1/t1;
  double t43 = 1.0+t2*t41;
  double t44 = 1/t43;
  double t47 = 4.0/3.0*t35/x*t39*y*t44;
  double t48 = t5*t29;
  double t49 = t1*t1;
  double t52 = t43*t43;
  double t53 = 1/t52;
  double t57 = t5*t39;
  double t59 = 1/t1/x;
  e1dx =-(t11*x+2.0/3.0*t15*x-2.0/3.0*t22*x)
              *t6*x*t29+t32*t1*t29-t36-t47+4.0/9.0*t48*t2/t49*t53+4.0/3.0*t57*y*t59*t44-4.0/3.0*t57*t2*y/t49/x*t53;
  e0dy = (t11*y+2.0/3.0*t15*y-2.0/3.0*t22*y)*t6*y*t29-t32*t2*t29+t36-t47-4.0/9.0*t48*t41*t53+4.0/3.0*t57*t59*t53*y;
}

class CustomExactSolution : public ExactSolutionVector
{
public:
  CustomExactSolution(Mesh* mesh) : ExactSolutionVector(mesh) {};
  ~CustomExactSolution() {};

  virtual scalar2 value(double x, double y) const 
  {
    scalar2 ex(0.0, 0.0);
    exact_sol_val(x, y,  ex.val[0], ex.val[1]);
    return ex;
  };

  virtual void derivatives (double x, double y, scalar2& dx, scalar2& dy) const 
  {
    scalar e1dx, e0dy;
    exact_sol_der(x, y, e1dx, e0dy);
    dx.val[0] = 0;
    dx.val[1] = e1dx;
    dy.val[0] = e0dy;
    dy.val[1] = 0;
    return;
  };
  
  virtual Ord ord(Ord x, Ord y) const 
  {
    return Ord(10);
  } 
};

/* Weak forms */

class CustomWeakForm : public WeakForm
{
public:
  CustomWeakForm(double mu_r, double kappa) : WeakForm(1)
  {
    cplx ii = cplx(0.0, 1.0);

    // Jacobian.
    add_matrix_form(new WeakFormsHcurl::DefaultJacobianCurlCurl(0, 0, HERMES_ANY, new HermesFunction(1.0/mu_r)));
    add_matrix_form(new WeakFormsHcurl::DefaultMatrixFormVol(0, 0, HERMES_ANY, new HermesFunction(-sqr(kappa))));
    add_matrix_form_surf(new WeakFormsHcurl::DefaultMatrixFormSurf(0, 0, HERMES_ANY, new HermesFunction(-kappa*ii)));

    // Residual.
    add_vector_form(new WeakFormsHcurl::DefaultResidualCurlCurl(0, HERMES_ANY, new HermesFunction(1.0/mu_r)));
    add_vector_form(new WeakFormsHcurl::DefaultResidualVol(0, HERMES_ANY, new HermesFunction(-sqr(kappa))));
    add_vector_form_surf(new WeakFormsHcurl::DefaultResidualSurf(0, HERMES_ANY, new HermesFunction(-kappa*ii)));
    add_vector_form_surf(new CustomVectorFormSurf());
  };

  class CustomVectorFormSurf : public WeakForm::VectorFormSurf
  {
  public:
    CustomVectorFormSurf()
              : WeakForm::VectorFormSurf(0) 
    {
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
    {
      scalar result = 0;
      for (int i = 0; i < n; i++) {
        double r = sqrt(e->x[i] * e->x[i] + e->y[i] * e->y[i]);
        double theta = atan2(e->y[i], e->x[i]);
        if (theta < 0) theta += 2.0*M_PI;
        double j13    = jv(-1.0/3.0, r),    j23    = jv(+2.0/3.0, r);
        double cost   = cos(theta),         sint   = sin(theta);
        double cos23t = cos(2.0/3.0*theta), sin23t = sin(2.0/3.0*theta);

        double Etau = e->tx[i] * (cos23t*sint*j13 - 2.0/(3.0*r)*j23*(cos23t*sint + sin23t*cost)) +
                      e->ty[i] * (-cos23t*cost*j13 + 2.0/(3.0*r)*j23*(cos23t*cost - sin23t*sint));

        result += wt[i] * cplx(cos23t*j23, -Etau) * ((v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
      }
      return -result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                    Geom<Ord> *e, ExtData<Ord> *ext) const 
    {
      return Ord(10);
    }
  };
};
