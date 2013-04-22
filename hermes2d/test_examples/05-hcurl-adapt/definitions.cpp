#include "hermes2d.h"

/* Exact solution */

#include "bessel.cpp"

static void exact_sol_val(double x, double y, std::complex<double>& e0, std::complex<double>& e1)
{
  double t1 = x*x;
  double t2 = y*y;
  double t4 = std::sqrt(t1+t2);
  double t5 = jv(-1.0/3.0,t4);
  double t6 = 1/t4;
  double t7 = jv(2.0/3.0,t4);
  double t11 = (t5-2.0/3.0*t6*t7)*t6;
  double t12 = std::atan2(y,x);
  if (t12 < 0) t12 += 2.0*M_PI;
  double t13 = 2.0/3.0*t12;
  double t14 = std::cos(t13);
  double t17 = std::sin(t13);
  double t18 = t7*t17;
  double t20 = 1/t1;
  double t23 = 1/(1.0+t2*t20);
  e0 = t11*y*t14-2.0/3.0*t18/x*t23;
  e1 = -t11*x*t14-2.0/3.0*t18*y*t20*t23;
}

static void exact_sol_der(double x, double y, std::complex<double>& e1dx, std::complex<double>& e0dy)
{
  double t1 = x*x;
  double t2 = y*y;
  double t3 = t1+t2;
  double t4 = std::sqrt(t3);
  double t5 = jv(2.0/3.0,t4);
  double t6 = 1/t4;
  double t7 = jv(-1.0/3.0,t4);
  double t11 = (-t5-t6*t7/3.0)*t6;
  double t14 = 1/t4/t3;
  double t15 = t14*t5;
  double t21 = t7-2.0/3.0*t6*t5;
  double t22 = 1/t3*t21;
  double t27 = std::atan2(y,x);
  if (t27 < 0) t27 += 2.0*M_PI;
  double t28 = 2.0/3.0*t27;
  double t29 = std::cos(t28);
  double t32 = t21*t14;
  double t35 = t21*t6;
  double t36 = t35*t29;
  double t39 = std::sin(t28);
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

class CustomExactSolution : public Hermes::Hermes2D::ExactSolutionVector<std::complex<double> >
{
public:
  CustomExactSolution(MeshSharedPtr mesh) : Hermes::Hermes2D::ExactSolutionVector<std::complex<double> >(mesh) {};
  ~CustomExactSolution() {};

  virtual Scalar2<std::complex<double> > value(double x, double y) const 
  {
    Scalar2<std::complex<double> >ex(0.0, 0.0);
#pragma omp critical (custom)
    {
      exact_sol_val(x, y,  ex[0], ex[1]);
    }
    return ex;
  };

  virtual void derivatives (double x, double y, Scalar2<std::complex<double> >& dx, Scalar2<std::complex<double> >& dy) const 
  {
    std::complex<double> e1dx, e0dy;
#pragma omp critical (custom)
    {
      exact_sol_der(x, y, e1dx, e0dy);
    }
    dx[0] = 0;
    dx[1] = e1dx;
    dy[0] = e0dy;
    dy[1] = 0;
    return;
  };

  virtual Hermes::Ord ord(Hermes::Ord x, Hermes::Ord y) const 
  {
    return Hermes::Ord(10);
  }
  
  virtual MeshFunction<std::complex<double> >* clone() const
  {
    return new CustomExactSolution(this->mesh);
  }
};

/* Weak forms */

class CustomWeakForm : public Hermes::Hermes2D::WeakForm<std::complex<double> >
{
public:
  CustomWeakForm(double mu_r, double kappa) : Hermes::Hermes2D::WeakForm<std::complex<double> >(1)
  {
    std::complex<double> ii = std::complex<double>(0.0, 1.0);

    // Jacobian.
    add_matrix_form(new WeakFormsHcurl::DefaultJacobianCurlCurl<std::complex<double> >(0, 0, HERMES_ANY, 1.0/mu_r));
    add_matrix_form(new WeakFormsHcurl::DefaultMatrixFormVol<std::complex<double> >(0, 0, HERMES_ANY, -sqr(kappa)));
    add_matrix_form_surf(new WeakFormsHcurl::DefaultMatrixFormSurf<std::complex<double> >(0, 0, HERMES_ANY, -kappa*ii));

    // Residual.
    add_vector_form(new WeakFormsHcurl::DefaultResidualCurlCurl<std::complex<double> >(0, HERMES_ANY, 1.0/mu_r));
    add_vector_form(new WeakFormsHcurl::DefaultResidualVol<std::complex<double> >(0, HERMES_ANY, -sqr(kappa)));
    add_vector_form_surf(new WeakFormsHcurl::DefaultResidualSurf<std::complex<double> >(0, HERMES_ANY, -kappa*ii));
    add_vector_form_surf(new CustomVectorFormSurf());
  };

  class CustomVectorFormSurf : public VectorFormSurf<std::complex<double> >
  {
  public:
    CustomVectorFormSurf()
      : VectorFormSurf<std::complex<double> >(0) 
    {
    }

    virtual std::complex<double> value(int n, double *wt, Func<std::complex<double> > *u_ext[], 
      Func<double> *v, Geom<double> *e, Func<std::complex<double> > **ext) const 
    {
      std::complex<double> result = 0;
#pragma omp critical (jv)
      {
        for (int i = 0; i < n; i++) {
          double r = std::sqrt(e->x[i] * e->x[i] + e->y[i] * e->y[i]);
          double theta = std::atan2(e->y[i], e->x[i]);
          if (theta < 0) theta += 2.0*M_PI;
          double j13    = jv(-1.0/3.0, r),    j23    = jv(+2.0/3.0, r);
          double cost   = std::cos(theta),         sint   = std::sin(theta);
          double cos23t = std::cos(2.0/3.0*theta), sin23t = std::sin(2.0/3.0*theta);

          double Etau = e->tx[i] * (cos23t*sint*j13 - 2.0/(3.0*r)*j23*(cos23t*sint + sin23t*cost)) +
            e->ty[i] * (-cos23t*cost*j13 + 2.0/(3.0*r)*j23*(cos23t*cost - sin23t*sint));

          result += wt[i] * std::complex<double>(cos23t*j23, -Etau) * ((v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
        }
      }
      return -result;
    }

    virtual Hermes::Ord ord(int n, double *wt, Func<Hermes::Ord> *u_ext[], Func<Hermes::Ord> *v,
      Geom<Hermes::Ord> *e, Func<Ord> **ext) const 
    {
      return Hermes::Ord(10);
    }

    virtual VectorFormSurf<std::complex<double> >* clone() const { return new CustomVectorFormSurf(); }
  };
};
