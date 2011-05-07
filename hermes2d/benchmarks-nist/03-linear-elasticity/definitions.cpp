#include "hermes2d.h"

/* Exact solution */

class CustomExactSolutionU : public ExactSolutionScalar
{
public:
  CustomExactSolutionU(Mesh* mesh, double E, double nu, double lambda, double Q) 
            : ExactSolutionScalar(mesh), E(E), nu(nu), lambda(lambda), Q(Q) {
    k = 3.0 - 4.0 * nu;
    mu = E / (2.0 * (1.0 + nu));
    A = -E * (1 - nu * nu)/(1 - 2 * nu);
    B = -E * (1 - nu * nu)/(2 - 2 * nu);
    C = -E * (1 - nu * nu)/((1 - 2 * nu) * (2 - 2 * nu));
    D = 1.0 / (2 * mu);

    u_F = (k - Q * (lambda + 1));
    v_F = (k + Q * (lambda + 1));
  };
  
  double get_angle(double y, double x) const {
    double theta = atan2(y, x);
    if (theta < 0)
      theta += 2 * M_PI;
    return theta;
  }

  double d_theta_dx(double x, double y) const {
    return -y/(x*x + y*y);
  }

  double d_theta_dxd_theta_dx(double x, double y) const {
    return 2*x*y/((x*x + y*y)*(x*x + y*y));
  }

  double d_theta_dy(double x, double y) const {
    return x/(x*x + y*y) ;
  }

  double d_theta_dyd_theta_dy(double x, double y) const {
    return -2*x*y/((x*x + y*y)*(x*x + y*y));
  }

  double d_theta_dxd_theta_dy(double x, double y) const {
    return (y*y - x*x)/((x*x + y*y)*(x*x + y*y));
  }

  double r(double x, double y) const {
    return pow((x*x + y*y), (lambda/2.0));  // r^labbda
  }

  double drdx(double x, double y) const {
    return lambda * x * pow((x*x + y*y), (lambda/2.0 - 1.0));
  }

  double drdxdrdx(double x, double y) const {
    return lambda * (pow((x*x + y*y), (lambda/2.0 - 1.0)) + (lambda - 2.0) 
           * x * x * pow((x*x + y*y), (lambda/2.0 - 2.0)));
  }

  double drdy(double x, double y) const {
    return lambda * y * pow((x*x + y*y), (lambda/2.0 - 1.0));
  }

  double drdydrdy(double x, double y) const {
    return lambda * (pow((x*x + y*y), (lambda/2.0 - 1.0)) + (lambda - 2.0) 
           * y * y * pow((x*x + y*y), (lambda/2.0 - 2.0)));
  }

  double drdxdrdy(double x, double y) const {
    return lambda * 2.0 * x * y * (lambda/2.0 - 1) * pow((x*x + y*y), 
           (lambda/2.0 - 2.0));
  }

  double u_r(double x, double y) const {
    return (u_F * cos(lambda * get_angle(y, x)) - lambda * cos((lambda - 2) 
           * get_angle(y, x)));
  }

  double du_rdx(double x, double y) const {
    return (u_F * (-1) * lambda * sin(lambda * get_angle(y, x)) * d_theta_dx(x, y)) 
           - (lambda * (-1) * (lambda - 2) * sin((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));
  }

  double du_rdxdu_rdx(double x, double y) const {
    return (u_F * (-1) * lambda * (sin(lambda * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
           + lambda * d_theta_dx(x, y) * d_theta_dx(x, y) * cos(lambda * get_angle(y, x)))) 
	   - (lambda * (-1) * (lambda - 2) * (sin((lambda - 2) * get_angle(y, x)) 
           * d_theta_dxd_theta_dx(x, y) + (lambda - 2) * d_theta_dx(x, y) * d_theta_dx(x, y) 
           * cos((lambda-2) * get_angle(y, x)) ));
  }

  double du_rdy(double x, double y) const {
    return (u_F * (-1) * lambda * sin(lambda * get_angle(y, x)) * d_theta_dy(x, y)) 
           - (lambda * (-1) * (lambda - 2) * sin((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
  }

  double du_rdydu_rdy(double x, double y) const {
    return (u_F * (-1) * lambda * (sin(lambda * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
           + lambda * d_theta_dy(x, y) * d_theta_dy(x, y) * cos(lambda * get_angle(y, x)))) 
	   - (lambda * (-1) * (lambda - 2) * (sin((lambda - 2) * get_angle(y, x)) 
           * d_theta_dyd_theta_dy(x, y) + (lambda - 2) * d_theta_dy(x, y) * d_theta_dy(x, y) 
           * cos((lambda-2) * get_angle(y, x)) ));
  }

  double du_rdxdu_rdy(double x, double y) const {
    return (u_F * (-1) * lambda * (sin(lambda * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
           + lambda * d_theta_dx(x, y) * d_theta_dy(x, y) * cos(lambda * get_angle(y, x)))) 
	   - (lambda * (-1) * (lambda - 2) * (sin((lambda - 2) * get_angle(y, x)) 
           * d_theta_dxd_theta_dy(x, y) + (lambda - 2) * d_theta_dx(x, y) * d_theta_dy(x, y) 
           * cos((lambda-2) * get_angle(y, x)) ));
  }

  double v_r(double x, double y) const {
    return (v_F * sin(lambda * get_angle(y, x)) + lambda * sin((lambda - 2) * get_angle(y, x)));
  }

  double dv_rdx(double x, double y) const {
    return (v_F * lambda * cos(lambda * get_angle(y, x)) * d_theta_dx(x, y)) 
           + (lambda * (lambda - 2) * cos((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));
  }

  double dv_rdxdv_rdx(double x, double y) const {
    return (v_F * lambda * (cos(lambda * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
           + lambda * d_theta_dx(x, y) * d_theta_dx(x, y) * (-1) * sin(lambda * get_angle(y, x)))) 
	   + (lambda * (lambda - 2) * (cos((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
           + (lambda - 2) * d_theta_dx(x, y) * d_theta_dx(x, y) * (-1) * sin((lambda-2) * get_angle(y, x)) ));
  }

  double dv_rdy(double x, double y) const {
    return (v_F * lambda * cos(lambda * get_angle(y, x)) * d_theta_dy(x, y)) + (lambda * (lambda - 2) 
           * cos((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
  }

  double dv_rdydv_rdy(double x, double y) const {
    return (v_F * lambda * (cos(lambda * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
           + lambda * d_theta_dy(x, y) * d_theta_dy(x, y) * (-1) * sin(lambda * get_angle(y, x)))) 
	   + (lambda * (lambda - 2) * (cos((lambda - 2) * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
           + (lambda - 2) * d_theta_dy(x, y) * d_theta_dy(x, y) * (-1) * sin((lambda-2) * get_angle(y, x)) ));
  }

  double dv_rdxdv_rdy(double x, double y) const {
    return (v_F * lambda * (cos(lambda * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
           + lambda * (-1) * d_theta_dx(x, y) * d_theta_dy(x, y) * sin(lambda * get_angle(y, x)))) 
	   + (lambda * (lambda - 2) * (cos((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
           + (lambda - 2) * (-1) * d_theta_dx(x, y) * d_theta_dy(x, y) * sin((lambda-2) * get_angle(y, x)) ));
  }

  // \partial^2 u / \partial x^2 
  double dudxdudx(double x, double y) const {
    return D * (drdxdrdx(x, y) * u_r(x, y) + 2 * drdx(x, y) * du_rdx(x, y) + r(x, y) * du_rdxdu_rdx(x, y) );
  }

  double dudydudy(double x, double y) const {
    return D * (drdydrdy(x, y) * u_r(x, y) + 2 * drdy(x, y) * du_rdy(x, y) + r(x, y) * du_rdydu_rdy(x, y) );
  }

  double dudxdudy(double x, double y) const {
    return D * (drdxdrdy(x, y) * u_r(x, y) + drdx(x, y) * du_rdy(x, y) + drdy(x, y) * du_rdx(x, y) + r(x, y) 
           * du_rdxdu_rdy(x, y) );
  }

  // \partial^2 v / \partial x^2
  double dvdxdvdx(double x, double y) const {
    return D * (drdxdrdx(x, y) * v_r(x, y) + 2 * drdx(x, y) * dv_rdx(x, y) + r(x, y) * dv_rdxdv_rdx(x, y) );
  }

  double dvdydvdy(double x, double y) const {
    return D * (drdydrdy(x, y) * v_r(x, y) + 2 * drdy(x, y) * dv_rdy(x, y) + r(x, y) * dv_rdydv_rdy(x, y) );
  }

  double dvdxdvdy(double x, double y) const {
    return D * (drdxdrdy(x, y) * v_r(x, y) + drdx(x, y) * dv_rdy(x, y) + drdy(x, y) * dv_rdx(x, y) + r(x, y) 
           * dv_rdxdv_rdy(x, y) );
  }

  // Exact solution u(x,y) and its derivatives.
  virtual double value (double x, double y) const {
    return D * r(x, y) * u_r(x, y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = D * drdx(x, y) * (u_F * cos(lambda * get_angle(y, x)) - lambda * cos((lambda - 2) * get_angle(y, x))) +
         D * r(x, y) * (u_F * (-1) * lambda * sin(lambda * get_angle(y, x)) * d_theta_dx(x, y)) - 
         D * r(x, y) * (lambda * (-1) * (lambda - 2) * sin((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));

    dy = D * drdy(x, y) * (u_F * cos(lambda * get_angle(y, x)) - lambda * cos((lambda - 2) * get_angle(y, x))) +
         D * r(x, y) * (u_F * (-1) * lambda * sin(lambda * get_angle(y, x)) * d_theta_dy(x, y)) - 
         D * r(x, y) * (lambda * (-1) * (lambda - 2) * sin((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
  }

  virtual Ord ord (Ord x, Ord y) const {
    return Ord(10);
  }

  double E, nu, lambda, Q, k, mu;
  double A, B, C, D, u_F, v_F;
};

class CustomExactSolutionV : public ExactSolutionScalar
{
public:
  CustomExactSolutionV(Mesh* mesh, double E, double nu, double lambda, double Q) 
            : ExactSolutionScalar(mesh), E(E), nu(nu), lambda(lambda), Q(Q) {
    k = 3.0 - 4.0 * nu;
    mu = E / (2.0 * (1.0 + nu));
    A = -E * (1 - nu * nu)/(1 - 2 * nu);
    B = -E * (1 - nu * nu)/(2 - 2 * nu);
    C = -E * (1 - nu * nu)/((1 - 2 * nu) * (2 - 2 * nu));
    D = 1.0 / (2 * mu);

    u_F = (k - Q * (lambda + 1));
    v_F = (k + Q * (lambda + 1));
  };
  
  double get_angle(double y, double x) const {
    double theta = atan2(y, x);
    if (theta < 0)
      theta += 2 * M_PI;
    return theta;
  }

  double d_theta_dx(double x, double y) const {
    return -y/(x*x + y*y);
  }

  double d_theta_dxd_theta_dx(double x, double y) const {
    return 2*x*y/((x*x + y*y)*(x*x + y*y));
  }

  double d_theta_dy(double x, double y) const {
    return x/(x*x + y*y) ;
  }

  double d_theta_dyd_theta_dy(double x, double y) const {
    return -2*x*y/((x*x + y*y)*(x*x + y*y));
  }

  double d_theta_dxd_theta_dy(double x, double y) const {
    return (y*y - x*x)/((x*x + y*y)*(x*x + y*y));
  }

  double r(double x, double y) const {
    return pow((x*x + y*y), (lambda/2.0));  // r^labbda
  }

  double drdx(double x, double y) const {
    return lambda * x * pow((x*x + y*y), (lambda/2.0 - 1.0));
  }

  double drdxdrdx(double x, double y) const {
    return lambda * (pow((x*x + y*y), (lambda/2.0 - 1.0)) + (lambda - 2.0) 
           * x * x * pow((x*x + y*y), (lambda/2.0 - 2.0)));
  }

  double drdy(double x, double y) const {
    return lambda * y * pow((x*x + y*y), (lambda/2.0 - 1.0));
  }

  double drdydrdy(double x, double y) const {
    return lambda * (pow((x*x + y*y), (lambda/2.0 - 1.0)) + (lambda - 2.0) 
           * y * y * pow((x*x + y*y), (lambda/2.0 - 2.0)));
  }

  double drdxdrdy(double x, double y) const {
    return lambda * 2.0 * x * y * (lambda/2.0 - 1) * pow((x*x + y*y), (lambda/2.0 - 2.0));
  }

  double u_r(double x, double y) const {
    return (u_F * cos(lambda * get_angle(y, x)) - lambda * cos((lambda - 2) * get_angle(y, x)));
  }

  double du_rdx(double x, double y) const {
    return (u_F * (-1) * lambda * sin(lambda * get_angle(y, x)) * d_theta_dx(x, y)) 
           - (lambda * (-1) * (lambda - 2) * sin((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));
  }

  double du_rdxdu_rdx(double x, double y) const {
    return (u_F * (-1) * lambda * (sin(lambda * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
           + lambda * d_theta_dx(x, y) * d_theta_dx(x, y) * cos(lambda * get_angle(y, x)))) 
	   - (lambda * (-1) * (lambda - 2) * (sin((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
           + (lambda - 2) * d_theta_dx(x, y) * d_theta_dx(x, y) * cos((lambda-2) * get_angle(y, x)) ));
  }

  double du_rdy(double x, double y) const {
    return (u_F * (-1) * lambda * sin(lambda * get_angle(y, x)) * d_theta_dy(x, y)) - (lambda * (-1) * 
           (lambda - 2) * sin((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
  }

  double du_rdydu_rdy(double x, double y) const {
    return (u_F * (-1) * lambda * (sin(lambda * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
           + lambda * d_theta_dy(x, y) * d_theta_dy(x, y) * cos(lambda * get_angle(y, x)))) 
	   - (lambda * (-1) * (lambda - 2) * (sin((lambda - 2) * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
           + (lambda - 2) * d_theta_dy(x, y) * d_theta_dy(x, y) * cos((lambda-2) * get_angle(y, x)) ));
  }

  double du_rdxdu_rdy(double x, double y) const {
    return (u_F * (-1) * lambda * (sin(lambda * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
           + lambda * d_theta_dx(x, y) * d_theta_dy(x, y) * cos(lambda * get_angle(y, x)))) 
	   - (lambda * (-1) * (lambda - 2) * (sin((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
           + (lambda - 2) * d_theta_dx(x, y) * d_theta_dy(x, y) * cos((lambda-2) * get_angle(y, x)) ));
  }

  double v_r(double x, double y) const {
    return (v_F * sin(lambda * get_angle(y, x)) + lambda * sin((lambda - 2) * get_angle(y, x)));
  }

  double dv_rdx(double x, double y) const {
    return (v_F * lambda * cos(lambda * get_angle(y, x)) * d_theta_dx(x, y)) + (lambda * (lambda - 2) 
           * cos((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));
  }

  double dv_rdxdv_rdx(double x, double y) const {
    return (v_F * lambda * (cos(lambda * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
           + lambda * d_theta_dx(x, y) * d_theta_dx(x, y) * (-1) * sin(lambda * get_angle(y, x)))) 
	   + (lambda * (lambda - 2) * (cos((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) 
           + (lambda - 2) * d_theta_dx(x, y) * d_theta_dx(x, y) * (-1) * sin((lambda-2) * get_angle(y, x)) ));
  }

  double dv_rdy(double x, double y) const {
    return (v_F * lambda * cos(lambda * get_angle(y, x)) * d_theta_dy(x, y)) + (lambda * (lambda - 2) 
           * cos((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
  }

  double dv_rdydv_rdy(double x, double y) const {
    return (v_F * lambda * (cos(lambda * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
           + lambda * d_theta_dy(x, y) * d_theta_dy(x, y) * (-1) * sin(lambda * get_angle(y, x)))) 
	   + (lambda * (lambda - 2) * (cos((lambda - 2) * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) 
           + (lambda - 2) * d_theta_dy(x, y) * d_theta_dy(x, y) * (-1) * sin((lambda-2) * get_angle(y, x)) ));
  }

  double dv_rdxdv_rdy(double x, double y) const {
    return (v_F * lambda * (cos(lambda * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
           + lambda * (-1) * d_theta_dx(x, y) * d_theta_dy(x, y) * sin(lambda * get_angle(y, x)))) 
	   + (lambda * (lambda - 2) * (cos((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) 
           + (lambda - 2) * (-1) * d_theta_dx(x, y) * d_theta_dy(x, y) * sin((lambda-2) * get_angle(y, x)) ));
  }

  // \partial^2 u / \partial x^2 
  double dudxdudx(double x, double y) const {
    return D * (drdxdrdx(x, y) * u_r(x, y) + 2 * drdx(x, y) * du_rdx(x, y) + r(x, y) * du_rdxdu_rdx(x, y) );
  }

  double dudydudy(double x, double y) const {
    return D * (drdydrdy(x, y) * u_r(x, y) + 2 * drdy(x, y) * du_rdy(x, y) + r(x, y) * du_rdydu_rdy(x, y) );
  }

  double dudxdudy(double x, double y) const {
    return D * (drdxdrdy(x, y) * u_r(x, y) + drdx(x, y) * du_rdy(x, y) + drdy(x, y) * du_rdx(x, y) + r(x, y) 
           * du_rdxdu_rdy(x, y) );
  }

  // \partial^2 v / \partial x^2
  double dvdxdvdx(double x, double y) const {
    return D * (drdxdrdx(x, y) * v_r(x, y) + 2 * drdx(x, y) * dv_rdx(x, y) + r(x, y) * dv_rdxdv_rdx(x, y) );
  }

  double dvdydvdy(double x, double y) const {
    return D * (drdydrdy(x, y) * v_r(x, y) + 2 * drdy(x, y) * dv_rdy(x, y) + r(x, y) * dv_rdydv_rdy(x, y) );
  }

  double dvdxdvdy(double x, double y) const {
    return D * (drdxdrdy(x, y) * v_r(x, y) + drdx(x, y) * dv_rdy(x, y) + drdy(x, y) * dv_rdx(x, y) + r(x, y) 
           * dv_rdxdv_rdy(x, y) );
  }

  // Exact solution v(x,y) and its derivatives.
  virtual double value (double x, double y) const {
    return D * r(x, y) * v_r(x, y);
  }

  virtual void derivatives (double x, double y, scalar& dx, scalar& dy) const {
    dx = D * drdx(x, y) * (v_F * sin(lambda * get_angle(y, x)) + lambda * sin((lambda - 2) * get_angle(y, x))) +
         D * r(x, y) * (v_F * lambda * cos(lambda * get_angle(y, x)) * d_theta_dx(x, y)) + 
         D * r(x, y) * (lambda * (lambda - 2) * cos((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));

    dy = D * drdy(x, y) * (v_F * sin(lambda * get_angle(y, x)) + lambda * sin((lambda - 2) * get_angle(y, x))) +
         D * r(x, y) * (v_F * lambda * cos(lambda * get_angle(y, x)) * d_theta_dy(x, y)) + 
         D * r(x, y) * (lambda * (lambda - 2) * cos((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
  }

  virtual Ord ord (Ord x, Ord y) const {
    return Ord(10);
  }

  double E, nu, lambda, Q, k, mu;
  double A, B, C, D, u_F, v_F;
};

/* Weak forms */

// We cannot use default Hermes linear elasticity forms since the equations 
// in Mitchell's paper are different (E, nu, lambda and mu do not satisfy 
// standard relations http://en.wikipedia.org/wiki/Linear_elasticity).
class CustomWeakFormElasticityNIST : public WeakForm
{
public:
  CustomWeakFormElasticityNIST(double E, double nu, double mu, double lambda) : WeakForm(2)
  {
    // Jacobian.
    add_matrix_form(new CustomMatrixFormVolElasticityNIST_0_0(E, nu));
    add_matrix_form(new CustomMatrixFormVolElasticityNIST_0_1(E, nu));
    add_matrix_form(new CustomMatrixFormVolElasticityNIST_1_1(E, nu));
    // Residual.
    add_vector_form(new CustomVectorFormVolElasticityNIST_0(E, nu));
    add_vector_form(new CustomVectorFormVolElasticityNIST_1(E, nu));
  }

private:
  class CustomMatrixFormVolElasticityNIST_0_0 : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVolElasticityNIST_0_0(double E, double nu) 
      : WeakForm::MatrixFormVol(0, 0, HERMES_ANY, HERMES_SYM) {
      A = -E * (1 - nu * nu)/(1 - 2 * nu);
      B = -E * (1 - nu * nu)/(2 - 2 * nu);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                 Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
    {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (A * u->dx[i] * v->dx[i] + B * u->dy[i] * v->dy[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
            Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      Ord result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (A * u->dx[i] * v->dx[i] + B * u->dy[i] * v->dy[i]);
      return result;
    }

    double A, B;
  };

  class CustomMatrixFormVolElasticityNIST_0_1 : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVolElasticityNIST_0_1(double E, double nu) 
      : WeakForm::MatrixFormVol(0, 1, HERMES_ANY, HERMES_SYM) { 
      C = -E * (1 - nu * nu)/((1 - 2 * nu) * (2 - 2 * nu));
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
    {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (C * u->dx[i] * v->dy[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      Ord result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (C * u->dx[i] * v->dy[i]);
      return result;
    }

    double C;
  };

  class CustomMatrixFormVolElasticityNIST_1_1 : public WeakForm::MatrixFormVol
  {
  public:
    CustomMatrixFormVolElasticityNIST_1_1(double E, double nu) 
      : WeakForm::MatrixFormVol(1, 1, HERMES_ANY, HERMES_SYM) { 
      A = -E * (1 - nu * nu)/(1 - 2 * nu);
      B = -E * (1 - nu * nu)/(2 - 2 * nu);
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
    {
      scalar result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (B * u->dx[i] * v->dx[i] + A * u->dy[i] * v->dy[i]);
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, 
                    Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      Ord result = 0;
      for (int i = 0; i < n; i++)
        result += wt[i] * (B * u->dx[i] * v->dx[i] + A * u->dy[i] * v->dy[i]);
      return result;
    }
  
    double A, B;
  };

  class CustomVectorFormVolElasticityNIST_0 : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVolElasticityNIST_0(double E, double nu) 
      : WeakForm::VectorFormVol(0, HERMES_ANY) {
      A = -E * (1 - nu * nu)/(1 - 2 * nu);
      B = -E * (1 - nu * nu)/(2 - 2 * nu);
      C = -E * (1 - nu * nu)/((1 - 2 * nu) * (2 - 2 * nu));
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[], 
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
    {
      scalar result = 0;
      for (int i = 0; i < n; i++) {
        // Contribution of matrix form 0, 0.
        result += wt[i] * (A * u_ext[0]->dx[i] * v->dx[i] + B * u_ext[0]->dy[i] * v->dy[i]);
        // Contribution of matrix form 0, 1.
        result += wt[i] * (C * u_ext[1]->dx[i] * v->dy[i]);
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      Ord result = 0;
      for (int i = 0; i < n; i++) {
        // Contribution of matrix form 0, 0.
        result += wt[i] * (A * u_ext[0]->dx[i] * v->dx[i] + B * u_ext[0]->dy[i] * v->dy[i]);
        // Contribution of matrix form 0, 1.
        result += wt[i] * (C * u_ext[1]->dx[i] * v->dy[i]);
      }
      return result;
    }

    double A, B, C;
  };

  class CustomVectorFormVolElasticityNIST_1 : public WeakForm::VectorFormVol
  {
  public:
    CustomVectorFormVolElasticityNIST_1(double E, double nu) 
      : WeakForm::VectorFormVol(0, HERMES_ANY) {
      A = -E * (1 - nu * nu)/(1 - 2 * nu);
      B = -E * (1 - nu * nu)/(2 - 2 * nu);
      C = -E * (1 - nu * nu)/((1 - 2 * nu) * (2 - 2 * nu));
    }

    virtual scalar value(int n, double *wt, Func<scalar> *u_ext[],  
                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
    {
      scalar result = 0;
      for (int i = 0; i < n; i++) {
        // Contribution of matrix form 1, 0.
        result += wt[i] * (C * u_ext[0]->dy[i] * v->dx[i]);
        // Contribution of matrix form 1, 1.
        result += wt[i] * (B * u_ext[1]->dx[i] * v->dx[i] + A * u_ext[1]->dy[i] * v->dy[i]);
      }
      return result;
    }

    virtual Ord ord(int n, double *wt, Func<Ord> *u_ext[], 
                    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
    {
      Ord result = 0;
      for (int i = 0; i < n; i++) {
        // Contribution of matrix form 1, 0.
        result += wt[i] * (C * u_ext[0]->dy[i] * v->dx[i]);
        // Contribution of matrix form 1, 1.
        result += wt[i] * (B * u_ext[1]->dx[i] * v->dx[i] + A * u_ext[1]->dy[i] * v->dy[i]);
      }
      return result;
    }

    double A, B, C;
  };

};

