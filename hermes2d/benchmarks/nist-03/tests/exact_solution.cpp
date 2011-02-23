const double A = -E * (1 - nu * nu)/(1 - 2 * nu);
const double B = -E * (1 - nu * nu)/(2 - 2 * nu);
const double C = -E * (1 - nu * nu)/((1 - 2 * nu) * (2 - 2 * nu));
const double D = 1.0 / (2 * G);

const double u_F = (k - Q * (lambda + 1));
const double v_F = (k + Q * (lambda + 1));

double get_angle(double y, double x)
{
  double theta = atan2(y, x);
  if (theta < 0)
    theta += 2 * M_PI;
  return theta;
}

double d_theta_dx(double x, double y)
{
  return -y/(x*x + y*y);
}

double d_theta_dxd_theta_dx(double x, double y)
{
  return 2*x*y/((x*x + y*y)*(x*x + y*y));
}

double d_theta_dy(double x, double y)
{
  return x/(x*x + y*y) ;
}

double d_theta_dyd_theta_dy(double x, double y)
{
  return -2*x*y/((x*x + y*y)*(x*x + y*y));
}

double d_theta_dxd_theta_dy(double x, double y)
{
  return (y*y - x*x)/((x*x + y*y)*(x*x + y*y));
}

double r(double x, double y)
{
  return pow((x*x + y*y), (lambda/2.0));  // r^labbda
}

double drdx(double x, double y)
{
  return lambda * x * pow((x*x + y*y), (lambda/2.0 - 1.0));
}

double drdxdrdx(double x, double y)
{
  return lambda * (pow((x*x + y*y), (lambda/2.0 - 1.0)) + (lambda - 2.0) * x * x * pow((x*x + y*y), (lambda/2.0 - 2.0)));
}

double drdy(double x, double y)
{
  return lambda * y * pow((x*x + y*y), (lambda/2.0 - 1.0));
}

double drdydrdy(double x, double y)
{
  return lambda * (pow((x*x + y*y), (lambda/2.0 - 1.0)) + (lambda - 2.0) * y * y * pow((x*x + y*y), (lambda/2.0 - 2.0)));
}

double drdxdrdy(double x, double y)
{
  return lambda * 2.0 * x * y * (lambda/2.0 - 1) * pow((x*x + y*y), (lambda/2.0 - 2.0));
}

double u_r(double x, double y)
{
  return (u_F * cos(lambda * get_angle(y, x)) - lambda * cos((lambda - 2) * get_angle(y, x)));
}

double du_rdx(double x, double y)
{
  return (u_F * (-1) * lambda * sin(lambda * get_angle(y, x)) * d_theta_dx(x, y)) - (lambda * (-1) * (lambda - 2) * sin((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));
}

double du_rdxdu_rdx(double x, double y)
{
  return (u_F * (-1) * lambda * (sin(lambda * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) + lambda * d_theta_dx(x, y) * d_theta_dx(x, y) * cos(lambda * get_angle(y, x)))) 
	     - (lambda * (-1) * (lambda - 2) * (sin((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) + (lambda - 2) * d_theta_dx(x, y) * d_theta_dx(x, y) * cos((lambda-2) * get_angle(y, x)) ));
}

double du_rdy(double x, double y)
{
  return (u_F * (-1) * lambda * sin(lambda * get_angle(y, x)) * d_theta_dy(x, y)) - (lambda * (-1) * (lambda - 2) * sin((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
}

double du_rdydu_rdy(double x, double y)
{
  return (u_F * (-1) * lambda * (sin(lambda * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) + lambda * d_theta_dy(x, y) * d_theta_dy(x, y) * cos(lambda * get_angle(y, x)))) 
	     - (lambda * (-1) * (lambda - 2) * (sin((lambda - 2) * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) + (lambda - 2) * d_theta_dy(x, y) * d_theta_dy(x, y) * cos((lambda-2) * get_angle(y, x)) ));
}

double du_rdxdu_rdy(double x, double y)
{
  return (u_F * (-1) * lambda * (sin(lambda * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) + lambda * d_theta_dx(x, y) * d_theta_dy(x, y) * cos(lambda * get_angle(y, x)))) 
	     - (lambda * (-1) * (lambda - 2) * (sin((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) + (lambda - 2) * d_theta_dx(x, y) * d_theta_dy(x, y) * cos((lambda-2) * get_angle(y, x)) ));
}

double v_r(double x, double y)
{
  return (v_F * sin(lambda * get_angle(y, x)) + lambda * sin((lambda - 2) * get_angle(y, x)));
}

double dv_rdx(double x, double y)
{
  return (v_F * lambda * cos(lambda * get_angle(y, x)) * d_theta_dx(x, y)) + (lambda * (lambda - 2) * cos((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));
}

double dv_rdxdv_rdx(double x, double y)
{
  return (v_F * lambda * (cos(lambda * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) + lambda * d_theta_dx(x, y) * d_theta_dx(x, y) * (-1) * sin(lambda * get_angle(y, x)))) 
	     + (lambda * (lambda - 2) * (cos((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dx(x, y) + (lambda - 2) * d_theta_dx(x, y) * d_theta_dx(x, y) * (-1) * sin((lambda-2) * get_angle(y, x)) ));
}

double dv_rdy(double x, double y)
{
  return (v_F * lambda * cos(lambda * get_angle(y, x)) * d_theta_dy(x, y)) + (lambda * (lambda - 2) * cos((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));
}

double dv_rdydv_rdy(double x, double y)
{
  return (v_F * lambda * (cos(lambda * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) + lambda * d_theta_dy(x, y) * d_theta_dy(x, y) * (-1) * sin(lambda * get_angle(y, x)))) 
	     + (lambda * (lambda - 2) * (cos((lambda - 2) * get_angle(y, x)) * d_theta_dyd_theta_dy(x, y) + (lambda - 2) * d_theta_dy(x, y) * d_theta_dy(x, y) * (-1) * sin((lambda-2) * get_angle(y, x)) ));
}

double dv_rdxdv_rdy(double x, double y)
{
  return (v_F * lambda * (cos(lambda * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) + lambda * (-1) * d_theta_dx(x, y) * d_theta_dy(x, y) * sin(lambda * get_angle(y, x)))) 
	     + (lambda * (lambda - 2) * (cos((lambda - 2) * get_angle(y, x)) * d_theta_dxd_theta_dy(x, y) + (lambda - 2) * (-1) * d_theta_dx(x, y) * d_theta_dy(x, y) * sin((lambda-2) * get_angle(y, x)) ));
}

// \partial^2 u / \partial x^2 
double dudxdudx(double x, double y)
{
  return D * (drdxdrdx(x, y) * u_r(x, y) + 2 * drdx(x, y) * du_rdx(x, y) + r(x, y) * du_rdxdu_rdx(x, y) );
}

double dudydudy(double x, double y)
{
  return D * (drdydrdy(x, y) * u_r(x, y) + 2 * drdy(x, y) * du_rdy(x, y) + r(x, y) * du_rdydu_rdy(x, y) );
}

double dudxdudy(double x, double y)
{
  return D * (drdxdrdy(x, y) * u_r(x, y) + drdx(x, y) * du_rdy(x, y) + drdy(x, y) * du_rdx(x, y) + r(x, y) * du_rdxdu_rdy(x, y) );
}

// \partial^2 v / \partial x^2
double dvdxdvdx(double x, double y)
{
  return D * (drdxdrdx(x, y) * v_r(x, y) + 2 * drdx(x, y) * dv_rdx(x, y) + r(x, y) * dv_rdxdv_rdx(x, y) );
}

double dvdydvdy(double x, double y)
{
  return D * (drdydrdy(x, y) * v_r(x, y) + 2 * drdy(x, y) * dv_rdy(x, y) + r(x, y) * dv_rdydv_rdy(x, y) );
}

double dvdxdvdy(double x, double y)
{
  return D * (drdxdrdy(x, y) * v_r(x, y) + drdx(x, y) * dv_rdy(x, y) + drdy(x, y) * dv_rdx(x, y) + r(x, y) * dv_rdxdv_rdy(x, y) );
}


// Exact solution u(x,y) and its derivatives.
static double u_fn(double x, double y) 
{
  return D * r(x, y) * u_r(x, y);
//	  (u_F * cos(lambda * get_angle(y, x)) - lambda * cos((lambda - 2) * get_angle(y, x)));
}

static double u_fndd(double x, double y, double& dx, double& dy)
{
  dx = D * drdx(x, y) * (u_F * cos(lambda * get_angle(y, x)) - lambda * cos((lambda - 2) * get_angle(y, x))) +
       D * r(x, y) * (u_F * (-1) * lambda * sin(lambda * get_angle(y, x)) * d_theta_dx(x, y)) - 
       D * r(x, y) * (lambda * (-1) * (lambda - 2) * sin((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));

  dy = D * drdy(x, y) * (u_F * cos(lambda * get_angle(y, x)) - lambda * cos((lambda - 2) * get_angle(y, x))) +
       D * r(x, y) * (u_F * (-1) * lambda * sin(lambda * get_angle(y, x)) * d_theta_dy(x, y)) - 
       D * r(x, y) * (lambda * (-1) * (lambda - 2) * sin((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));

  return u_fn(x, y);
}

// Exact solution v(x,y) and its derivatives.
static double v_fn(double x, double y)
{
  return D * r(x, y) * (v_F * sin(lambda * get_angle(y, x)) + lambda * sin((lambda - 2) * get_angle(y, x)));
}

static double v_fndd(double x, double y, double& dx, double& dy)
{

  dx = D * drdx(x, y) * (v_F * sin(lambda * get_angle(y, x)) + lambda * sin((lambda - 2) * get_angle(y, x))) +
       D * r(x, y) * (v_F * lambda * cos(lambda * get_angle(y, x)) * d_theta_dx(x, y)) + 
       D * r(x, y) * (lambda * (lambda - 2) * cos((lambda - 2) * get_angle(y, x)) * d_theta_dx(x, y));

  dy = D * drdy(x, y) * (v_F * sin(lambda * get_angle(y, x)) + lambda * sin((lambda - 2) * get_angle(y, x))) +
       D * r(x, y) * (v_F * lambda * cos(lambda * get_angle(y, x)) * d_theta_dy(x, y)) + 
       D * r(x, y) * (lambda * (lambda - 2) * cos((lambda - 2) * get_angle(y, x)) * d_theta_dy(x, y));

  return v_fn(x, y);
}

