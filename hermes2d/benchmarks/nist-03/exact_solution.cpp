double D = 1.0 / (2 * G);
double u_F = (k - Q * (lambda + 1));
double v_F = (k + Q * (lambda + 1));

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

double d_theta_dy(double x, double y)
{
  return x/(x*x + y*y) ;
}

double r(double x, double y)
{
  return pow((x*x + y*y), (lambda/2.0));  // r^labbda
}

double drdx(double x, double y)
{
  return lambda * x * pow((x*x + y*y), (lambda/2.0 - 1.0));
}

double drdy(double x, double y)
{
  return lambda * y * pow((x*x + y*y), (lambda/2.0 - 1.0));
}

// Exact solution u(x,y) and its derivatives.
static double u_fn(double x, double y) 
{
  return D * r(x, y) * (u_F * cos(lambda * get_angle(y, x)) - lambda * cos((lambda - 2) * get_angle(y, x)));
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

