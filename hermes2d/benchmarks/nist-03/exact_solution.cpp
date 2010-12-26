double D = 1.0 / (2 * G);
double u_F = (k - Q * (lambda + 1));
double v_F = (k + Q * (lambda + 1));

static double theta(double x, double y)
{
  return asin(y/pow((x*x + y*y), (1.0/2.0)));
}

static double d_theta_dx(double x, double y)
{
  return y / (x * x + y * y);
}

static double d_theta_dy(double x, double y)
{
  return -x ;
}

static double r(double x, double y)
{
  return pow((x*x + y*y), (lambda/2.0));  // r^labbda
}

static double drdx(double x, double y)
{
  return lambda * x * pow((x*x + y*y), (lambda/2.0 - 1.0));
}

static double drdy(double x, double y)
{
  return lambda * y * pow((x*x + y*y), (lambda/2.0 - 1.0));
}

// Exact solution u(x,y) and its derivatives.
static double u_fn(double x, double y) 
{
  return D * r(x, y) * (u_F * cos(lambda * theta(x, y)) - lambda * cos((lambda - 2) * theta(x, y)));
}

static double u_fndd(double x, double y, double& dx, double& dy)
{
  dx = D * drdx(x, y) * (u_F * cos(lambda * theta(x, y)) - lambda * cos((lambda - 2) * theta(x, y))) +
       D * r(x, y) * (u_F * (-1) * lambda * sin(lambda * theta(x, y)) * d_theta_dx(x, y)) - 
       D * r(x, y) * (lambda * (-1) * (lambda - 2) * sin((lambda - 2) * theta(x, y)) * d_theta_dx(x, y));

  dy = D * drdy(x, y) * (u_F * cos(lambda * theta(x, y)) - lambda * cos((lambda - 2) * theta(x, y))) +
       D * r(x, y) * (u_F * (-1) * lambda * sin(lambda * theta(x, y)) * d_theta_dy(x, y)) - 
       D * r(x, y) * (lambda * (-1) * (lambda - 2) * sin((lambda - 2) * theta(x, y)) * d_theta_dy(x, y));

  return u_fn(x, y);
}

// Exact solution v(x,y) and its derivatives.
static double v_fn(double x, double y)
{
  return D * r(x, y) * (v_F * sin(lambda * theta(x, y)) + lambda * sin((lambda - 2) * theta(x, y)));
}

static double v_fndd(double x, double y, double& dx, double& dy)
{
  dx = D * drdx(x, y) * (v_F * sin(lambda * theta(x, y)) + lambda * sin((lambda - 2) * theta(x, y))) +
       D * r(x, y) * (v_F * lambda * cos(lambda * theta(x, y)) * d_theta_dx(x, y)) - 
       D * r(x, y) * (lambda * (lambda - 2) * sin((lambda - 2) * theta(x, y)) * d_theta_dx(x, y));

  dy = D * drdy(x, y) * (v_F * sin(lambda * theta(x, y)) + lambda * sin((lambda - 2) * theta(x, y))) +
       D * r(x, y) * (v_F * lambda * cos(lambda * theta(x, y)) * d_theta_dy(x, y)) - 
       D * r(x, y) * (lambda * (lambda - 2) * sin((lambda - 2) * theta(x, y)) * d_theta_dy(x, y));

  return v_fn(x, y);
}

