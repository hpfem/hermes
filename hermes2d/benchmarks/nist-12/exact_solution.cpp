double get_angle(double y, double x)
{
  double theta = atan2(y, x);
  if (theta < 0)
    theta += 2 * M_PI;
  return theta;
}

static double fn(double x, double y)
{
  double ALPHA_C = (M_PI/ OMEGA_C);

  return exp(-ALPHA_P * (pow((x - X_P), 2) + pow((y - Y_P), 2)))
         + (pow(sqrt(x*x + y*y), ALPHA_C) * sin(ALPHA_C * get_angle(y, x)))
         + atan(ALPHA_W * (sqrt(pow(x - X_W, 2) + pow(y - Y_W, 2)) - R_0))
         + exp(-(1 + y) / EPSILON);
}

static double fndd(double x, double y, double& dx, double& dy)
{
  //For a more elegant showing please execute file "generate_diff_f_x.py" or "generate_diff_f_y.py"
  double a_P = -ALPHA_P * ( (x - X_P) * (x - X_P) + (y - Y_P) * (y - Y_P));

  double ALPHA_C = (M_PI/ OMEGA_C);
  double a_C = sqrt(x*x + y*y);
  double b_C = pow(a_C, (ALPHA_C - 1.0));
  double c_C = pow(a_C, ALPHA_C);
  double d_C = ((y*y)/(x*x) + 1.0 );

  double a_W = pow(x - X_W, 2);
  double b_W = pow(y - Y_W, 2);
  double c_W = sqrt(a_W + b_W);
  double d_W = (ALPHA_W * x - (ALPHA_W * X_W));
  double e_W = (ALPHA_W * y - (ALPHA_W * Y_W));
  double f_W = (pow(ALPHA_W * c_W - (ALPHA_W * R_0), 2) + 1.0);

  dx = -exp(a_P) * (2 * ALPHA_P * (x - X_P))
       + (((ALPHA_C* x* sin(ALPHA_C * get_angle(y,x)) *b_C)/a_C) - ((ALPHA_C *y *cos(ALPHA_C * get_angle(y, x)) * c_C)/(pow(x, 2.0) *d_C)))
       + (d_W / (c_W * f_W));
  dy = -exp(a_P) * (2 * ALPHA_P * (y - Y_P))
       + (((ALPHA_C* cos(ALPHA_C* get_angle(y, x)) *c_C)/(x * d_C)) + ((ALPHA_C* y* sin(ALPHA_C* get_angle(y, x)) *b_C)/a_C))
       + (e_W / (c_W * f_W))
       + (-1) * (1.0 / EPSILON) * exp(-(1 + y) / EPSILON); 

  return fn(x, y);
}
