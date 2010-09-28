double get_angle(double y, double x)
{
  double theta = atan2(y, x);
  if (theta < 0)
    theta += 2 * M_PI;
  return theta;
}

static double fn(double x, double y)
{
  if (PROB_PARAM == 0){
   OMEGA = ((5.0 * M_PI)/ 4.0);
   ALPHA = (M_PI/ OMEGA);
   }
else if (PROB_PARAM == 1){
   OMEGA = ((3.0 * M_PI)/ 2.0);
   ALPHA = (M_PI/ OMEGA);
   }
else if (PROB_PARAM == 2){
   OMEGA = ((7.0 * M_PI)/ 4.0);
   ALPHA = (M_PI/ OMEGA);
   }
else{
   OMEGA = (2.0 * M_PI);
   ALPHA = (M_PI/ OMEGA);
   }

  return (pow(sqrt(x*x + y*y), ALPHA) * sin(ALPHA * get_angle(y, x)));
}

static double fndd(double x, double y, double& dx, double& dy)
{
  if (PROB_PARAM == 0){
   OMEGA = ((5.0 * M_PI)/ 4.0);
   ALPHA = (M_PI/ OMEGA);
   }
else if (PROB_PARAM == 1){
   OMEGA = ((3.0 * M_PI)/ 2.0);
   ALPHA = (M_PI/ OMEGA);
   }
else if (PROB_PARAM == 2){
   OMEGA = ((7.0 * M_PI)/ 4.0);
   ALPHA = (M_PI/ OMEGA);
   }
else{
   OMEGA = (2.0 * M_PI);
   ALPHA = (M_PI/ OMEGA);
   }

  double a = sqrt(x*x + y*y);
  double b = pow(a, (ALPHA - 1.0));
  double c = pow(a, ALPHA);
  double d = ((y*y)/(x*x) + 1.0 );

  dx = (((ALPHA* x* sin(ALPHA * get_angle(y,x)) *b)/a) - ((ALPHA *y *cos(ALPHA * get_angle(y, x)) * c)/(pow(x, 2.0) *d)));
  dy = (((ALPHA* cos(ALPHA* get_angle(y, x)) *c)/(x * d)) + ((ALPHA* y* sin(ALPHA* get_angle(y, x)) *b)/a));

  return fn(x, y);
}

