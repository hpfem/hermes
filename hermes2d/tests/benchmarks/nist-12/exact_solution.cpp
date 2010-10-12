//double get_angle(double y, double x )<--  ignore this for now, but might come in handy later if test doesn't work out.
//{
//  double theta = atan2(y, x);
//  if (theta < 0)
//    theta += 2 * M_PI;
//  return theta;
//}


static double fn(double x, double y)
{
  double a = pow(x - X_w, 2);
  double b = pow(y - Y_w, 2);
  double c = pow(x - X_p, 2);
  double d = pow(y - Y_p, 2);

  return (pow(sqrt(x*x + y*y), M_PI/OMEGA)* sin(atan2(y, x)* (M_PI/OMEGA)) + atan(ALPHA_w* (sqrt(a+b) - R_0)) + exp(-ALPHA_p* (c+d)) + exp(-(1+y)/ EPSILON));
}

static double fndd(double x, double y, double& dx, double& dy)
{
  //For a more elegant showing please execute file "generate_diff_f_x.py" or "generate_diff_f_y.py"

  dx = 2*x*sin(2*atan(y/x)/3)/(3*pow((pow(x,2) + pow(y,2)),(2.0/3.0))) - 2*y*pow((pow(x,2) + pow(y,2)),(1.0/3.0))*cos(2*atan(y/x)/3)/(3*pow(x,2)*(1 + pow(y,2)/pow(x,2))) + 200*x/((1 + pow((150 - 200*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0))),2))*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0))) + (-2000*x + 500*pow(5,(1.0/2.0)))*exp(-1000*pow((1.0/4.0 + y),2) - 1000*pow((x - pow(5,(1.0/2.0))/4),2));

  dy = 200*(3.0/4.0 + y)/((1 + pow((150 - 200*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0))),2))*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0))) + 2*y*sin(2*atan(y/x)/3)/(3*pow((pow(x,2) + pow(y,2)),(2.0/3.0))) + 2*pow((pow(x,2) + pow(y,2)),(1.0/3.0))*cos(2*atan(y/x)/3)/(3*x*(1 + pow(y,2)/pow(x,2))) - 100*exp(-100 - 100*y) - (500 + 2000*y)*exp(-1000*pow((1.0/4.0 + y),2) - 1000*pow((x - pow(5,(1.0/2.0))/4),2));

  return fn(x, y);
}
