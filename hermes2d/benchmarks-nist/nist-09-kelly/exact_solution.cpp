static double fn(double x, double y)
{
  if (PROB_PARAM == 0){
   ALPHA = 20;
   X_LOC = -0.05;
   Y_LOC = -0.05;
   R_ZERO = 0.7;
   }
else if (PROB_PARAM == 1){
   ALPHA = 1000;
   X_LOC = -0.05;
   Y_LOC = -0.05;
   R_ZERO = 0.7;
   }
else if (PROB_PARAM == 2){
   ALPHA = 1000;
   X_LOC = 1.5;
   Y_LOC = 0.25;
   R_ZERO = 0.92;
   }
else{
   ALPHA = 50;
   X_LOC = 0.5;
   Y_LOC = 0.5;
   R_ZERO = 0.25;
   }
 
  return atan(ALPHA * (sqrt(pow(x - X_LOC, 2) + pow(y - Y_LOC, 2)) - R_ZERO));
}

static double fndd(double x, double y, double& dx, double& dy)
{ 
  if (PROB_PARAM == 0){
   ALPHA = 20;
   X_LOC = -0.05;
   Y_LOC = -0.05;
   R_ZERO = 0.7;
   }
else if (PROB_PARAM == 1){
   ALPHA = 1000;
   X_LOC = -0.05;
   Y_LOC = -0.05;
   R_ZERO = 0.7;
   }
else if (PROB_PARAM == 2){
   ALPHA = 1000;
   X_LOC = 1.5;
   Y_LOC = 0.25;
   R_ZERO = 0.92;
   }
else{
   ALPHA = 50;
   X_LOC = 0.5;
   Y_LOC = 0.5;
   R_ZERO = 0.25;
   }

  double a = pow(x - X_LOC, 2);
  double b = pow(y - Y_LOC, 2);
  double c = sqrt(a + b);
  double d = (ALPHA*x - (ALPHA * X_LOC));
  double e = (ALPHA*y - (ALPHA * Y_LOC));
  double f = (pow(ALPHA*c - (ALPHA * R_ZERO), 2) + 1.0);

  dx = (d/(c * f));
  dy = (e/(c * f));

  return fn(x, y);
}
