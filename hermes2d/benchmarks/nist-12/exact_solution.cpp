class ExactSolutionNIST12 : public ExactSolution1D
{
public:
  ExactSolutionNIST12(Mesh* mesh, double alpha_p, double x_p, double y_p, double alpha_w, double x_w, double y_w, double omega_c, double r_0, double epsilon)
    : ExactSolution1D(mesh), alpha_p(alpha_p), x_p(x_p), y_p(y_p), alpha_w(alpha_w), x_w(x_w), y_w(y_w), omega_c(omega_c), r_0(r_0), epsilon(epsilon) { };

  double get_angle(double y, double x) {
    double theta = atan2(y, x);
    if (theta < 0)
      theta += 2 * M_PI;
    return theta;
  }
  template<typename Real>
  Real rhs(Real x, Real y) {
    //For more elegant showing please execute file "generate_rhs.py" 

    Real a_P = (-alpha_p * pow((x - x_p), 2) - alpha_p * pow((y - y_p), 2));

    Real a_W = pow(x - x_w, 2);
    Real b_W = pow(y - y_w, 2);
    Real c_W = sqrt(a_W + b_W);
    Real d_W = ((alpha_w * x - (alpha_w * x_w)) * (2 * x - (2 * x_w)));
    Real e_W = ((alpha_w * y - (alpha_w * y_w)) * (2 * y - (2 * y_w)));
    Real f_W = (pow(alpha_w * c_W - (alpha_w * r_0), 2) + 1.0);
    Real g_W = (alpha_w * c_W - (alpha_w * r_0));

    return 4 * exp(a_P) * alpha_p * (alpha_p * (x - x_p) * (x - x_p) + alpha_p * (y - y_p) * (y - y_p) - 1)
           + ((alpha_w/(c_W * f_W)) - (d_W/(2 * pow(a_W + b_W, 1.5) * f_W)) - ((alpha_w * d_W * g_W)/((a_W + b_W) * pow(f_W, 2))) 
           + (alpha_w/(c_W * f_W)) - (e_W/(2 * pow(a_W + b_W, 1.5) * f_W)) - ((alpha_w * e_W * g_W)/((a_W + b_W) * pow(f_W, 2))))
           + (1.0 / epsilon) * (1.0 / epsilon) * exp(-(1 + y) / epsilon);  
  }

  double fn(double x, double y) {
    double ALPHA_C = (M_PI/ omega_c);

    return exp(-alpha_p * (pow((x - x_p), 2) + pow((y - y_p), 2)))
           + (pow(sqrt(x*x + y*y), ALPHA_C) * sin(ALPHA_C * get_angle(y, x)))
           + atan(alpha_w * (sqrt(pow(x - x_w, 2) + pow(y - y_w, 2)) - r_0))
           + exp(-(1 + y) / epsilon);
  };

  // Function representing an exact one-dimension valued solution.
  virtual scalar exact_function (double x, double y, scalar& dx, scalar& dy) {
    double a_P = -alpha_p * ( (x - x_p) * (x - x_p) + (y - y_p) * (y - y_p));

    double ALPHA_C = (M_PI/ omega_c);
    double a_C = sqrt(x*x + y*y);
    double b_C = pow(a_C, (ALPHA_C - 1.0));
    double c_C = pow(a_C, ALPHA_C);
    double d_C = ((y*y)/(x*x) + 1.0 );

    double a_W = pow(x - x_w, 2);
    double b_W = pow(y - y_w, 2);
    double c_W = sqrt(a_W + b_W);
    double d_W = (alpha_w * x - (alpha_w * x_w));
    double e_W = (alpha_w * y - (alpha_w * y_w));
    double f_W = (pow(alpha_w * c_W - (alpha_w * r_0), 2) + 1.0);

    dx = -exp(a_P) * (2 * alpha_p * (x - x_p))
         + (((ALPHA_C* x* sin(ALPHA_C * get_angle(y,x)) *b_C)/a_C) - ((ALPHA_C *y *cos(ALPHA_C * get_angle(y, x)) * c_C)/(pow(x, 2.0) *d_C)))
         + (d_W / (c_W * f_W));
    dy = -exp(a_P) * (2 * alpha_p * (y - y_p))
         + (((ALPHA_C* cos(ALPHA_C* get_angle(y, x)) *c_C)/(x * d_C)) + ((ALPHA_C* y* sin(ALPHA_C* get_angle(y, x)) *b_C)/a_C))
         + (e_W / (c_W * f_W))
         + (-1) * (1.0 / epsilon) * exp(-(1 + y) / epsilon); 

    return fn(x, y);
  };

  // Members.
  double alpha_p;
  double x_p;
  double y_p;
  double alpha_w;
  double x_w;
  double y_w;
  double omega_c;
  double r_0;
  double epsilon;
};