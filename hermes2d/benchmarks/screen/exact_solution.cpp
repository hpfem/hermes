// Fresnel integral
extern "C" void fresnl( double xxa, double *ssa, double *cca );

// exact solution
scalar Fn(double u)
{
  double s, c;
  fresnl(sqrt(2/M_PI) * u, &s , &c);
  scalar fres = cplx(c,-s);
  scalar a = cplx(0.0, M_PI/4);
  scalar b = cplx(0.0, u*u);
  return 0.5*sqrt(M_PI) * exp(b) * (exp(-a) - sqrt(2.0)*(fres));
}

scalar Fder(double u)
{
  scalar a = cplx(0.0, M_PI/4);
  scalar b = cplx(0.0, u*u);
  scalar d = cplx(0.0, 2.0*u);
  double s, c;
  fresnl(sqrt(2/M_PI) * u, &s , &c);
  scalar fres = cplx(c,-s);
  scalar fresder = exp(-b);

  return 0.5*sqrt(M_PI) * exp(b) * ( d * (exp(-a) - sqrt(2.0)*(fres)) - sqrt(2.0)*fresder*sqrt(2.0/M_PI) );
}

scalar Fder2(double u)
{
  scalar a = cplx(0.0, M_PI/4);
  scalar i = cplx(0.0,1.0);
  scalar b = cplx(0.0, u*u);
  scalar d = cplx(0.0, 2.0*u);
  double s, c;
  fresnl(sqrt(2/M_PI) * u, &s , &c);
  scalar fres = cplx(c,-s);
  scalar fresder = exp(-b);
  scalar fresder2 = exp(-b)*(-2.0 * i * u);

  return 2.0 * u * i * Fder(u) +
         0.5 * sqrt(M_PI) * exp(b) *
          ( 2.0 * i * (exp(-a) - sqrt(2.0)*(fres)) + d * (-sqrt(2.0)*fresder*sqrt(2.0/M_PI)) - sqrt(2.0) * fresder2 * sqrt(2.0/M_PI) );
}

scalar der_Hr(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = cplx(0.0, M_PI/4 - k*r);
  scalar i = cplx(0.0,1.0);
  return 1/sqrt(M_PI) * exp(a) *
        ( (-i*k)*(Fn(sqrt(2*k*r)*sin(t/2 - M_PI/8)) + Fn(sqrt(2*k*r)*sin(t/2 + M_PI/8))) +
        (Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8))*(sqrt(k)/sqrt(2*r)*sin(t/2 - M_PI/8)) +
         Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8))*(sqrt(k)/sqrt(2*r)*sin(t/2 + M_PI/8))));
}

scalar der_Hrr(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = cplx(0.0, M_PI/4 - k*r);
  scalar i = cplx(0.0,1.0);
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar b1 = (sqrt(k/(2*r))*sin(t/2 - M_PI/8));
  scalar b2 = (sqrt(k/(2*r))*sin(t/2 + M_PI/8));
  return -i*k*der_Hr(x,y) + 1/sqrt(M_PI) * exp(a) *
        ( (-i*k)*(f1_d*b1 + f2_d*b2) +
        ( f1_d2*b1*b1 + f2_d2*b2*b2) +
          f1_d*(-0.5*sqrt(k/(2*r*r*r))*sin(t/2 - M_PI/8))  + f2_d*(-0.5*sqrt(k/(2*r*r*r))*sin(t/2 + M_PI/8)));
}

scalar der_Hrt(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar i = cplx(0.0,1.0);
  scalar a = cplx(0.0, M_PI/4 - k*r);
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar b1 = (sqrt(k)/sqrt(2*r)*sin(t/2 - M_PI/8));
  scalar b2 = (sqrt(k)/sqrt(2*r)*sin(t/2 + M_PI/8));
  scalar c1 = (sqrt(k*r)/sqrt(2.0)*cos(t/2 - M_PI/8));
  scalar c2 = (sqrt(k*r)/sqrt(2.0)*cos(t/2 + M_PI/8));
  return 1/sqrt(M_PI) * exp(a) *
        ( (-i*k)*(f1_d*c1 + f2_d*c2) +
        ( f1_d2*b1*c1 + f2_d2*b2*c2) +
          f1_d*(0.5*sqrt(k/(2*r))*cos(t/2 - M_PI/8))  + f2_d*(0.5*sqrt(k/(2*r))*cos(t/2 + M_PI/8)));
}

scalar der_Ht(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = cplx(0.0, M_PI/4 - k*r);
  return 1/sqrt(M_PI) * exp(a) *
         (Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8))*(sqrt(k*r/2)*cos(t/2 - M_PI/8)) +
          Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8))*(sqrt(k*r/2)*cos(t/2 + M_PI/8)));
}

scalar der_Htr(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar i = cplx(0.0,1.0);
  scalar a = cplx(0.0, M_PI/4 - k*r);
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar b1 = (sqrt(k)/sqrt(2*r)*sin(t/2 - M_PI/8));
  scalar b2 = (sqrt(k)/sqrt(2*r)*sin(t/2 + M_PI/8));
  scalar c1 = (sqrt(k*r)/sqrt(2.0)*cos(t/2 - M_PI/8));
  scalar c2 = (sqrt(k*r)/sqrt(2.0)*cos(t/2 + M_PI/8));
  return -i*k*der_Ht(x,y) + 1/sqrt(M_PI) * exp(a) *
         ((f1_d2*b1*c1 + f2_d2*b2*c2) +
          f1_d*(0.5*sqrt(k/(2*r))*cos(t/2 - M_PI/8))  + f2_d*(0.5*sqrt(k/(2*r))*cos(t/2 + M_PI/8)));
}

scalar der_Htt(double x, double y)
{
  double r = sqrt(x*x + y*y);
  double t = atan2(y,x);
  scalar a = cplx(0.0, M_PI/4 - k*r);
  scalar f1_d = Fder(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d = Fder(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar f1_d2 = Fder2(sqrt(2*k*r)*sin(t/2 - M_PI/8));
  scalar f2_d2 = Fder2(sqrt(2*k*r)*sin(t/2 + M_PI/8));
  scalar c1 = (sqrt(k*r/(2))*cos(t/2 - M_PI/8));
  scalar c2 = (sqrt(k*r/(2))*cos(t/2 + M_PI/8));
  return 1/sqrt(M_PI) * exp(a) *
         ((f1_d2*c1*c1 + f2_d2*c2*c2) +
          f1_d*(-0.5*sqrt(k*r/2)*sin(t/2 - M_PI/8))  + f2_d*(-0.5*sqrt(k*r/2)*sin(t/2 + M_PI/8)));
}

scalar exact0(double x, double y, scalar& dx, scalar& dy)
{
  double r = sqrt(x*x + y*y);
  double theta = atan2(y,x);
  scalar Hr = der_Hr(x,y);
  scalar Ht = der_Ht(x,y);
  scalar i = cplx(0.0,1.0);
  return  -i * (Hr * y/r + Ht * x/(r*r));
}

scalar exact1(double x, double y, scalar& dx, scalar& dy)
{
  double r = sqrt(x*x + y*y);
  double theta = atan2(y,x);
  scalar Hr = der_Hr(x,y);
  scalar Ht = der_Ht(x,y);
  scalar i = cplx(0.0,1.0);
  return  i * ( Hr * x/r - Ht * y/(r*r));
}

static void exact_sol(double x, double y, scalar& u0, scalar& u1, scalar& u1dx, scalar& u0dy)
{
  scalar dx,dy;
  u0 = exact0(x,y,dx,dy);
  u1 = exact1(x,y,dx,dy);

  scalar Hr = der_Hr(x,y);
  scalar Ht = der_Ht(x,y);
  scalar Hrr = der_Hrr(x,y);
  scalar Hrt = der_Hrt(x,y);
  scalar Htr = der_Htr(x,y);
  scalar Htt = der_Htt(x,y);

  double r = sqrt(x*x + y*y);
  double theta = atan2(y,x);
  scalar i = cplx(0.0,1.0);

  u1dx =  i * (( Hrr * x/r + Hrt * (-y/(r*r))) * x/r     + Hr * (y*y)/(r*r*r) -
               ((Htr * x/r + Htt * (-y/(r*r))) * y/(r*r) + Ht * (-2.0*x*y/(r*r*r*r))));
  u0dy = -i * (( Hrr * y/r + Hrt *   x/(r*r))  * y/r     + Hr * (x*x)/(r*r*r) +
                (Htr * y/r + Htt *   x/(r*r))  * x/(r*r) + Ht * (-2.0*x*y/(r*r*r*r)));
}

scalar2& exact(double x, double y, scalar2& dx, scalar2& dy)
{
  static scalar2 ex;
  exact_sol(x,y, ex[0], ex[1], dx[1], dy[0]);

  dx[0] = 0.0; // not important
  dy[1] = 0.0; // not important

  return ex;
}
