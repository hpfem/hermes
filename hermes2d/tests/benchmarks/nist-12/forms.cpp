template<typename Real, typename Scalar>
Scalar bilinear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_grad_u_grad_v<Real, Scalar>(n, wt, u, v);
}

template<typename Real>
Real rhs(Real x, Real y)
{
  //For more elegant showing please execute file "generate_rhs.py" 

  return 4*sin(2*atan(y/x)/3)/(3*pow((pow(x,2) + pow(y,2)),(2.0/3.0))) - 8*pow(x,2)*sin(2*atan(y/x)/3)/(9*pow((pow(x,2) + pow(y,2)),(5.0/3.0))) - 8*pow(y,2)*sin(2*atan(y/x)/3)/(9*pow((pow(x,2) + pow(y,2)),(5.0/3.0))) + 80000*pow(x,2)*(150 - 200*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0)))/(pow((1 + pow((150 - 200*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0))),2)),2)*(pow(x,2) + pow((3.0/4.0 + y),2))) + 80000*pow((3.0/4.0 + y),2)*(150 - 200*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0)))/(pow((1 + pow((150 - 200*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0))),2)),2)*(pow(x,2) + pow((3.0/4.0 + y),2))) - 4*pow((pow(x,2) + pow(y,2)),(1.0/3.0))*sin(2*atan(y/x)/3)/(9*pow(x,2)*pow((1 + pow(y,2)/pow(x,2)),2)) - 4*y*pow((pow(x,2) + pow(y,2)),(1.0/3.0))*cos(2*atan(y/x)/3)/(3*pow(x,3)*pow((1 + pow(y,2)/pow(x,2)),2)) - 4*pow(y,3)*pow((pow(x,2) + pow(y,2)),(1.0/3.0))*cos(2*atan(y/x)/3)/(3*pow(x,5)*pow((1 + pow(y,2)/pow(x,2)),2)) - 4*pow(y,2)*pow((pow(x,2) + pow(y,2)),(1.0/3.0))*sin(2*atan(y/x)/3)/(9*pow(x,4)*pow((1 + pow(y,2)/pow(x,2)),2)) + 4*y*pow((pow(x,2) + pow(y,2)),(1.0/3.0))*cos(2*atan(y/x)/3)/(3*pow(x,3)*(1 + pow(y,2)/pow(x,2))) - 200*pow(x,2)/((1 + pow((150 - 200*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0))),2))*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(3.0/2.0))) - 200*pow((3.0/4.0 + y),2)/((1 + pow((150 - 200*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0))),2))*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(3.0/2.0))) + 400/((1 + pow((150 - 200*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0))),2))*pow((pow(x,2) + pow((3.0/4.0 + y),2)),(1.0/2.0))) + 10000*exp(-100 - 100*y) - 4000*exp(-1000*pow((1.0/4.0 + y),2) - 1000*pow((x - pow(5,(1.0/2.0))/4),2)) + pow((500 + 2000*y),2)*exp(-1000*pow((1.0/4.0 + y),2) - 1000*pow((x - pow(5,(1.0/2.0))/4),2)) + pow((-2000*x + 500*pow(5,(1.0/2.0))),2)*exp(-1000*pow((1.0/4.0 + y),2) - 1000*pow((x - pow(5,(1.0/2.0))/4),2));
}

template<typename Real, typename Scalar>
Scalar linear_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return -int_F_v<Real, Scalar>(n, wt, rhs, v, e);
}
