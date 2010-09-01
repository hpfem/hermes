// K (van Genuchten).
double K(double h)
{
  double alpha;
  double n;
  double m;
  alpha = ALPHA;
  n = N;
  m = M;
  if (h < 0) return 
      K_S*pow((1 + pow((-alpha*h),n)),(-m/2))*pow((1 -
     pow((-alpha*h),(m*n))*pow((1 + pow((-alpha*h),n)),(-m))),2) ;
  else return K_S;    
}

// dK/dh (van Genuchten).
double dKdh(double h)
{
  double alpha;
  double n;
  double m;
  alpha = ALPHA;
  n = N;
  m = M;
  if (h < 0) return 
	K_S*pow((1 + pow((-alpha*h),n)),(-m/2))*(1 -
	pow((-alpha*h),(m*n))*pow((1 +
	pow((-alpha*h),n)),(-m)))*(-2*m*n*pow((-alpha*h),(m*n))*pow((1 +
	pow((-alpha*h),n)),(-m))/h +
	2*m*n*pow((-alpha*h),n)*pow((-alpha*h),(m*n))*pow((1 +
	pow((-alpha*h),n)),(-m))/(h*(1 + pow((-alpha*h),n)))) -
	K_S*m*n*pow((-alpha*h),n)*pow((1 + pow((-alpha*h),n)),(-m/2))*pow((1 -
	pow((-alpha*h),(m*n))*pow((1 + pow((-alpha*h),n)),(-m))),2)/(2*h*(1 +
	pow((-alpha*h),n))) ;
  else return 0;
}

// ddK/dhh (van Genuchten).
double ddKdhh(double h)
{
  double alpha;
  double n;
  double m;
  alpha = ALPHA;
  n = N;
  m = M;
  if (h < 0) return 
      K_S*pow((1 + pow((-alpha*h),n)),(-m/2))*(1 -
      pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m)))*(-2*pow(m,2)*pow(n,2)*pow((-alpha*h),(m*n))*pow((1
      + pow((-alpha*h),n)),(-m))/pow(h,2) +
      2*m*n*pow((-alpha*h),(m*n))*pow((1 + pow((-alpha*h),n)),(-m))/pow(h,2)
      - 2*m*n*pow((-alpha*h),n)*pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))/(pow(h,2)*(1 + pow((-alpha*h),n))) -
      2*m*pow(n,2)*pow((-alpha*h),(2*n))*pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))/(pow(h,2)*pow((1 + pow((-alpha*h),n)),2)) -
      2*pow(m,2)*pow(n,2)*pow((-alpha*h),(2*n))*pow((-alpha*h),(m*n))*pow((1
      + pow((-alpha*h),n)),(-m))/(pow(h,2)*pow((1 + pow((-alpha*h),n)),2)) +
      2*m*pow(n,2)*pow((-alpha*h),n)*pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))/(pow(h,2)*(1 + pow((-alpha*h),n))) +
      4*pow(m,2)*pow(n,2)*pow((-alpha*h),n)*pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))/(pow(h,2)*(1 + pow((-alpha*h),n)))) +
      K_S*pow((1 + pow((-alpha*h),n)),(-m/2))*(-m*n*pow((-alpha*h),(m*n))*pow((1
      + pow((-alpha*h),n)),(-m))/h +
      m*n*pow((-alpha*h),n)*pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))/(h*(1 +
      pow((-alpha*h),n))))*(-2*m*n*pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))/h +
      2*m*n*pow((-alpha*h),n)*pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))/(h*(1 + pow((-alpha*h),n)))) -
      K_S*m*n*pow((-alpha*h),n)*pow((1 + pow((-alpha*h),n)),(-m/2))*(1 -
      pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m)))*(-2*m*n*pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))/h +
      2*m*n*pow((-alpha*h),n)*pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))/(h*(1 + pow((-alpha*h),n))))/(h*(1 +
      pow((-alpha*h),n))) + K_S*m*n*pow((-alpha*h),n)*pow((1 +
      pow((-alpha*h),n)),(-m/2))*pow((1 - pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))),2)/(2*pow(h,2)*(1 + pow((-alpha*h),n))) +
      K_S*m*pow(n,2)*pow((-alpha*h),(2*n))*pow((1 +
      pow((-alpha*h),n)),(-m/2))*pow((1 - pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))),2)/(2*pow(h,2)*pow((1 +
      pow((-alpha*h),n)),2)) - K_S*m*pow(n,2)*pow((-alpha*h),n)*pow((1 +
      pow((-alpha*h),n)),(-m/2))*pow((1 - pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))),2)/(2*pow(h,2)*(1 + pow((-alpha*h),n))) +
      K_S*pow(m,2)*pow(n,2)*pow((-alpha*h),(2*n))*pow((1 +
      pow((-alpha*h),n)),(-m/2))*pow((1 - pow((-alpha*h),(m*n))*pow((1 +
      pow((-alpha*h),n)),(-m))),2)/(4*pow(h,2)*pow((1 +
      pow((-alpha*h),n)),2)) ;

  else return 0;
}

// C (van Genuchten).
double C(double h)
{
  if (h < 0) return 
   STORATIVITY*pow((1 + pow((-ALPHA*h),N)),(-M))*(THETA_S - THETA_R)/THETA_S - 
   M*N*pow((-ALPHA*h),N)*pow((1 + pow((-ALPHA*h),N)),(-M))*(THETA_S - THETA_R)/(h*(1 + pow((-ALPHA*h),N)));
  else return STORATIVITY;    
}

// dC/dh (van Genuchten).
double dCdh(double h)
{
  if (h < 0) return
    M*N*pow((-ALPHA*h),N)*pow((1 + pow((-ALPHA*h),N)),(-M))*(THETA_S - THETA_R)/(pow(h,2)*(1 + pow((-ALPHA*h),N))) + M*pow(N,2)*pow((-ALPHA*h),(2*N))*pow((1 + 
    pow((-ALPHA*h),N)),(-M))*(THETA_S - THETA_R)/(pow(h,2)*pow((1 + pow((-ALPHA*h),N)),2)) + pow(M,2)*pow(N,2)*pow((-ALPHA*h),(2*N))*pow((1 + pow((-ALPHA*h),N)),(-M))*(THETA_S - THETA_R)/(pow(h,2)*
    pow((1 + pow((-ALPHA*h),N)),2)) - M*pow(N,2)*pow((-ALPHA*h),N)*pow((1 + pow((-ALPHA*h),N)),(-M))*(THETA_S - THETA_R)/(pow(h,2)*(1 + pow((-ALPHA*h),N))) - 
    M*N*STORATIVITY*pow((-ALPHA*h),N)*pow((1 + pow((-ALPHA*h),N)),(-M))*(THETA_S - THETA_R)/(THETA_S*h*(1 + pow((-ALPHA*h),N))) ;
  else return 0;    
}
