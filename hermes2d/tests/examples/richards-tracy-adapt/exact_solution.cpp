// value of the boundary condition along the upper edge
double boundary_value(double x, double& dhdx) {
//  a...intent(in)... material parameter
//  H_R.intent(in)... pressure head if material is very dry
//  x.. intent(in)... x-coordinate of the boundary point
//  dhdx..intent(out).. derivative of the boundary condition value
//  h...intent(out)...  value of the boundary condition
  dhdx = (1-exp(ALPHA*H_R)*M_PI*cos(M_PI*x/ALPHA))/(ALPHA*ALPHA*(exp(A*H_R)
         +(1-exp(ALPHA*H_R))*sin(M_PI*x/ALPHA))); 
  return 1/ALPHA*log(exp(ALPHA*H_R) + (1-exp(ALPHA*H_R))*sin(M_PI*x/ALPHA)); 
}

// exact solution and its derivative
double exact_sol(double x, double z, double& dhdx, double& dhdz) {

// !--------i/o variables-------------------------
//       !> coordinates of the desired point
//       real(kind=rkind), dimension(:), intent(in) :: x,z
//       !> simulation time to plot
//       real(kind=rkind), intent(in)               :: TIME
//       !> material parameter
//       real(kind=rkind), dimension(:), intent(in) :: ALPHA
//       !> "a very dry" media pressure head
//       real(kind=rkind, intent(in)                :: H_R
//       !> saturated water content(porosity) and residual water content
//       real(kind=rkind), intent(in)               :: THETA_R, THETA_S
//       !> saturated hydraulic conductivity
//       real(kind=rkind), intent(in)               :: K_S
//       !> domain length(x_3) + width(x_1)  
//       real(kind=rkind), intent(in)               :: L, A 
//       !> pressure head at the desired point
//       real(kind=rkind), intent(out)              :: h
//       !> solution derivations
//       real(kind=rkind), intent(out)              :: dhdx, dhdz

  double reps = 1e-9 ; // computer accuracy
  double ho = 1 - exp(ALPHA*H_R) ;
  double beta = sqrt(ALPHA*ALPHA/4 + (M_PI/A)*(M_PI/A)) ;
  double c = ALPHA*(THETA_S - THETA_R)/K_S ;
  double hss = ho*sin(M_PI*x/A)*exp(ALPHA/2*(L-z))*sinh(beta*z)/sinh(beta*L) ;
  double sum = 0 ;

  // fourier series for the function value
  for (int i=1; i >= 0; i++) {
    double lambda = i*M_PI/L ;
    double gamma = 1/c*(beta*beta + lambda*lambda) ;
    double tmp = pow((double)-1,(double)i)*lambda/gamma*sin(lambda*z)*exp(-gamma*TIME) ;
    sum += tmp ;
    if (fabs(1/c*lambda/gamma * exp(-gamma*TIME)) < reps*reps && i > 100) break;
  }  

  double phi = 2*ho/(L*c)*sin(M_PI*x/A)*exp(ALPHA/2*(L-z))*sum ;
  double hbar = phi + hss ;

  // result for the function value
  double h = 1/ALPHA*log(exp(ALPHA*H_R)+hbar) ; // !!!!! function h !!!!

  double a_const = exp(ALPHA*H_R) ;
  double f_xz = 2*ho/(L*c)*sin(M_PI*x/A)*exp(ALPHA/2*(L-z)) ; 
  double g_xz = sum ;
  double k_xz = ho*sin(M_PI*x/A)*exp(ALPHA/2*(L-z))*sinh(beta*z)/sinh(beta*L) ;    
  double dfdx = 2*exp(1/2*ALPHA*(L-z))*ho*M_PI*cos(M_PI*x/A)/(A*c*L) ;
  double dfdz = -ALPHA*exp(1/2*ALPHA*(L-z))*ho*sin(M_PI*x/A)/(c*L) ;
  double sum2 = 0 ;
     
  // fourier series for the derivatives
  for (int i=1; i >= 0; i++) {
    double lambda = i*M_PI/L;
    double gamma = 1/c*(beta*beta + lambda*lambda);
    double tmp  = pow((double)-1,(double)i)*i*exp(-gamma*TIME)*lambda*lambda*cos(lambda*z)/gamma;
    sum2 += tmp;
    if (fabs(lambda*lambda*exp(-gamma*TIME)/gamma) < reps*reps && i > 100) break;
  }

  double dgdz = sum2 ;
  double dgdx = 0 ;
  double dkdx = exp(1/2*ALPHA*(L-z))*ho*M_PI*cos(M_PI*x/A)*sinh(beta*z)/(A*sinh(beta*L)) ;
  double dkdz = beta*exp(1/2*ALPHA*(L-z))*ho*cosh(beta*z)*sin(M_PI*x/A)/sinh(beta*L)
                - ALPHA*exp(1/2*ALPHA*(L-z))*ho*sin(M_PI*x/A)*sinh(beta*z)/(2*sinh(beta*L)) ;
    
  // result for derivatives
  dhdx = 1/ALPHA*1/(a_const + f_xz*g_xz +k_xz)*(dfdx*g_xz + dgdx*f_xz + dkdx) ;
  dhdz = 1/ALPHA*1/(a_const + f_xz*g_xz +k_xz)*(dfdz*g_xz + dgdz*f_xz + dkdz) ;
    
  return h;
}
