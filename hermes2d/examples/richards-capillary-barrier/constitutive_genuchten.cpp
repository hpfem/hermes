
double horner(double *pol, double x, int n){
  double px=0.0;
  for (int i=0; i<n; i++){
//     printf("polhor %i \n", (n-1-i));
    px = px*x + pol[n-1-i] ;
  }
//   printf("%lf ", px);
  return px;
}

// K (van Genuchten).
double K(double h, int layer)
{
  double value ;
  int location ;
  if (h > LOW_LIMIT && h < 0 && POLYNOMIALS_READY == true){
// 	  double tmp=0.0;
// 	  for (int i=0; i<8; i++){
// 	    tmp += POLYNOMIALS[0][layer][i]*pow(h,i) ;
// 	  }
//           printf("%lf %lf %i %lf %lf \n", horner(POLYNOMIALS[0][layer], h, (6+NUM_OF_INSIDE_PTS)), h, tmp);
	   
//     printf("pred horner %i  \n", layer);
    return horner(POLYNOMIALS[layer][0], h, (6+NUM_OF_INSIDE_PTS));
//     return 0;

  }
  else
  {
    if (CONSTITUTIVE_TABLES_READY == false || h < -15000.0)  {
      ALPHA = ALPHA_vals[layer] ;
      N = N_vals[layer] ;
      M = M_vals[layer] ;
      K_S = K_S_vals[layer] ;
      THETA_R = THETA_R_vals[layer] ;
      THETA_S = THETA_S_vals[layer] ;
      STORATIVITY = STORATIVITY_vals[layer] ;
      if (h < 0) return 
	  K_S*(pow(1 - pow(-(ALPHA*h),M*N)/
	  pow(1 + pow(-(ALPHA*h),N),M),2)/
      pow(1 + pow(-(ALPHA*h),N),M/2.)) ;
      else return K_S;    
    }
    else if (h<0) {
      location = -int(h*100) ;
      value = (K_TABLE[layer][location+1] - K_TABLE[layer][location])*(-100*h-location)+K_TABLE[layer][location] ;
      return value ;
    }
    else return K_S_vals[layer] ;
  }
}

// dK/dh (van Genuchten).
double dKdh(double h, int layer)
{
  double value ;
  int location ;
  
  if (h > LOW_LIMIT && h < 0 && POLYNOMIALS_READY == true){
    return horner(POLYNOMIALS[layer][1], h, (5+NUM_OF_INSIDE_PTS));
  }
  else
  {
    if (CONSTITUTIVE_TABLES_READY == false || h < -15000.0)  {
      ALPHA = ALPHA_vals[layer] ;
      N = N_vals[layer] ;
      M = M_vals[layer] ;
      K_S = K_S_vals[layer] ;
      THETA_R = THETA_R_vals[layer] ;
      THETA_S = THETA_S_vals[layer] ;
      STORATIVITY = STORATIVITY_vals[layer] ;
      if (h < 0) return 
	  K_S*((ALPHA*pow(-(ALPHA*h),-1 + N)*
	  pow(1 + pow(-(ALPHA*h),N),-1 - M/2.)*
	  pow(1 - pow(-(ALPHA*h),M*N)/
	  pow(1 + pow(-(ALPHA*h),N),M),2)*M*N)/2. + 
	  (2*(1 - pow(-(ALPHA*h),M*N)/
	  pow(1 + pow(-(ALPHA*h),N),M))*
	  (-(ALPHA*pow(-(ALPHA*h),-1 + N + M*N)*
	  pow(1 + pow(-(ALPHA*h),N),-1 - M)*M*N) + 
	  (ALPHA*pow(-(ALPHA*h),-1 + M*N)*M*N)/
	  pow(1 + pow(-(ALPHA*h),N),M)))/
	  pow(1 + pow(-(ALPHA*h),N),M/2.));
      else return 0;
    }
    else if (h < 0) {
      location = -int(h*100) ;
      value = (dKdh_TABLE[layer][location+1] - dKdh_TABLE[layer][location])*(-100*h-location)+dKdh_TABLE[layer][location] ;
      return value ;
    }
    else return 0 ;
  }
}

// ddK/dhh (van Genuchten).
double ddKdhh(double h, int layer)
{

  int location ;
  double value ;
  
//     if (h < -15001) {
//       printf ("strange values of solution \n") ;
//       exit(1) ;
//     }
  
  
  if (h > LOW_LIMIT && h < 0 && POLYNOMIALS_READY == true){
    return horner(POLYNOMIALS[layer][2], h, (4+NUM_OF_INSIDE_PTS));
  }
  else
  {
    if (CONSTITUTIVE_TABLES_READY == false || h < -15000.0)  {
      ALPHA = ALPHA_vals[layer] ;
      N = N_vals[layer] ;
      M = M_vals[layer] ;
      K_S = K_S_vals[layer] ;
      THETA_R = THETA_R_vals[layer] ;
      THETA_S = THETA_S_vals[layer] ;
      STORATIVITY = STORATIVITY_vals[layer] ;
      //consider (cuts off singularity of this strange function. but still consider revision)
//       if (h > -0.5 && h < 0) {
// 	h = -0.5 ;
//       }
//       else
//       {
      if (h  < 0 ) return 
	    K_S*( -(pow(ALPHA,2)*pow(-(ALPHA*h),-2 + N)*
	  pow(1 + pow(-(ALPHA*h),N),-1 - M/2.)*
	  pow(1 - pow(-(ALPHA*h),M*N)/
	  pow(1 + pow(-(ALPHA*h),N),M),2)*M*(-1 + N)*N)/2.\
	    - (pow(ALPHA,2)*pow(-(ALPHA*h),-2 + 2*N)*
	    pow(1 + pow(-(ALPHA*h),N),-2 - M/2.)*
	  pow(1 - pow(-(ALPHA*h),M*N)/
	      pow(1 + pow(-(ALPHA*h),N),M),2)*(-1 - M/2.)*M*
	  pow(N,2))/2. + 2*ALPHA*pow(-(ALPHA*h),-1 + N)*
	    pow(1 + pow(-(ALPHA*h),N),-1 - M/2.)*
	  (1 - pow(-(ALPHA*h),M*N)/
	    pow(1 + pow(-(ALPHA*h),N),M))*M*N*
	  (-(ALPHA*pow(-(ALPHA*h),-1 + N + M*N)*
	    pow(1 + pow(-(ALPHA*h),N),-1 - M)*M*N) + 
	  (ALPHA*pow(-(ALPHA*h),-1 + M*N)*M*N)/
	  pow(1 + pow(-(ALPHA*h),N),M)) + 
	  (2*pow(-(ALPHA*pow(-(ALPHA*h),-1 + N + M*N)*
	      pow(1 + pow(-(ALPHA*h),N),-1 - M)*M*N) + 
	    (ALPHA*pow(-(ALPHA*h),-1 + M*N)*M*N)/
	    pow(1 + pow(-(ALPHA*h),N),M),2))/
	    pow(1 + pow(-(ALPHA*h),N),M/2.) + 
	  (2*(1 - pow(-(ALPHA*h),M*N)/
	    pow(1 + pow(-(ALPHA*h),N),M))*
	    (pow(ALPHA,2)*pow(-(ALPHA*h),-2 + 2*N + M*N)*
	    pow(1 + pow(-(ALPHA*h),N),-2 - M)*(-1 - M)*M*
	    pow(N,2) + pow(ALPHA,2)*
	    pow(-(ALPHA*h),-2 + N + M*N)*
	    pow(1 + pow(-(ALPHA*h),N),-1 - M)*pow(M,2)*
	    pow(N,2) - (pow(ALPHA,2)*
	      pow(-(ALPHA*h),-2 + M*N)*M*N*(-1 + M*N))/
	    pow(1 + pow(-(ALPHA*h),N),M) + 
	    pow(ALPHA,2)*pow(-(ALPHA*h),-2 + N + M*N)*
	    pow(1 + pow(-(ALPHA*h),N),-1 - M)*M*N*
	    (-1 + N + M*N)))/pow(1 + pow(-(ALPHA*h),N),M/2.)) ;

      else return 0;
    }
    else if (h < 0) {
      location = -int(h*100) ;
  //     printf("%lf %i \n ", h, location) ;
  //     exit(1) ;
      value = (ddKdhh_TABLE[layer][location+1] - ddKdhh_TABLE[layer][location])*(-100*h-location)+ddKdhh_TABLE[layer][location] ;
      return value ;
    }
    else return 0 ;
  }	
}

// C (van Genuchten).
double C(double h, int layer)
{
  int location ;
  double value ;
  
  if (CONSTITUTIVE_TABLES_READY == false || h < -15000.0)  {
    ALPHA = ALPHA_vals[layer] ;
    N = N_vals[layer] ;
    M = M_vals[layer] ;
    K_S = K_S_vals[layer] ;
    THETA_R = THETA_R_vals[layer] ;
    THETA_S = THETA_S_vals[layer] ;
    STORATIVITY = STORATIVITY_vals[layer] ;
    if (h < 0) return 
    ALPHA*pow(-(ALPHA*h),-1 + N)*
      pow(1 + pow(-(ALPHA*h),N),-1 - M)*M*N*(THETA_S - THETA_R) + 
    (STORATIVITY*((THETA_S - THETA_R)/pow(1 + pow(-(ALPHA*h),N),M) + 
	  THETA_R))/THETA_S;
    else return STORATIVITY; 
  }
  else if (h<0) {
    location = -int(h*100) ;
    value = (C_TABLE[layer][location+1] - C_TABLE[layer][location])*(-100*h-location)+C_TABLE[layer][location] ;
    return value ;
  }
  else return STORATIVITY_vals[layer] ;
}

// dC/dh (van Genuchten).
double dCdh(double h, int layer)
{
  int location ;
  double value ;
    
  if (CONSTITUTIVE_TABLES_READY == false|| h < -15000.0)  {
    ALPHA = ALPHA_vals[layer] ;
    N = N_vals[layer] ;
    M = M_vals[layer] ;
    K_S = K_S_vals[layer] ;
    THETA_R = THETA_R_vals[layer] ;
    THETA_S = THETA_S_vals[layer] ;
    STORATIVITY = STORATIVITY_vals[layer] ;
    if (h < 0) return
      -(pow(ALPHA,2)*pow(-(ALPHA*h),-2 + N)*
	pow(1 + pow(-(ALPHA*h),N),-1 - M)*M*(-1 + N)*N*
	(THETA_S - THETA_R)) - pow(ALPHA,2)*pow(-(ALPHA*h),-2 + 2*N)*
      pow(1 + pow(-(ALPHA*h),N),-2 - M)*(-1 - M)*M*pow(N,2)*
      (THETA_S - THETA_R) + (ALPHA*pow(-(ALPHA*h),-1 + N)*
	pow(1 + pow(-(ALPHA*h),N),-1 - M)*M*N*STORATIVITY*
	(THETA_S - THETA_R))/THETA_S;
    else return 0; 
  }
  else if (h<0) {
    location = -int(h*100) ;
    value = (dCdh_TABLE[layer][location+1] - dCdh_TABLE[layer][location])*(-100*h-location)+dCdh_TABLE[layer][location] ;
    return value ;
  }
  else return 0;
      
}


int gem_full(double** A, double* b, double* X, int n){
  int i,j,k;
  
  double** aa;
  double dotproduct, tmp;
  aa = new double*[n];
  
  for (i=0; i<n; i++){
    aa[i] = new double[n+1];
  }
  
  for (i=0;i<n; i++){
    for (j=0; j<n; j++){
      aa[i][j] = A[i][j] ;
    }
    aa[i][n] = b[i];
  }
  
//      for (i=0; i<7; i++){
//     printf("%lf %lf %lf %lf %lf %lf %lf %lf \n", aa[i][0], aa[i][1], aa[i][2], aa[i][3], aa[i][4], aa[i][5], aa[i][6], aa[i][7]) ;
//    
//   }
  


  for (j=0; j<(n-1); j++){
    for (i=j+1; i<n; i++){
    tmp = aa[i][j]/aa[j][j];

      for (k=0; k<(n+1); k++){
// 	            printf("tmp je %lf \n", tmp);

	aa[i][k] = aa[i][k] - tmp*aa[j][k] ;
      }
    }
//       for (i=0; i<7; i++){
//     printf("%lf %lf %lf %lf %lf %lf %lf %lf \n", aa[i][0], aa[i][1], aa[i][2], aa[i][3], aa[i][4], aa[i][5], aa[i][6], aa[i][7]) ;
//    
//   }
  }
  

//    printf("---------------- \n");
//   for (i=0; i<7; i++){
//     printf("%lf %lf %lf %lf %lf %lf %lf %lf \n", aa[i][0], aa[i][1], aa[i][2], aa[i][3], aa[i][4], aa[i][5], aa[i][6], aa[i][7]) ;
//   }
//     printf("---------------- \n");
  for (i=n-1; i>-1; i--){
    dotproduct=0.0;
    for (j=i+1; j<n; j++){
      dotproduct = dotproduct + aa[i][j]*X[j] ;
    }
    X[i] = (aa[i][n]-dotproduct)/aa[i][i] ;
  }
    
//     printf("%lf %lf %lf %lf %lf %lf %lf \n", X[0], X[1], X[2], X[3], X[4], X[5], X[6]) ;
  
  delete []aa;
  return 0;
}

double dot_product(double *v1, double *v2, int n)
{
  double result = 0;
  for (int i=0; i<n; i++){
    result = result + v1[i]*v2[i] ;
  }
  return result ; 
}

int multiply_vct(double *vct1, double a, int n, double *vct2){
  for (int i=0; i<n; i++){
    vct2[i]=a*vct1[i] ;
  }
  return 0 ;
}

int sum_vct(double *vct1, double a1,  double *vct2, double a2, double*vct3, int n){
  for (int i=0; i<n; i++){
    vct3[i] = a1*vct1[i]+a2*vct2[i];
  }
  return 0;
}


  
