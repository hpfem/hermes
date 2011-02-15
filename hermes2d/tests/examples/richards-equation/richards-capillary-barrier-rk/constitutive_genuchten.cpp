// This function uses Horner scheme to efficiently evaluate values of 
// polynomials in approximating K(h) function close to full saturation.
double horner(double *pol, double x, int n){
  double px=0.0;
  for (int i=0; i<n; i++){
    px = px*x + pol[n-1-i];
  }
  return px;
}

// K (van Genuchten).
double K(double h, int layer)
{
  double value ;
  int location ;
  if (h > LOW_LIMIT && h < 0 && POLYNOMIALS_READY && CONSTITUTIVE_TABLE_METHOD == 1){
    return horner(POLYNOMIALS[layer][0], h, (6+NUM_OF_INSIDE_PTS));
  }
  else
  {
    if (CONSTITUTIVE_TABLE_METHOD == 0 || !CONSTITUTIVE_TABLES_READY || h < TABLE_LIMIT)  {
      ALPHA = ALPHA_vals[layer];
      N = N_vals[layer];
      M = M_vals[layer];
      K_S = K_S_vals[layer];
      THETA_R = THETA_R_vals[layer];
      THETA_S = THETA_S_vals[layer];
      STORATIVITY = STORATIVITY_vals[layer];
      if (h < 0) return 
	  K_S*(pow(1 - pow(-(ALPHA*h), M*N)/
	  pow(1 + pow(-(ALPHA*h), N), M),2)/
      pow(1 + pow(-(ALPHA*h), N), M/2.));
      else return K_S;    
    }
    else if (h<0 && CONSTITUTIVE_TABLE_METHOD == 1) {
      location = -int(h/TABLE_PRECISION) ;
      value = (K_TABLE[layer][location+1] - K_TABLE[layer][location])*(-h/TABLE_PRECISION-location)+K_TABLE[layer][location];
      return value;
    }
    else if (h<0 && CONSTITUTIVE_TABLE_METHOD ==2) {
      location = POL_SEARCH_HELP[int(-h)]; 
      return horner(K_POLS[location][layer][0], h, 6);
    }
    else return K_S_vals[layer];
  }
}

// dK/dh (van Genuchten).
double dKdh(double h, int layer)
{
  double value ;
  int location ;
  
  if (h > LOW_LIMIT && h < 0 && POLYNOMIALS_READY && CONSTITUTIVE_TABLE_METHOD == 1){
    return horner(POLYNOMIALS[layer][1], h, (5+NUM_OF_INSIDE_PTS));
  }
  else
  {
    if (CONSTITUTIVE_TABLE_METHOD == 0 || !CONSTITUTIVE_TABLES_READY || h < TABLE_LIMIT)  {
      ALPHA = ALPHA_vals[layer];
      N = N_vals[layer];
      M = M_vals[layer];
      K_S = K_S_vals[layer];
      THETA_R = THETA_R_vals[layer];
      THETA_S = THETA_S_vals[layer];
      STORATIVITY = STORATIVITY_vals[layer];
      if (h < 0) return 
	  K_S*((ALPHA*pow(-(ALPHA*h), -1 + N)*
	  pow(1 + pow(-(ALPHA*h), N), -1 - M/2.)*
	  pow(1 - pow(-(ALPHA*h), M*N)/
	  pow(1 + pow(-(ALPHA*h), N), M), 2)*M*N)/2. + 
	  (2*(1 - pow(-(ALPHA*h), M*N)/
	  pow(1 + pow(-(ALPHA*h), N), M))*
	  (-(ALPHA*pow(-(ALPHA*h), -1 + N + M*N)*
	  pow(1 + pow(-(ALPHA*h), N), -1 - M)*M*N) + 
	  (ALPHA*pow(-(ALPHA*h), -1 + M*N)*M*N)/
	  pow(1 + pow(-(ALPHA*h), N), M)))/
	  pow(1 + pow(-(ALPHA*h), N), M/2.));
      else return 0;
    }
    else if (h < 0 && CONSTITUTIVE_TABLE_METHOD == 1) {
      location = -int(h/TABLE_PRECISION);
      value = (dKdh_TABLE[layer][location+1] - dKdh_TABLE[layer][location])*(-h/TABLE_PRECISION-location)+dKdh_TABLE[layer][location] ;
      return value ;
    }
    else if (h<0 && CONSTITUTIVE_TABLE_METHOD ==2) {
      location = POL_SEARCH_HELP[int(-h)]; 
      return horner(K_POLS[location][layer][1], h, 5);
    }
    else return 0 ;
  }
}

// ddK/dhh (van Genuchten).
double ddKdhh(double h, int layer)
{

  int location ;
  double value ;
    
  
  if (h > LOW_LIMIT && h < 0 && POLYNOMIALS_READY && CONSTITUTIVE_TABLE_METHOD == 1){
    return horner(POLYNOMIALS[layer][2], h, (4+NUM_OF_INSIDE_PTS));
  }
  else
  {
    if (CONSTITUTIVE_TABLE_METHOD == 0 || !CONSTITUTIVE_TABLES_READY || h < TABLE_LIMIT) {
      ALPHA = ALPHA_vals[layer];
      N = N_vals[layer];
      M = M_vals[layer];
      K_S = K_S_vals[layer];
      THETA_R = THETA_R_vals[layer];
      THETA_S = THETA_S_vals[layer];
      STORATIVITY = STORATIVITY_vals[layer];

      if (h  < 0 ) return 
	    K_S*( -(pow(ALPHA, 2)*pow(-(ALPHA*h), -2 + N)*
	  pow(1 + pow(-(ALPHA*h), N), -1 - M/2.)*
	  pow(1 - pow(-(ALPHA*h), M*N)/
	  pow(1 + pow(-(ALPHA*h), N), M), 2)*M*(-1 + N)*N)/2.\
	    - (pow(ALPHA, 2)*pow(-(ALPHA*h), -2 + 2*N)*
	    pow(1 + pow(-(ALPHA*h), N), -2 - M/2.)*
	  pow(1 - pow(-(ALPHA*h), M*N)/
	      pow(1 + pow(-(ALPHA*h), N), M), 2)*(-1 - M/2.)*M*
	  pow(N, 2))/2. + 2*ALPHA*pow(-(ALPHA*h), -1 + N)*
	    pow(1 + pow(-(ALPHA*h), N), -1 - M/2.)*
	  (1 - pow(-(ALPHA*h), M*N)/
	    pow(1 + pow(-(ALPHA*h), N), M))*M*N*
	  (-(ALPHA*pow(-(ALPHA*h), -1 + N + M*N)*
	    pow(1 + pow(-(ALPHA*h), N), -1 - M)*M*N) + 
	  (ALPHA*pow(-(ALPHA*h), -1 + M*N)*M*N)/
	  pow(1 + pow(-(ALPHA*h), N), M)) + 
	  (2*pow(-(ALPHA*pow(-(ALPHA*h), -1 + N + M*N)*
	      pow(1 + pow(-(ALPHA*h), N), -1 - M)*M*N) + 
	    (ALPHA*pow(-(ALPHA*h), -1 + M*N)*M*N)/
	    pow(1 + pow(-(ALPHA*h), N), M), 2))/
	    pow(1 + pow(-(ALPHA*h), N), M/2.) + 
	  (2*(1 - pow(-(ALPHA*h), M*N)/
	    pow(1 + pow(-(ALPHA*h), N), M))*
	    (pow(ALPHA, 2)*pow(-(ALPHA*h), -2 + 2*N + M*N)*
	    pow(1 + pow(-(ALPHA*h), N), -2 - M)*(-1 - M)*M*
	    pow(N, 2) + pow(ALPHA, 2)*
	    pow(-(ALPHA*h), -2 + N + M*N)*
	    pow(1 + pow(-(ALPHA*h), N), -1 - M)*pow(M, 2)*
	    pow(N, 2) - (pow(ALPHA, 2)*
	      pow(-(ALPHA*h), -2 + M*N)*M*N*(-1 + M*N))/
	    pow(1 + pow(-(ALPHA*h), N), M) + 
	    pow(ALPHA, 2)*pow(-(ALPHA*h), -2 + N + M*N)*
	    pow(1 + pow(-(ALPHA*h), N), -1 - M)*M*N*
	    (-1 + N + M*N)))/pow(1 + pow(-(ALPHA*h), N), M/2.));

      else return 0;
    }
    else if (h < 0 && CONSTITUTIVE_TABLE_METHOD == 1) {
      location = -int(h/TABLE_PRECISION) ;
      value = (ddKdhh_TABLE[layer][location+1] - 
              ddKdhh_TABLE[layer][location])*(-h/TABLE_PRECISION-location)+ddKdhh_TABLE[layer][location];
      return value ;
    }
    else if (h<0 && CONSTITUTIVE_TABLE_METHOD == 2) {
      location = POL_SEARCH_HELP[int(-h)]; 
      return horner(K_POLS[location][layer][2], h, 4);
    }
    else return 0;
  }	
}

// C (van Genuchten).
double C(double h, int layer)
{
  int location ;
  double value ;
  
  if (CONSTITUTIVE_TABLE_METHOD == 0 || !CONSTITUTIVE_TABLES_READY || h < TABLE_LIMIT)  {
    ALPHA = ALPHA_vals[layer];
    N = N_vals[layer];
    M = M_vals[layer];
    K_S = K_S_vals[layer];
    THETA_R = THETA_R_vals[layer];
    THETA_S = THETA_S_vals[layer];
    STORATIVITY = STORATIVITY_vals[layer];
    if (h < 0) return 
    ALPHA*pow(-(ALPHA*h), -1 + N)*
      pow(1 + pow(-(ALPHA*h), N), -1 - M)*M*N*(THETA_S - THETA_R) + 
    (STORATIVITY*((THETA_S - THETA_R)/pow(1 + pow(-(ALPHA*h), N), M) + 
	  THETA_R))/THETA_S;
    else return STORATIVITY; 
  }
  else if (h<0 && CONSTITUTIVE_TABLE_METHOD == 1) {
    location = -int(h/TABLE_PRECISION);
    value = (C_TABLE[layer][location+1] - C_TABLE[layer][location])*(-h/TABLE_PRECISION-location)+C_TABLE[layer][location];
    return value;
  }
  else if (h<0 && CONSTITUTIVE_TABLE_METHOD == 2) {
      location = POL_SEARCH_HELP[int(-h)]; 
      return horner(C_POLS[location][layer][0], h, 4);
  }
  else return STORATIVITY_vals[layer];
}

// dC/dh (van Genuchten).
double dCdh(double h, int layer)
{
  int location ;
  double value ;
    
  if (CONSTITUTIVE_TABLE_METHOD == 0 || !CONSTITUTIVE_TABLES_READY || h < TABLE_LIMIT) {
    ALPHA = ALPHA_vals[layer];
    N = N_vals[layer];
    M = M_vals[layer];
    K_S = K_S_vals[layer];
    THETA_R = THETA_R_vals[layer];
    THETA_S = THETA_S_vals[layer];
    STORATIVITY = STORATIVITY_vals[layer];
    if (h < 0) return
      -(pow(ALPHA, 2)*pow(-(ALPHA*h), -2 + N)*
	pow(1 + pow(-(ALPHA*h), N), -1 - M)*M*(-1 + N)*N*
	(THETA_S - THETA_R)) - pow(ALPHA, 2)*pow(-(ALPHA*h), -2 + 2*N)*
      pow(1 + pow(-(ALPHA*h), N), -2 - M)*(-1 - M)*M*pow(N, 2)*
      (THETA_S - THETA_R) + (ALPHA*pow(-(ALPHA*h), -1 + N)*
	pow(1 + pow(-(ALPHA*h), N), -1 - M)*M*N*STORATIVITY*
	(THETA_S - THETA_R))/THETA_S;
    else return 0; 
  }
  else if (h<0 && CONSTITUTIVE_TABLE_METHOD == 1) {
    location = -int(h/TABLE_PRECISION) ;
    value = (dCdh_TABLE[layer][location+1] - 
            dCdh_TABLE[layer][location])*(-h/TABLE_PRECISION-location)+dCdh_TABLE[layer][location];
    return value ;
  }
  else if (h<0 && CONSTITUTIVE_TABLE_METHOD == 2) {
      location = POL_SEARCH_HELP[int(-h)]; 
      return horner(C_POLS[location][layer][1], h, 3);
  }
  else return 0;
      
}







  
