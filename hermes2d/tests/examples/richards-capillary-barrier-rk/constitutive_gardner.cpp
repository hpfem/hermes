// K (Gardner).
double K(double h)
{
  if (h < 0) return K_S*exp(ALPHA*h);
  else return K_S;    
}

// dK/dh (Gardner).
double dKdh(double h)
{
  if (h < 0) return K_S*ALPHA*exp(ALPHA*h);
  else return 0;
}

// ddK/dhh (Gardner).
double ddKdhh(double h)
{
  if (h < 0) return K_S*ALPHA*ALPHA*exp(ALPHA*h);
  else return 0;
}

// C (Gardner).
double C(double h)
{
  if (h < 0) return ALPHA*(THETA_S - THETA_R)*exp(ALPHA*h);
  else return ALPHA*(THETA_S - THETA_R); 
//   else return STORATIVITY; 
}

// dC/dh (Gardner).
double dCdh(double h)
{
  if (h < 0) return ALPHA*(THETA_S - THETA_R)*ALPHA*exp(ALPHA*h);
  else return 0;    
}

