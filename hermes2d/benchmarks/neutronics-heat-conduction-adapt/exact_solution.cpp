scalar phi_exact(double x, double y, double& dx, double& dy)
{
  dx = CF*PHI_FTIME<double>()*(sin(M_PI*x/LX)*sin(M_PI*y/LY)/LX*y/LY
         + M_PI/LX*cos(M_PI*x/LX)*sin(M_PI*y/LY)*x/LX*y/LY);
  dy = CF*PHI_FTIME<double>()*(sin(M_PI*x/LX)*sin(M_PI*y/LY)*x/LX/LY
         + sin(M_PI*x/LX)*M_PI/LY*cos(M_PI*y/LY)*x/LX*y/LY);
  return CF*PHI_FTIME<double>()*sin(M_PI*x/LX)*sin(M_PI*y/LY)*x/LX*y/LY;
}

scalar T_exact(double x, double y, double& dx, double& dy)
{
  dx = CT*T_FTIME<double>()*M_PI/LX*cos(M_PI*x/LX)*sin(M_PI*y/LY);
  dy = CT*T_FTIME<double>()*sin(M_PI*x/LX)*M_PI/LY*cos(M_PI*y/LY);
  return CT*T_FTIME<double>()*sin(M_PI*x/LX)*sin(M_PI*y/LY);
}
