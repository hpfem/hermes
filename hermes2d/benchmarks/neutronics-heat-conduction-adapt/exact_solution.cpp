scalar phi_exact(double x, double y, double& dx, double& dy)
{
  dx = CF*PHI_FTIME(x,y)*(sin(PI*x/LX)*sin(PI*y/LY)/LX*y/LY
		     + PI/LX*cos(PI*x/LX)*sin(PI*y/LY)*x/LX*y/LY);
  dy = CF*PHI_FTIME(x,y)*(sin(PI*x/LX)*sin(PI*y/LY)*x/LX/LY
		     + sin(PI*x/LX)*PI/LY*cos(PI*y/LY)*x/LX*y/LY);
  return CF*PHI_FTIME(x,y)*sin(PI*x/LX)*sin(PI*y/LY)*x/LX*y/LY;
}

scalar T_exact(double x, double y, double& dx, double& dy)
{
  dx = CT*T_FTIME(x,y)*PI/LX*cos(PI*x/LX)*sin(PI*y/LY);
  dy = CT*T_FTIME(x,y)*sin(PI*x/LX)*PI/LY*cos(PI*y/LY);
  return CT*T_FTIME(x,y)*sin(PI*x/LX)*sin(PI*y/LY);
}
