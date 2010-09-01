scalar phi_exact(double x, double y, double& dx, double& dy)
{
  dx = CF*(1+exp(rF*TIME))*sin(PI*x/LX)*sin(PI*y/LY)/LX*y/LY
       + CF*(1+exp(rF*TIME))*PI/LX*cos(PI*x/LX)*sin(PI*y/LY)*x/LX*y/LY;
  dy = CF*(1+exp(rF*TIME))*sin(PI*x/LX)*sin(PI*y/LY)*x/LX/LY
       + CF*(1+exp(rF*TIME))*sin(PI*x/LX)*PI/LY*cos(PI*y/LY)*x/LX*y/LY;
  return CF*(1+exp(rF*TIME))*sin(PI*x/LX)*sin(PI*y/LY)*x/LX*y/LY;
}

scalar T_exact(double x, double y, double& dx, double& dy)
{
  dx = CT*(1+tanh(rT*TIME))*PI/LX*(PI*x/LX)*sin(PI*y/LY);
  dy = CT*(1+tanh(rT*TIME))*sin(PI*x/LX)*PI/LY*sin(PI*y/LY);
  return CT*(1+tanh(rT*TIME))*sin(PI*x/LX)*sin(PI*y/LY);
}
