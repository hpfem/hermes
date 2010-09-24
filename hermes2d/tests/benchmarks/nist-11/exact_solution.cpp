static double fn(double x, double y)
{
  double theta = atan2(y,x);
  if (theta < 0) theta = theta + 2.*M_PI;
  double r = sqrt(x*x + y*y);

  double mu;
  if (theta <= M_PI/2.) {
    mu = cos((M_PI/2. - SIGMA)*TAU) * cos((theta - M_PI/2. + RHO)*TAU);
  }
  else {
    if (theta <= M_PI) {
      mu = cos(RHO*TAU) * cos((theta - M_PI + SIGMA)*TAU);
    }
    else {
      if (theta <= 3.*M_PI/2.) {
        mu = cos(SIGMA*TAU) * cos((theta - M_PI - RHO)*TAU);
      }
      else {
        mu = cos((M_PI/2. - RHO)*TAU) * cos((theta - 3.*M_PI/2. - SIGMA)*TAU);
      }
    }
  }

  return pow(r, TAU) * mu;
}

static double fndd(double x, double y, double& dx, double& dy)
{
  double theta = atan2(y,x);
  if (theta < 0) theta = theta + 2*M_PI;
  double r = sqrt(x*x + y*y);
  // x-derivative
  if (theta <= M_PI/2.) {
    dx = TAU*x*pow(r, (2.*(-1 + TAU/2.))) *
    cos((M_PI/2. - SIGMA)*TAU) *
    cos(TAU*(-M_PI/2. + RHO + theta)) +
    (TAU*y*pow(r, TAU)*cos((M_PI/2. - SIGMA)*TAU) *
    sin(TAU*(-M_PI/2. + RHO + theta))/(r*r));
  }
  else {
    if (theta <= M_PI) {
      dx = TAU*x * pow(r, (2.*(-1 + TAU/2.))) * cos(RHO*TAU) *
      cos(TAU*(-M_PI + SIGMA + theta)) +
      (TAU*y * pow(r, TAU) * cos(RHO*TAU) *
      sin(TAU*(-M_PI + SIGMA + theta))/(r*r));
    }
    else {
      if (theta <= 3.*M_PI/2.) {
        dx = TAU*x * pow(r, (2.*(-1 + TAU/2.))) * cos(SIGMA*TAU) *
        cos(TAU*(-M_PI - RHO + theta)) +
	(TAU*y * pow(r, TAU) * cos(SIGMA*TAU) *
	 sin(TAU*(-M_PI - RHO + theta))/(r*r));
      }
      else {
        dx = TAU*x* pow(r, (2*(-1 + TAU/2.))) *
        cos((M_PI/2. - RHO)*TAU) *
	cos(TAU*(-3.*M_PI/2. - SIGMA + theta)) +
	(TAU*y*pow(r, TAU) * cos((M_PI/2. - RHO)*TAU) *
	sin(TAU*(-3.*M_PI/2. - SIGMA + theta))/(r*r));
      }
    }
  }
  // y-derivative
  if (theta <= M_PI/2.) {
    dy = TAU*y * pow(r, (2*(-1 + TAU/2.))) *
      cos((M_PI/2. - SIGMA)*TAU) *
      cos(TAU*(-M_PI/2. + RHO + theta)) -
      (TAU * pow(r, TAU) * cos((M_PI/2. - SIGMA)*TAU) *
       sin(TAU*(-M_PI/2. + RHO + theta))*x/(r*r));
  }
  else {
    if (theta <= M_PI) {
      dy = TAU*y* pow(r, (2*(-1 + TAU/2.))) * cos(RHO*TAU) *
	cos(TAU*(-M_PI + SIGMA + theta)) -
        (TAU * pow(r, TAU) * cos(RHO*TAU) *
	 sin(TAU*(-M_PI + SIGMA + theta))*x/(r*r));
    }
    else {
      if (theta <= 3.*M_PI/2.) {
        dy = TAU*y * pow(r, (2*(-1 + TAU/2.))) * cos(SIGMA*TAU) *
	  cos(TAU*(-M_PI - RHO + theta)) -
	  (TAU * pow(r, TAU) * cos(SIGMA*TAU) *
	   sin(TAU*(-M_PI - RHO + theta))*x/(r*r));
      }
      else {
        dy = TAU*y * pow(r, (2*(-1 + TAU/2.))) *
        cos((M_PI/2. - RHO)*TAU) *
	  cos(TAU*(-3.*M_PI/2. - SIGMA + theta)) -
	  (TAU * pow(r, TAU) * cos((M_PI/2. - RHO)*TAU) *
	   sin(TAU*((-3.*M_PI)/2. - SIGMA + theta))*x/(r*r));
      }
    }
  }
  return fn(x,y);
}
