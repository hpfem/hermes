// Exact solution and its derivatives.
double exact(double x, double y, double &dx, double &dy)
{
	dx = (1- 2*x) * y * (1 - y);
	dy = (1- 2*y) * x * (1 - x);
	return  x * y * (1-x) * (1-y);
}

template<typename Real>
Real u(Real x, Real y)  {  return (x - x*x) * (y - y*y);  }
template<typename Real>
Real dudx(Real x, Real y)  {  return (1- 2*x) * y * (1 - y);  }
template<typename Real>
Real dudy(Real x, Real y)  {  return (1- 2*y) * x * (1 - x);  }

template<typename Real>
Real dudxx(Real x, Real y)  {  return -2.0 * (y-y*y);  }
template<typename Real>
Real dudyy(Real x, Real y)  {  return -2.0 * (x-x*x);  }
template<typename Real>
Real dudxy(Real x, Real y)  {  return (1- 2*y) * (1 - 2*x);  }

template<typename Real>
Real k(Real x, Real y)  {  return 1.0 / sqrt(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)));  }
template<typename Real>
Real kx(Real x, Real y)  {  return -0.5 * pow(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)), -1.5) *
                 (2.0 * dudx(x,y) * dudxx(x,y) + 2.0 * dudy(x,y) * dudxy(x,y));  }
template<typename Real>
Real ky(Real x, Real y)  {  return -0.5 * pow(1.0 + sqr(dudx(x,y)) + sqr(dudy(x,y)), -1.5) *
                 (2.0 * dudx(x,y) * dudxy(x,y) + 2.0 * dudy(x,y) * dudyy(x,y));  }

template<typename Real>
Real f(Real x, Real y)
{  return - kx(x,y) * dudx(x,y) - ky(x,y) * dudy(x,y) - k(x,y) * (dudxx(x,y) + dudyy(x,y)); }
