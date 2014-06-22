#include "coupling.h"

CouplingErrorFormVelocity::CouplingErrorFormVelocity(VelocityComponent component, double c_p)
  : MatrixFormVol<double>(component, 4), component(component), c_p(c_p) 
{
}

double CouplingErrorFormVelocity::value(int n, double *wt, Func<double> *u_ext[], Func<double> *sln_i, Func<double> *sln_j, Geom<double> *e, Func<double>* *ext) const
{
  double result = 0.;
  if(component == velX)
  {
    for (int point_i = 0; point_i < n; point_i++)
    {
      result += wt[point_i] * sln_i->val[point_i] * sln_j->dx[point_i] * c_p;
    }
  }
  else
  {
    for (int point_i = 0; point_i < n; point_i++)
    {
      result += wt[point_i] * sln_i->val[point_i] * sln_j->dy[point_i] * c_p;
    }
  }
  return result;
}

Ord CouplingErrorFormVelocity::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *sln_i, Func<Ord> *sln_j, Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  if(component == velX)
  {
    for (int point_i = 0; point_i < n; point_i++)
    {
      result += wt[point_i] * sln_i->val[point_i] * sln_j->dx[point_i] * c_p;
    }
  }
  else
  {
    for (int point_i = 0; point_i < n; point_i++)
    {
      result += wt[point_i] * sln_i->val[point_i] * sln_j->dy[point_i] * c_p;
    }
  }
  return result;
}

MatrixFormVol<double>* CouplingErrorFormVelocity::clone() const
{
  return new CouplingErrorFormVelocity(*this);
}

CouplingErrorFormTemperature::CouplingErrorFormTemperature(VelocityComponent component, double c_p)
  : MatrixFormVol<double>(4, component), component(component), c_p(c_p) 
{
}

double CouplingErrorFormTemperature::value(int n, double *wt, Func<double> *u_ext[], Func<double> *sln_j, Func<double> *sln_i, Geom<double> *e, Func<double>* *ext) const
{
  double result = 0.;
  if(component == velX)
  {
    for (int point_i = 0; point_i < n; point_i++)
    {
      result += wt[point_i] * sln_i->val[point_i] * sln_j->dx[point_i] * c_p;
    }
  }
  else
  {
    for (int point_i = 0; point_i < n; point_i++)
    {
      result += wt[point_i] * sln_i->val[point_i] * sln_j->dy[point_i] * c_p;
    }
  }
  return result;
}

Ord CouplingErrorFormTemperature::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *sln_j, 
                                   Func<Ord> *sln_i, Geom<Ord> *e, Func<Ord>* *ext) const
{
  Ord result = Ord(0);
  if(component == velX)
  {
    for (int point_i = 0; point_i < n; point_i++)
    {
      result += wt[point_i] * sln_i->val[point_i] * sln_j->dx[point_i] * c_p;
    }
  }
  else
  {
    for (int point_i = 0; point_i < n; point_i++)
    {
      result += wt[point_i] * sln_i->val[point_i] * sln_j->dy[point_i] * c_p;
    }
  }
  return result;
}

MatrixFormVol<double>* CouplingErrorFormTemperature::clone() const
{
  return new CouplingErrorFormTemperature(*this);
}