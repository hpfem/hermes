#include "exact.h"
#include <cstdlib>

SemiAnalyticSolution::SemiAnalyticSolution(std::string file) : NA(false)
{
  std::string exact_solution_not_available;
  exact_solution_not_available.append("Exact solution could not be read from the supplied file:");
  exact_solution_not_available.append(file);
  
  std::ifstream ifs(file.c_str());
  
  if (ifs.fail())
  {
    warning(exact_solution_not_available.c_str());
    NA = true;
    return;
  }
  
  std::string ns, xs, ys, us, ws;
  
  std::getline(ifs, ns);
  n = strtoul(ns.c_str(), NULL, 0);
  x.reserve(n);
  y.reserve(n);
  u.reserve(n);
  w.reserve(n);
  
  while(!std::getline(ifs, xs, '\t').eof())
  {
    x.push_back(strtold(xs.c_str(), NULL));
    std::getline(ifs, ys, '\t');
    y.push_back(strtold(ys.c_str(), NULL));
    std::getline(ifs, us, '\t');
    u.push_back(strtold(us.c_str(), NULL));
    std::getline(ifs, ws);
    w.push_back(strtold(ws.c_str(), NULL));
  }
  
  if (x.size() != n || y.size() != n || u.size() != n || w.size() != n)
  {
    warning(exact_solution_not_available.c_str());
    NA = true;
  }
  
  ifs.close();
}

double SemiAnalyticSolution::get_l2_norm()
{
  if (NA) return -1;
  
  long double res = 0.0;
  
  for (unsigned int i = 0; i < n; i++)
    res += w[i] * u[i] * u[i];
  
  return sqrt(res);
}

double SemiAnalyticSolution::get_l2_rel_err(Solution *approx_sln)
{
  if (NA) return -1;
  
  long double res = 0.0, nrm = 0.0;
  
  for (unsigned int i = 0; i < n; i++)
  {
    nrm += w[i] * u[i] * u[i];
    res += w[i] * pow(u[i] - approx_sln->get_pt_value(x[i],y[i]),2);
  }
  
  return sqrt(res/nrm);
}