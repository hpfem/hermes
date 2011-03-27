#include <hermes2d.h>
#include <cstdlib>

// Class that represents the semi-analytic solution on the rectangle [0,1]x[0,1].
// Contains coordinates of the Gauss-Kronrod quadrature nodes and the corresponding
// function values and quadrature weights, as read from the file supplied to the
// constructor.
class SemiAnalyticSolution
{  
  std::vector<long double> x; //  x coordinate.
  std::vector<long double> y; //  y coordinate.
  std::vector<long double> u; //  u(x,y).
  std::vector<long double> w; //  Gauss-Kronrod quadrature weights.
  
  unsigned long int n;  // Number of quadrature points (size of vectors x,y,u,w).
  
  bool NA;  // Indicates that the exact solution could not be read from the input file.
  
  public:
    SemiAnalyticSolution(std::string file);
    
    // The following two functions use the quadrature specified by the input file to
    // calculate integrals over the rectangle [0,1]x[0,1].
    //
    // Calculates L2-norm of the exact solution.
    double get_l2_norm();     
    //
    // Calculates L2-norm of the relative difference between
    // the exact solution and the one computed by Hermes.  
    double get_l2_rel_err(Solution *approx_sln);  
};


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


    







// Inner product between two 2D vectors.
template<typename Real>
inline Real dot2(Real x1, Real y1, Real x2, Real y2)
{
  return x1*x2 + y1*y2;
}

// Weak form of the source term.
template<typename Real, typename Scalar>
Scalar source_liform(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  return int_F_v<Real, Scalar>(n, wt, F<Real>, v, e);
}

//////////////////////////////////////////////////// CONTINUOUS APPROXIMATION ////////////////////////////////////////////////////

// Weak forms:

template<typename Real, typename Scalar>
Scalar cg_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++)
  {
    Real a = fn_a<Real>(e->x[i], e->y[i]);
    Real b = fn_b<Real>(e->x[i], e->y[i]);
    Real c = fn_c<Real>(e->x[i], e->y[i]);
    result += wt[i] * v->val[i] * (dot2<Real>(a, b, u->dx[i], u->dy[i]) + c * u->val[i]);
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar cg_boundary_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++)
  {
    Real a = fn_a<Real>(e->x[i], e->y[i]);
    Real b = fn_b<Real>(e->x[i], e->y[i]);
    Real beta_dot_n = dot2<Real>(a, b, e->nx[i], e->ny[i]);
    
    if (beta_dot_n < 0)
      result += -wt[i] * u->val[i] * v->val[i] * beta_dot_n;
  }
  return result;
}
template<>
Ord cg_boundary_biform(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  Ord result = 0;
  
  for (int i = 0; i < n; i++) 
  {
    Ord a = fn_a<Ord>(e->x[i], e->y[i]);
    Ord b = fn_b<Ord>(e->x[i], e->y[i]);
    Ord beta_dot_n = dot2<Ord>(a, b, e->nx[i], e->ny[i]);
    result += wt[i] * u->val[i] * beta_dot_n * v->val[i];
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar cg_boundary_liform(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++) 
  {
    Real x = e->x[i], y = e->y[i];
    Real a = fn_a<Real>(x, y);
    Real b = fn_b<Real>(x, y);
    Real beta_dot_n = dot2<Real>(a, b, e->nx[i], e->ny[i]);
    
    if (beta_dot_n < 0)    // inflow
    {
      Scalar g = essential_bc_values<Real, Scalar>(e->edge_marker, x, y);
      result += -wt[i] * beta_dot_n * g * v->val[i];
    }
  }
  
  return result;
}
template<>
Ord cg_boundary_liform(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  Ord x = e->x[0], y = e->y[0];
  Ord a = fn_a<Ord>(x, y);
  Ord b = fn_b<Ord>(x, y);
  Ord beta_dot_n = dot2<Ord>(a, b, e->nx[0], e->ny[0]);
  Ord g = essential_bc_values<Ord,Ord>(e->edge_marker, x, y);

  return beta_dot_n * g * v->val[0];
}

// Empirical streamline upwind Petrov-Galerkin stabilization.
template<typename Real, typename Scalar>
Scalar stabilization_biform_supg(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Real h_e = e->diam;
  Scalar result = 0;
  Real norm_a_sq = 0.;
  Real norm_b_sq = 0.;
  for (int i=0; i < n; i++) 
  {
    Real a = fn_a<Real>(e->x[i], e->y[i]);
    Real b = fn_b<Real>(e->x[i], e->y[i]);
    Real c = fn_c<Real>(e->x[i], e->y[i]);
    Real f = F<Real>(e->x[i], e->y[i]);
    
    Real R = dot2<Real>(a, b, u->dx[i], u->dy[i]) + c * u->val[i] - f;    
    result += wt[i] * dot2<Real>(a, b, v->dx[i], v->dy[i]) * R;
    norm_a_sq += 0.5 * wt[i] * sqr(a);
    norm_b_sq += 0.5 * wt[i] * sqr(b);
  }
  
  return result * sqr(h_e)/(4*(norm_a_sq + norm_b_sq));
}


//////////////////////////////////////////////////// DISCONTINUOUS APPROXIMATION ////////////////////////////////////////////////////

// Scalar average, vector jump.

#define AVG(w)      ( 0.5 * (w->get_val_central(i) + w->get_val_neighbor(i)) )

#define JUMP(w)     w->get_val_central(i)*e->nx[i] - w->get_val_neighbor(i)*e->nx[i],\
                    w->get_val_central(i)*e->ny[i] - w->get_val_neighbor(i)*e->ny[i] 

// Weak forms:

template<typename Real, typename Scalar>
Scalar dg_volumetric_biform(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
  {
    Real a = fn_a<Real>(e->x[i], e->y[i]);
    Real b = fn_b<Real>(e->x[i], e->y[i]);
    Real c = fn_c<Real>(e->x[i], e->y[i]);
    result += wt[i] * u->val[i] * ( c * v->val[i] - dot2<Real>(a, b, v->dx[i], v->dy[i]) );
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar dg_interface_biform(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  Real theta = 0.5;   // Stabilization parameter. Standard upwind scheme is obtained for theta = 0.5.
    
  for (int i = 0; i < n; i++) 
  {
    Real a = fn_a<Real>(e->x[i], e->y[i]);
    Real b = fn_b<Real>(e->x[i], e->y[i]);
    Real beta_dot_n = dot2<Real>(a, b, e->nx[i], e->ny[i]);
    result += wt[i] * AVG(u) * dot2<Real>(a, b, JUMP(v));
    result += wt[i] * theta * magn(beta_dot_n) * dot2<Real>(JUMP(u), JUMP(v));
  }
  
  return result;
}

template<typename Real, typename Scalar>
Scalar dg_boundary_biform(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++) 
  {
    Real a = fn_a<Real>(e->x[i], e->y[i]);
    Real b = fn_b<Real>(e->x[i], e->y[i]);
    Real beta_dot_n = dot2<Real>(a, b, e->nx[i], e->ny[i]);
    if (beta_dot_n >= 0)   // outflow
      result += wt[i] * u->val[i] * beta_dot_n * v->val[i];
  }
  
  return result;
}
template<>
Ord dg_boundary_biform(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  Ord result = 0;
  
  for (int i = 0; i < n; i++) 
  {
    Ord a = fn_a<Ord>(e->x[i], e->y[i]);
    Ord b = fn_b<Ord>(e->x[i], e->y[i]);
    Ord beta_dot_n = dot2<Ord>(a, b, e->nx[i], e->ny[i]);
    result += wt[i] * u->val[i] * beta_dot_n * v->val[i];
  }
  
  return result;
}

template<typename Real, typename Scalar>
Scalar dg_boundary_liform(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++) 
  {
    Real x = e->x[i], y = e->y[i];
    Real a = fn_a<Real>(x, y);
    Real b = fn_b<Real>(x, y);
    Real beta_dot_n = dot2<Real>(a, b, e->nx[i], e->ny[i]);
    
    if (beta_dot_n < 0)    // inflow
    {
      Scalar g = essential_bc_values<Real, Scalar>(e->edge_marker, x, y);
      result += -wt[i] * beta_dot_n * g * v->val[i];
    }
  }
  
  return result;
}
template<>
Ord dg_boundary_liform(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  Ord x = e->x[0], y = e->y[0];
  Ord a = fn_a<Ord>(x, y);
  Ord b = fn_b<Ord>(x, y);
  Ord beta_dot_n = dot2<Ord>(a, b, e->nx[0], e->ny[0]);
  Ord g = essential_bc_values<Ord,Ord>(e->edge_marker, x, y);
  
  return beta_dot_n * g * v->val[0];
}
