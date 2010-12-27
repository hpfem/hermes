// Inner product between two 2D vectors.
template<typename Real>
inline Real dot2(Real x1, Real y1, Real x2, Real y2)
{
  return x1*x2 + y1*y2;
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
    result += wt[i] * u->val[i] * (c * v->val[i] - dot2<Real>(a, b, v->dx[i], v->dy[i]));
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar cg_biform_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++)
  {
    Real a = fn_a<Real>(e->x[i], e->y[i]);
    Real b = fn_b<Real>(e->x[i], e->y[i]);
    Real beta_dot_n = dot2<Real>(a, b, e->nx[i], e->ny[i]);
    result += wt[i] * u->val[i] * v->val[i] * beta_dot_n;
  }
  return result;
}

// Empirical streamline upwind Petrov-Galerkin stabilization.
template<typename Real, typename Scalar>
Scalar stabilization_biform_supg(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  //Real h_e = e->diam;
  double h_e = MIN_DIAM;
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
    norm_a_sq += wt[i] * sqr(a);
    norm_b_sq += wt[i] * sqr(b);  
   }
  
  return result * sqr(h_e)/sqrt(norm_a_sq + norm_b_sq);
}
//////////////////////////////////////////////////// DISCONTINUOUS APPROXIMATION ////////////////////////////////////////////////////

// Scalar average, vector jump.

#define AVG(w)      ( 0.5 * (w->get_val_central(i) + w->get_val_neighbor(i)) )

#define JUMP(w)     w->get_val_central(i)*e->nx[i] - w->get_val_neighbor(i)*e->nx[i],\
                    w->get_val_central(i)*e->ny[i] - w->get_val_neighbor(i)*e->ny[i] 

// Stabilization of the advection term according to ref. [1].
template<typename Real>
inline Real dg_advection_stabilizer(Real a, Real b, Geom<Real>* surf, int pt)
{
  double theta = 1.0;
  switch (method)
  {
    case DG_UPWIND:
      return 0.5 * abs(dot2<Real>(a, b, surf->nx[pt], surf->ny[pt]));
      break;
    case DG_JUMP_STAB:
      return theta * abs(dot2<Real>(a, b, surf->nx[pt], surf->ny[pt])); 
      break;
  }
}

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
  
  for (int i = 0; i < n; i++) 
  {
    Real a = fn_a<Real>(e->x[i], e->y[i]);
    Real b = fn_b<Real>(e->x[i], e->y[i]);
    result += wt[i] * AVG(u) * dot2<Real>(a, b, JUMP(v));
    result += wt[i] * dg_advection_stabilizer<Real>(a, b, e, i) * dot2<Real>(JUMP(u), JUMP(v));
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