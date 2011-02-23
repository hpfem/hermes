//////////////////////////////////////////////////// CONTINUOUS APPROXIMATION ////////////////////////////////////////////////////
//
// References:
//
// R. Codina:
// Comparison of some finite element methods for solving the diffusion-convection-reaction equation.
// Comput. Meth. Appl. Mech. Engrg. 156 (1998) 185-210.
//
// F. Shakib, T.J.R. Hughes, Z. Johan: 
// A new finite element formulation for computational fluid dynamics: X. The compressible Euler and Navier-Stokes equations. 
// Comput. Methods Appl. Mech. Engrg. 89 (1991) 141-219.
//

// Bilinear form.
template<typename Real, typename Scalar>
Scalar cg_biform(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i=0; i < n; i++)
  {
    result += wt[i] * (EPSILON * (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i])
                               - (B1 * u->val[i] * v->dx[i] + B2 * u->val[i] * v->dy[i])
                      );
  }
  return result;
}

// Streamline upwind Petrov-Galerkin stabilization, stabilization parameter according to Codina (derived from nodal exactness).
template<typename Real, typename Scalar>
Scalar stabilization_biform_supg(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  Real h_e = e->diam;
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    double b_norm = sqrt(B1*B1 + B2*B2);
    Real Pe = b_norm * h_e / (2*EPSILON);
    Real alpha = 1 + 2./(exp(2*Pe) - 1) - 1./Pe;  // coth(Pe)-1/Pe
    Real tau = alpha*h_e / (2*b_norm);
    
    result += wt[i] * tau * (B1 * v->dx[i] + B2 * v->dy[i])
                          * (B1 * u->dx[i] + B2 * u->dy[i] - EPSILON * u->laplace[i]);
  }
  return result;
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h if you want to use second derivatives of shape functions in weak forms.");
#endif  
}

// Galerkin least-squares stabilization, stabilization parameter according to Codina.
template<typename Real, typename Scalar>
Scalar stabilization_biform_gls(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  Real h_e = e->diam;
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    double b_norm = sqrt(B1*B1 + B2*B2);
    Real Pe = b_norm * h_e / (2*EPSILON);
    double C1 = (P_INIT == 1) ? 1./3. : 1./9.;
    double C2 = (P_INIT == 1) ? 1.    : 1./2.;
    Real tau = std::min(C1*Pe, C2)*h_e / (2*b_norm);
    
    result += wt[i] * tau * (B1 * v->dx[i] + B2 * v->dy[i] - EPSILON * v->laplace[i])
                          * (B1 * u->dx[i] + B2 * u->dy[i] - EPSILON * u->laplace[i]);
  }
  return result;
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h if you want to use second derivatives of shape functions in weak forms.");
#endif
}

template<>
Ord stabilization_biform_gls(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(  (B1 * v->dx[0] + B2 * v->dy[0] - EPSILON * v->laplace[0])
             * (B1 * u->dx[0] + B2 * u->dy[0] - EPSILON * u->laplace[0]) );
}

// Subgrid scale stabilization, stabilization parameter from discrete max. principle, according to Codina.
//
template<typename Real, typename Scalar>
Scalar stabilization_biform_sgs(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  Real h_e = e->diam;
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    double b_norm = sqrt(B1*B1 + B2*B2);   
    Real tau = 1. / (4*EPSILON/sqr(h_e) + 2*b_norm/h_e);
    result += wt[i] * tau * (B1 * v->dx[i] + B2 * v->dy[i] + EPSILON * v->laplace[i])
                          * (B1 * u->dx[i] + B2 * u->dy[i] - EPSILON * u->laplace[i]);
  }
  return result;
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h if you want to use second derivatives of shape functions in weak forms.");
#endif
}

// Subgrid scale stabilization, stabilization parameter according to Shakib.
//
template<typename Real, typename Scalar>
Scalar stabilization_biform_sgs_alt(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
#ifdef H2D_SECOND_DERIVATIVES_ENABLED
  Real h_e = e->diam;
  Scalar result = 0;
  for (int i=0; i < n; i++) {
    double b_norm = sqrt(B1*B1 + B2*B2);   
    Real tau = 1. / sqrt(9*sqr(4*EPSILON/sqr(h_e)) + sqr(2*b_norm/h_e));
    result += wt[i] * tau * (B1 * v->dx[i] + B2 * v->dy[i] + EPSILON * v->laplace[i])
                          * (B1 * u->dx[i] + B2 * u->dy[i] - EPSILON * u->laplace[i]);
  }
  return result;
#else
  error("Define H2D_SECOND_DERIVATIVES_ENABLED in h2d_common.h if you want to use second derivatives of shape functions in weak forms.");
#endif
}



//////////////////////////////////////////////////// DISCONTINUOUS APPROXIMATION ////////////////////////////////////////////////////

// Flux definition.

template<typename Real>
inline Real calculate_a_dot_v(Real x, Real y, Real vx, Real vy) 
{
  return B1*vx + B2*vy;
}

template<typename Real, typename Scalar>
inline Scalar upwind_flux(Real u_cent, Real u_neib, Real a_dot_n)
{
  Scalar eff_val;
  if (a_dot_n > 0) eff_val = u_cent;
  else if (a_dot_n < 0) eff_val = u_neib;
  else eff_val = 0.5*(u_cent + u_neib);
  return a_dot_n * eff_val;
}

template<>
inline Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n)
{
  return a_dot_n * (u_cent + u_neib); 
}

// Weak forms.

#define JUMP(w)       ( w->get_val_central(i) - w->get_val_neighbor(i) )
#define AVG_GRAD(w)   ( 0.5 * ( (w->get_dx_central(i) + w->get_dx_neighbor(i))*e->nx[i] + \
                                (w->get_dy_central(i) + w->get_dy_neighbor(i))*e->ny[i] ) )

//--- ADVECTION

template<typename Real, typename Scalar>
Scalar dg_volumetric_biform_advection(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += -wt[i] * u->val[i] * calculate_a_dot_v<Real>(e->x[i], e->y[i], v->dx[i], v->dy[i]);
  return result;
}

template<typename Real, typename Scalar>
Scalar dg_interface_biform_advection(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;

  for (int i = 0; i < n; i++) {
    Real a_dot_n = calculate_a_dot_v<Real>(e->x[i], e->y[i], e->nx[i], e->ny[i]);
    result += wt[i] * upwind_flux<Real,Scalar>(u->get_val_central(i), u->get_val_neighbor(i), a_dot_n) * JUMP(v);
  }
  
  return result;
}

template<typename Real, typename Scalar>
Scalar dg_boundary_biform_advection(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++) {
    Real a_dot_n = calculate_a_dot_v<Real>(e->x[i], e->y[i], e->nx[i], e->ny[i]);
    result += wt[i] * upwind_flux<Real,Scalar>(u->val[i], 0, a_dot_n) * v->val[i];
  }
  
  return result;
}

template<typename Real, typename Scalar>
Scalar dg_boundary_liform_advection(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  
  for (int i = 0; i < n; i++) {
    Real x = e->x[i], y = e->y[i];
    Real a_dot_n = calculate_a_dot_v<Real>(x, y, e->nx[i], e->ny[i]);
    result += -wt[i] * upwind_flux<Real,Scalar>(0, essential_bc_values<Real,Scalar>(e->edge_marker, x, y), a_dot_n) * v->val[i];
  }
  
  return result;
}

//--- DIFFUSION AND PENALTIES
//
// SIPG (Symmetric Interior Penalty Galerkin): theta = -1, C_W nontrivial, optimal convergence rate
// IIPG (Incomplete Interior Penalty Galerkin): theta = 0, C_W nontrivial, sub-optimal convergence rate
// NIPG (Nonsymmetric Interior Penalty Galerkin): theta = 1, C_W arbitrary, sub-optimal convergence rate
// (I haven't searched the literature for a proof of these claims yet, it is from a mini-seminar lectured by Vit Dolejsi).
//
// Dirichlet b.c. are weakly imposed through the boundary biforms, no linear form is thus present 
// (there is neither a source term nor a Neumann boundary in this benchmark problem).
//
const int theta = 1; // Using NIPG.
const int C_W = 1;

template<typename Real, typename Scalar>
Scalar dg_volumetric_biform_diffusion(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  for (int i = 0; i < n; i++)
    result += wt[i] * EPSILON * ( u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i] );
  return result;
}

template<typename Real, typename Scalar>
Scalar dg_interface_biform_diffusion(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  //Real sigma = 2 * C_W / (e->diam + e->get_neighbor_diam());
  Real edge_len = 0.;
  for (int i = 0; i < n; i++)
    edge_len += wt[i];
  
  Real sigma = C_W * EPSILON / (0.5*edge_len);
  
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * EPSILON * ( -AVG_GRAD(u) * JUMP(v) + theta * AVG_GRAD(v) * JUMP(u) ); // diffusion
    result += wt[i] * sigma * JUMP(u) * JUMP(v);                          // interior discontinuity penalization
  }
  return result;
}

template<>
Ord dg_interface_biform_diffusion(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  int i = 0;
  return AVG_GRAD(u) * JUMP(v) + theta * AVG_GRAD(v) * JUMP(u) + JUMP(u) * JUMP(v);
}

template<typename Real, typename Scalar>
Scalar dg_boundary_biform_diffusion(int n, double *wt, Func<Real> *u_ext[], Func<Real> *u, Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  //Real sigma = C_W * EPSILON / e->diam;
  Real edge_len = 0.;
  for (int i = 0; i < n; i++)
    edge_len += wt[i];
  
  Real sigma = C_W * EPSILON / (0.5*edge_len);
  
  for (int i = 0; i < n; i++)
  {
    result += wt[i] * EPSILON * ( -( u->dx[i]*e->nx[i] + u->dy[i]*e->ny[i] ) * v->val[i]
                                  + theta * ( v->dx[i]*e->nx[i] + v->dy[i]*e->ny[i] ) * (u->val[i]) );   
    result += wt[i] * sigma * ( u->val[i] ) * v->val[i];                    
  }
  return result;
}

template<typename Real, typename Scalar>
Scalar dg_boundary_liform_diffusion(int n, double *wt, Func<Real> *u_ext[], Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext)
{
  Scalar result = 0;
  //Real sigma = C_W * EPSILON / e->diam;
  Real edge_len = 0.;
  for (int i = 0; i < n; i++)
    edge_len += wt[i];
  
  Real sigma = C_W * EPSILON / (0.5*edge_len);
  
  for (int i = 0; i < n; i++)
  {
    Scalar u_dir = essential_bc_values<Real,Scalar>(e->edge_marker, e->x[i], e->y[i]);
    result += wt[i] * EPSILON * theta * ( v->dx[i]*e->nx[i] + v->dy[i]*e->ny[i] ) * u_dir;
    result += wt[i] * sigma * u_dir * v->val[i];                        
  }
  return result;
}