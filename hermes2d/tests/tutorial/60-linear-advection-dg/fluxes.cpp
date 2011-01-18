template<typename Real>
inline Real calculate_a_dot_v(Real x, Real y, Real vx, Real vy) 
{
  Real norm = std::max(1e-12, sqrt(sqr(x) + sqr(y)));
  return -y/norm*vx + x/norm*vy;
}

template<>
inline Ord calculate_a_dot_v(Ord x, Ord y, Ord vx, Ord vy) 
{
  return Ord(10);
}

template<typename Real, typename Scalar>
inline Scalar upwind_flux(Real u_cent, Real u_neib, Real a_dot_n)
{
  return a_dot_n * (a_dot_n >= 0 ? u_cent : u_neib); 
}

template<>
inline Ord upwind_flux(Ord u_cent, Ord u_neib, Ord a_dot_n)
{
  return a_dot_n * (u_cent + u_neib); 
}
