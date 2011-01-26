template<typename real, typename scalar>
    scalar bilinear_form_left(int n, double *wt, Func<real> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e, ExtData<real> *data) 
    {
     scalar result = 0;
     for (int i=0; i < n; i++) {
       double x = e->x[i];
       double y=e->y[i];
       double z=e->z[i];
       result += (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i] + u->dz[i]*v->dz[i] + V(x,y,z)*u->val[i]*v->val[i])*wt[i];
     }
     return result;
    }

Ord bilinear_form_left_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return u->val[0] * v->val[0] * e->x[0] * e->x[0] * e->x[0]; // returning the sum of the degrees of the basis
                                                    // and test function plus two
}

    
template<typename real, typename scalar>
    scalar bilinear_form_right(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e, ExtData<scalar> *data) 
    {
      return int_u_v<real, scalar>(n, wt, u, v, e);
    }
