int P=4;
template<typename real, typename scalar>
scalar bilinear_form_left(int n, double *wt, Func<real> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e, ExtData<real> *ext) 
    {
     scalar result = 0;
     Func<scalar>* pot = ext->fn[0];
     for (int i=0; i < n; i++) {
       result += (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i] + u->dz[i]*v->dz[i] + pot->val[i]*u->val[i]*v->val[i])*wt[i];
     }
     return result;
    }


template<typename real, typename scalar>
scalar  bilinear_form_right(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, 
			    Func<real> *v, Geom<real> *e, ExtData<scalar> *ext)
{
  scalar result = 0;
   for (int i=0; i < n; i++) {
    result +=u->val[i]*v->val[i]* wt[i];
  }
  return result;
}
Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(2*P+2);
}


