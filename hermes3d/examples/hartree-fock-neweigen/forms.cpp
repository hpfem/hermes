template<typename real, typename scalar>
scalar bilinear_form_laplace(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, 
                     Func<real> *v, Geom<real> *e, ExtData<scalar> *ext)
{
  scalar result = 0;
  for (int i=0; i < n; i++) {
    result += (u->dx[i]*v->dx[i]+  u->dy[i]*v->dy[i]+ u->dz[i]*v->dz[i])*wt[i];
  }
  return result;

}



Ord bilinear_form_ord1(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return u->val[0] * v->val[0] ; // returning the sum of the degrees of the basis and the test function 
}



template<typename real, typename scalar>
scalar bilinear_form_left(int n, double *wt, Func<real> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e, ExtData<real> *ext) 
    {
     scalar result = 0;
     Func<scalar>* pot = ext->fn[0];
     Func<scalar>* wfun = ext->fn[1];
     Func<scalar>* coul_pot = ext->fn[2];
     for (int i=0; i < n; i++) {
       result += (u->dx[i]*v->dx[i] + u->dy[i]*v->dy[i] + u->dz[i]*v->dz[i] + (pot->val[i]+coul_pot->val[i])*u->val[i]*v->val[i])*wfun->val[i]*wt[i];
     }
     return result;
    }


template<typename real, typename scalar>
scalar bilinear_form_coul_pot(int n, double *wt, Func<real> *u_ext[], Func<real> *u, Func<real> *v, Geom<real> *e, ExtData<real> *ext) 
    {
     scalar result = 0;
     Func<scalar>* wfun  = ext->fn[0];
     Func<scalar>* coul_pot = ext->fn[1];
     for (int i=0; i < n; i++) {
       result += coul_pot->val[i]*u->val[i]*v->val[i]*wfun->val[i]*wt[i];
     }
     return result;
    }


    
template<typename real, typename scalar>
scalar  bilinear_form_right(int n, double *wt, Func<scalar> *u_ext[], Func<real> *u, 
		       Func<real> *v, Geom<real> *e, ExtData<scalar> *ext)
{
  scalar result = 0;
  Func<scalar>* wfun = ext->fn[0];
  //info("n=%d",n);
   for (int i=0; i < n; i++) {
    result +=wfun->val[i]*u->val[i]*v->val[i]* wt[i];
  }
  return result;
}
Ord bilinear_form_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                      Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(2*P);
}


template<typename real, typename scalar>
scalar  linear_form_poisson(int n, double *wt, Func<scalar> *u_ext[], Func<real> *v, Geom<real> *e, ExtData<scalar> *data)
{
  scalar result = 0;
  Func<scalar>* phi = data->fn[0];
  Func<scalar>* wfun = data->fn[1];
  for (int i=0; i < n; i++) {
    result +=8*PI*wfun->val[i]*phi->val[i]*phi->val[i]*v->val[i]* wt[i];
    // 8*pi required on rhs of poisson equation ( 4PI) and factor 2 for rydberg units
  }
  return result;
}

Ord linear_form_poisson_ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext)
{
  return Ord(2*P);
}
