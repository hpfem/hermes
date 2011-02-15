// Integral errors/norms in the axisymmetric coordinate system
 
// function used to combine the contributions from solution components to the total error 
double error_total( double (*efn)(MeshFunction*, MeshFunction*, RefMap*, RefMap*), 
                    double (*nfn)(MeshFunction*, RefMap*), 
                    Hermes::vector<Solution*>& slns1, 
                    Hermes::vector<Solution*>& slns2 )
{
  double error = 0.0, norm = 0.0;

  std::vector<Solution*>::iterator it1, it2;
  for (it1=slns1.begin(), it2=slns2.begin(); it1 < slns1.end(); it1++, it2++) {
    assert(it2 < slns2.end());
    error += sqr(calc_abs_error(efn, *it1, *it2));
    if (nfn) norm += sqr(calc_norm(nfn, *it2));
  }
  
  /* Alternatively, without iterators:
  double error = 0.0, norm = 0.0;

  for (int i=0; i < slns1.size(); i++) {
    error += sqr(calc_abs_error(efn, slns1[i], slns2[i]));
    if (nfn) norm += sqr(calc_norm(nfn, slns2[i]));
  }
  */

  return (nfn ? sqrt(error/norm) : sqrt(error));
}

//// H1 space //////////////////////////////////////////////////////////////////////////////////////

// Function used to calculate error in H1 norm.
double error_fn_h1_axisym(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2*std::max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o);
  sln2->set_quad_order(o);

  scalar *uval, *vval, *dudx, *dudy, *dvdx, *dvdy;
  uval = sln1->get_fn_values();
  vval = sln2->get_fn_values();
  sln1->get_dx_dy_values(dudx, dudy);
  sln2->get_dx_dy_values(dvdx, dvdy);

  double* x = ru->get_phys_x(o);
  double result = 0.0;
  h1_integrate_expression(x[i] * (sqr(uval[i] - vval[i]) +
                                  sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i])) );
  return 2*M_PI*result;
}

// Function used to calculate H1 norm of the solution.
double norm_fn_h1_axisym(MeshFunction* sln, RefMap* ru)
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 * sln->get_fn_order() + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o);

  scalar *uval, *dudx, *dudy;
  uval = sln->get_fn_values();
  sln->get_dx_dy_values(dudx, dudy);

  double* x = ru->get_phys_x(o);
  double result = 0.0;
  h1_integrate_expression(x[i] * (sqr(uval[i]) + sqr(dudx[i]) + sqr(dudy[i])) );
  return 2*M_PI*result;
}

double h1_error_axisym(MeshFunction* sln1, MeshFunction* sln2)
{
  double error = calc_abs_error(error_fn_h1_axisym, sln1, sln2);
  double norm = calc_norm(norm_fn_h1_axisym, sln2);
  return error/norm;
}

double h1_norm_axisym(MeshFunction* sln)
{
  return calc_norm(norm_fn_h1_axisym, sln);
}


//// L2 space //////////////////////////////////////////////////////////////////////////////////////

// function used to calculate error in L2 norm
double error_fn_l2_axisym(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2*std::max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o, H2D_FN_VAL);
  sln2->set_quad_order(o, H2D_FN_VAL);

  scalar *uval, *vval;
  uval = sln1->get_fn_values();
  vval = sln2->get_fn_values();

	double* x = ru->get_phys_x(o);
  double result = 0.0;
  h1_integrate_expression(x[i]*sqr(uval[i] - vval[i]));
  return 2*M_PI*result;
}

// function used to calculate L2 norm of the solution
double norm_fn_l2_axisym(MeshFunction* sln, RefMap* ru)
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 *sln->get_fn_order() + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o, H2D_FN_VAL);

  scalar* uval = sln->get_fn_values();

	double* x = ru->get_phys_x(o);
  double result = 0.0;
  h1_integrate_expression(x[i]*sqr(uval[i]));
  return 2*M_PI*result;
}


double l2_error_axisym(MeshFunction* sln1, MeshFunction* sln2)
{
  double error = calc_abs_error(error_fn_l2_axisym, sln1, sln2);
  double norm = calc_norm(norm_fn_l2_axisym, sln2);
  return error/norm;
}

double l2_norm_axisym(MeshFunction* sln)
{
  return calc_norm(norm_fn_l2_axisym, sln);
}

