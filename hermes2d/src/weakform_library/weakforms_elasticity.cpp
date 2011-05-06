#include "../hermes2d.h"
#ifndef H2D_COMPLEX
namespace WeakFormsElasticity 
{
  DefaultJacobianElasticity_0_0::DefaultJacobianElasticity_0_0
    (unsigned int i, unsigned int j, double lambda, double mu)
    : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) 
  {
  }

  DefaultJacobianElasticity_0_0::DefaultJacobianElasticity_0_0
    (unsigned int i, unsigned int j, std::string area, double lambda, double mu)
    : WeakForm::MatrixFormVol(i, j, area, HERMES_SYM), lambda(lambda), mu(mu) 
  {
  }

  template<typename Real, typename Scalar>
  Scalar DefaultJacobianElasticity_0_0::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                    Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
  {
    return (lambda + 2*mu) * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
                        mu * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
  }

  scalar DefaultJacobianElasticity_0_0::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                                      Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
  {
    return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord DefaultJacobianElasticity_0_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                                 Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }


  DefaultJacobianElasticity_0_1::DefaultJacobianElasticity_0_1
    (unsigned int i, unsigned int j, double lambda, double mu)
    : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) 
  {
  }
  
  DefaultJacobianElasticity_0_1::DefaultJacobianElasticity_0_1
    (unsigned int i, unsigned int j, std::string area, double lambda, double mu)
    : WeakForm::MatrixFormVol(i, j, area, HERMES_SYM), lambda(lambda), mu(mu) 
  {
  }

  template<typename Real, typename Scalar>
  Scalar DefaultJacobianElasticity_0_1::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                      Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
  {
    return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u, v) +
                mu * int_dudx_dvdy<Real, Scalar>(n, wt, u, v);
  }

  scalar DefaultJacobianElasticity_0_1::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
  {
    return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord DefaultJacobianElasticity_0_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
          Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
  {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }


  DefaultResidualElasticity_0_0::DefaultResidualElasticity_0_0
    (unsigned int i, double lambda, double mu)
    : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu) 
  {
  }

  DefaultResidualElasticity_0_0::DefaultResidualElasticity_0_0
    (unsigned int i, std::string area, double lambda, double mu)
    : WeakForm::VectorFormVol(i,  area), lambda(lambda), mu(mu) 
  {
  }

  template<typename Real, typename Scalar>
  Scalar DefaultResidualElasticity_0_0::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                      Geom<Real> *e, ExtData<Scalar> *ext) const
  {
    return (2*mu + lambda) * int_dudx_dvdx<Real, Scalar>(n, wt, u_ext[0], v) +
                mu * int_dudy_dvdy<Real, Scalar>(n, wt, u_ext[0], v);
  }

  scalar DefaultResidualElasticity_0_0::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                        Geom<double> *e, ExtData<scalar> *ext) const
  {
    return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
  }

  Ord DefaultResidualElasticity_0_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                  Geom<Ord> *e, ExtData<Ord> *ext) const
  {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }
  

  DefaultResidualElasticity_0_1::DefaultResidualElasticity_0_1
    (unsigned int i, double lambda, double mu)
    : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu) 
  {
  }
    
  DefaultResidualElasticity_0_1::DefaultResidualElasticity_0_1
    (unsigned int i, std::string area, double lambda, double mu)
    : WeakForm::VectorFormVol(i,  area), lambda(lambda), mu(mu) 
  {
  }

  template<typename Real, typename Scalar>
  Scalar DefaultResidualElasticity_0_1::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                      Geom<Real> *e, ExtData<Scalar> *ext) const
  {
    return lambda * int_dudy_dvdx<Real, Scalar>(n, wt, u_ext[1], v) +
                mu * int_dudx_dvdy<Real, Scalar>(n, wt, u_ext[1], v);
  }

  scalar DefaultResidualElasticity_0_1::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                        Geom<double> *e, ExtData<scalar> *ext) const
  {
    return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
  }

  Ord DefaultResidualElasticity_0_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                  Geom<Ord> *e, ExtData<Ord> *ext) const
  {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }


  DefaultResidualElasticity_1_0::DefaultResidualElasticity_1_0
    (unsigned int i, double lambda, double mu)
    : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu) 
  {
  }
    
  DefaultResidualElasticity_1_0::DefaultResidualElasticity_1_0
    (unsigned int i, std::string area, double lambda, double mu)
    : WeakForm::VectorFormVol(i,  area), lambda(lambda), mu(mu) 
  {
  }

  template<typename Real, typename Scalar>
  Scalar DefaultResidualElasticity_1_0::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                      Geom<Real> *e, ExtData<Scalar> *ext) const
  {
    return mu * int_dudy_dvdx<Real, Scalar>(n, wt, u_ext[0], v) +
            lambda * int_dudx_dvdy<Real, Scalar>(n, wt, u_ext[0], v);
  }

  scalar DefaultResidualElasticity_1_0::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                        Geom<double> *e, ExtData<scalar> *ext) const
  {
    return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
  }

  Ord DefaultResidualElasticity_1_0::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                  Geom<Ord> *e, ExtData<Ord> *ext) const
  {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }
  
    
  DefaultResidualElasticity_1_1::DefaultResidualElasticity_1_1
    (unsigned int i, double lambda, double mu)
    : WeakForm::VectorFormVol(i), lambda(lambda), mu(mu) 
  {
  }
    
  DefaultResidualElasticity_1_1::DefaultResidualElasticity_1_1
    (unsigned int i, std::string area, double lambda, double mu)
    : WeakForm::VectorFormVol(i,  area), lambda(lambda), mu(mu) 
  {
  }

  template<typename Real, typename Scalar>
  Scalar DefaultResidualElasticity_1_1::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v, 
                      Geom<Real> *e, ExtData<Scalar> *ext) const
  {
    return (2*mu + lambda) * int_dudy_dvdy<Real, Scalar>(n, wt, u_ext[1], v) +
                mu * int_dudx_dvdx<Real, Scalar>(n, wt, u_ext[1], v);
  }

  scalar DefaultResidualElasticity_1_1::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v, 
                        Geom<double> *e, ExtData<scalar> *ext) const
  {
    return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
  }

  Ord DefaultResidualElasticity_1_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v, 
                  Geom<Ord> *e, ExtData<Ord> *ext) const
  {
      return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }

   
  DefaultJacobianElasticity_1_1::DefaultJacobianElasticity_1_1
    (unsigned int i, unsigned int j, double lambda, double mu)
    : WeakForm::MatrixFormVol(i, j, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) 
  {
  }
    
  DefaultJacobianElasticity_1_1::DefaultJacobianElasticity_1_1
    (unsigned int i, unsigned int j, std::string area, double lambda, double mu)
    : WeakForm::MatrixFormVol(i, j, area, HERMES_SYM), lambda(lambda), mu(mu) 
  {
  }

  template<typename Real, typename Scalar>
  Scalar DefaultJacobianElasticity_1_1::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                      Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const
  {
    return              mu * int_dudx_dvdx<Real, Scalar>(n, wt, u, v) +
            (lambda + 2*mu) * int_dudy_dvdy<Real, Scalar>(n, wt, u, v);
  }

  scalar DefaultJacobianElasticity_1_1::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
  {
    return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord DefaultJacobianElasticity_1_1::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) const
  {
      return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }


  DefaultJacobianElasticity_00_11::DefaultJacobianElasticity_00_11
    (Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, double lambda, double mu)
    : WeakForm::MultiComponentMatrixFormVol(coordinates, HERMES_ANY, HERMES_SYM), lambda(lambda), mu(mu) 
  {
  }
    
  DefaultJacobianElasticity_00_11::DefaultJacobianElasticity_00_11
    (Hermes::vector<std::pair<unsigned int, unsigned int> >coordinates, std::string area, double lambda, double mu)
    : WeakForm::MultiComponentMatrixFormVol(coordinates, area, HERMES_SYM), lambda(lambda), mu(mu) 
  {
  }

  template<typename Real, typename Scalar>
  void DefaultJacobianElasticity_00_11::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u, Func<Real> *v,
                    Geom<Real> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const 
  {
    Scalar result_0_0 = 0;
    Scalar result_1_1 = 0;
    for (int i = 0; i < n; i++) {
      result_0_0 += wt[i] * ((lambda + 2*mu) * u->dx[i] * v->dx[i] + mu * u->dy[i] * v->dy[i]);
      result_1_1 += wt[i] * (mu * u->dx[i] * v->dx[i] + (lambda + 2*mu) * u->dy[i] * v->dy[i]);
    }
    result.push_back(result_0_0);
    result.push_back(result_1_1);
  }

  void DefaultJacobianElasticity_00_11::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                      Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const 
  {
    matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext, result);
  }

  Ord DefaultJacobianElasticity_00_11::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Hermes::vector<Ord> result;
    matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext, result);

    // Choose the maximum order.
    Ord to_return = result[0];
    if(result[1] > to_return)
      to_return = result[1];
    return to_return;
  }


  DefaultResidualElasticity_00_11::DefaultResidualElasticity_00_11
    (Hermes::vector<unsigned int> coordinates, double lambda, double mu)
    : WeakForm::MultiComponentVectorFormVol(coordinates), lambda(lambda), mu(mu) 
  {
  }
    
  DefaultResidualElasticity_00_11::DefaultResidualElasticity_00_11
    (Hermes::vector<unsigned int> coordinates, std::string area, double lambda, double mu)
    : WeakForm::MultiComponentVectorFormVol(coordinates, area), lambda(lambda), mu(mu) 
  {
  }

  template<typename Real, typename Scalar>
  void DefaultResidualElasticity_00_11::vector_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *v,
                      Geom<Real> *e, ExtData<Scalar> *ext, Hermes::vector<Scalar>& result) const 
  {
    Scalar result_0 = 0;
    Scalar result_1 = 0;
    for (int i = 0; i < n; i++) {
      result_0 += wt[i] * ((lambda + 2*mu) * u_ext[0]->dx[i] * v->dx[i] + mu * u_ext[0]->dy[i] * v->dy[i]);
      result_1 += wt[i] * (mu * u_ext[1]->dx[i] * v->dx[i] + (lambda + 2*mu) * u_ext[1]->dy[i] * v->dy[i]);
    }
    result.push_back(result_0);
    result.push_back(result_1);
  }

  void DefaultResidualElasticity_00_11::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                      Geom<double> *e, ExtData<scalar> *ext, Hermes::vector<scalar>& result) const 
  {
    vector_form<double, scalar>(n, wt, u_ext, v, e, ext, result);
  }

  Ord DefaultResidualElasticity_00_11::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                      Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Hermes::vector<Ord> result;
    vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext, result);

    // Choose the maximum order.
    Ord to_return = result[0];
    if(result[1] > to_return)
      to_return = result[1];
    return to_return;
  }
};
#endif