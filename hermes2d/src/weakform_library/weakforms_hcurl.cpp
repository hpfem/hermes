#include "../hermes2d.h"

namespace WeakFormsHcurl {

  DefaultLinearCurlCurl::DefaultLinearCurlCurl(int i, int j, scalar coeff, SymFlag sym)
    : WeakForm::MatrixFormVol(i, j, HERMES_ANY, sym), coeff(coeff) 
  {
  }
  DefaultLinearCurlCurl::DefaultLinearCurlCurl(int i, int j, std::string area, scalar coeff, SymFlag sym)
    : WeakForm::MatrixFormVol(i, j, area, sym), coeff(coeff) 
  {
  }

  template<typename Real, typename Scalar>
  Scalar DefaultLinearCurlCurl::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                      Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
  {
    return coeff * int_curl_e_curl_f<Real, Scalar>(n, wt, u, v);
  }

  scalar DefaultLinearCurlCurl::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                        Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
  {
    return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord DefaultLinearCurlCurl::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                  Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  WeakForm::MatrixFormVol* DefaultLinearCurlCurl::clone() 
  {
    return new DefaultLinearCurlCurl(*this);
  }
     
      
  DefaultLinearMass::DefaultLinearMass(int i, int j, scalar coeff, SymFlag sym)
    : WeakForm::MatrixFormVol(i, j, HERMES_ANY, sym), coeff(coeff) { }
  DefaultLinearMass::DefaultLinearMass(int i, int j, std::string area, scalar coeff, SymFlag sym)
    : WeakForm::MatrixFormVol(i, j, area, sym), coeff(coeff) { }

  template<typename Real, typename Scalar>
  Scalar DefaultLinearMass::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                      Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
    return coeff * int_e_f<Real, Scalar>(n, wt, u, v);
  }

  scalar DefaultLinearMass::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                        Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const {
    return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord DefaultLinearMass::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) const {
    return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  WeakForm::MatrixFormVol* DefaultLinearMass::clone() {
    return new DefaultLinearMass(*this);
  }
     
      
  DefaultMatrixFormSurf::DefaultMatrixFormSurf(int i, int j, scalar coeff)
    : WeakForm::MatrixFormSurf(i, j), coeff(coeff) { }
  DefaultMatrixFormSurf::DefaultMatrixFormSurf(int i, int j, std::string area, scalar coeff)
    : WeakForm::MatrixFormSurf(i, j, area), coeff(coeff) { }

  template<typename Real, typename Scalar>
  Scalar DefaultMatrixFormSurf::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                      Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const {
    return coeff * int_e_tau_f_tau<Real, Scalar>(n, wt, u, v, e);
  }

  scalar DefaultMatrixFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                        Geom<double> *e, ExtData<scalar> *ext) const {
    return matrix_form(n, wt, u_ext, u, v, e, ext);
  }

  Ord DefaultMatrixFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                  Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const {
    return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  WeakForm::MatrixFormSurf* DefaultMatrixFormSurf::clone() {
    return new DefaultMatrixFormSurf(*this);
  }

      
  DefaultVectorFormConst::DefaultVectorFormConst(int i, scalar coeff0, scalar coeff1)
    : WeakForm::VectorFormVol(i), coeff0(coeff0), coeff1(coeff1) { }
  DefaultVectorFormConst::DefaultVectorFormConst(int i, std::string area, scalar coeff0, scalar coeff1)
      : WeakForm::VectorFormVol(i, area), coeff0(coeff0), coeff1(coeff1) { }

  scalar DefaultVectorFormConst::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                        Geom<double> *e, ExtData<scalar> *ext) const {
    scalar int_v0 = 0, int_v1 = 0;
    for (int i = 0; i < n; i++) int_v0 += wt[i] * v->val0[i];
    for (int i = 0; i < n; i++) int_v1 += wt[i] * v->val1[i];
    return coeff0 * int_v0 + coeff1 * int_v1;
  }

  Ord DefaultVectorFormConst::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
          Geom<Ord> *e, ExtData<Ord> *ext) const {
    Ord int_v0 = 0, int_v1 = 0;
    for (int i = 0; i < n; i++) int_v0 += wt[i] * v->val0[i];
    for (int i = 0; i < n; i++) int_v1 += wt[i] * v->val1[i];
    return coeff0 * int_v0 + coeff1 * int_v1;
  }

  WeakForm::VectorFormVol* DefaultVectorFormConst::clone() {
    return new DefaultVectorFormConst(*this);
  }
};