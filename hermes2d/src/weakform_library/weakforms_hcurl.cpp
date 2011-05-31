#include "../hermes2d.h"

namespace WeakFormsHcurl 
{
  DefaultMatrixFormVol<Scalar>::DefaultMatrixFormVol<Scalar>
    (int i, int j, std::string area, Scalar const_coeff, 
    DefaultFunction<Scalar>* f_coeff, SymFlag sym, 
    GeomType gt)
    : MatrixFormVol<Scalar>(i, j, area, sym), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultMatrixFormVol<Scalar>::DefaultMatrixFormVol<Scalar>
    (int i, int j, Hermes::vector<std::string> areas,Scalar const_coeff, 
    DefaultFunction<Scalar>* f_coeff, SymFlag sym, GeomType gt)
    : MatrixFormVol<Scalar>(i, j, areas, sym), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultMatrixFormVol<Scalar>::~DefaultMatrixFormVol<Scalar>() 
  {
    if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
  };

  Scalar DefaultMatrixFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
    Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const 
  {
    Scalar result = 0;
    if (gt == HERMES_PLANAR) {
      result = const_coeff * int_e_f<double, Scalar>(n, wt, u, v);
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  Ord DefaultMatrixFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      result = const_coeff * int_e_f<Ord, Ord>(n, wt, u, v);
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  MatrixFormVol<Scalar>* DefaultMatrixFormVol<Scalar>::clone() 
  {
    return new DefaultMatrixFormVol<Scalar>(*this);
  }


  DefaultJacobianCurlCurl::DefaultJacobianCurlCurl(int i, int j, std::string area, Scalar const_coeff,
    CubicSpline* c_spline,
    SymFlag sym, GeomType gt)
    : MatrixFormVol<Scalar>(i, j, area, sym), 
    idx_j(j), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
  {
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  };

  DefaultJacobianCurlCurl::DefaultJacobianCurlCurl(int i, int j, Hermes::vector<std::string> areas, 
    Scalar const_coeff, CubicSpline* c_spline,
    SymFlag sym, GeomType gt)
    : MatrixFormVol<Scalar>(i, j, areas, sym),
    idx_j(j), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
  {
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultJacobianCurlCurl::~DefaultJacobianCurlCurl() 
  {
    if (spline_coeff != HERMES_DEFAULT_SPLINE) delete spline_coeff;
  };

  Scalar DefaultJacobianCurlCurl::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
    Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const 
  {
    Scalar result = 0;
    if (gt == HERMES_PLANAR) {
      result = const_coeff * int_curl_e_curl_f<double, Scalar>(n, wt, u, v);
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  Ord DefaultJacobianCurlCurl::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      result = const_coeff * int_curl_e_curl_f<Ord, Ord>(n, wt, u, v);
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  MatrixFormVol<Scalar>* DefaultJacobianCurlCurl::clone() 
  {
    return new DefaultJacobianCurlCurl(*this);
  }


  DefaultVectorFormVol<Scalar>::DefaultVectorFormVol<Scalar>(int i, std::string area, Scalar const_coeff0, Scalar const_coeff1,
    DefaultFunction<Scalar>* f_coeff0, DefaultFunction<Scalar>* f_coeff1,
    GeomType gt)
    : VectorFormVol<Scalar>(i, area), const_coeff0(const_coeff0), const_coeff1(const_coeff1),
    function_coeff0(f_coeff0), function_coeff1(f_coeff1), gt(gt)
  { 
    // If f_coeff0 is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff0 == HERMES_DEFAULT_FUNCTION) this->function_coeff0 = new DefaultFunction(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
    if (f_coeff1 == HERMES_DEFAULT_FUNCTION) this->function_coeff1 = new DefaultFunction(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultVectorFormVol<Scalar>::DefaultVectorFormVol<Scalar>(int i, Hermes::vector<std::string> areas, 
    Scalar const_coeff0, Scalar const_coeff1,
    DefaultFunction<Scalar>* f_coeff0, DefaultFunction<Scalar>* f_coeff1,
    GeomType gt)
    : VectorFormVol<Scalar>(i, areas), const_coeff0(const_coeff0), const_coeff1(const_coeff1),
    function_coeff0(f_coeff0), function_coeff1(f_coeff1), gt(gt)
  { 
    // If f_coeff0 is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff0 == HERMES_DEFAULT_FUNCTION) this->function_coeff0 = new DefaultFunction(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
    if (f_coeff1 == HERMES_DEFAULT_FUNCTION) this->function_coeff1 = new DefaultFunction(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultVectorFormVol<Scalar>::~DefaultVectorFormVol<Scalar>() 
  {
    if (function_coeff0 != HERMES_DEFAULT_FUNCTION) delete function_coeff0;
    if (function_coeff1 != HERMES_DEFAULT_FUNCTION) delete function_coeff1;
  };

  Scalar DefaultVectorFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
    Geom<double> *e, ExtData<Scalar> *ext) const 
  {
    Scalar int_v0 = 0, int_v1 = 0;
    for (int i = 0; i < n; i++) int_v0 += wt[i] * v->val0[i];
    for (int i = 0; i < n; i++) int_v1 += wt[i] * v->val1[i];
    return const_coeff0 * int_v0 + const_coeff1 * int_v1;
  }

  Ord DefaultVectorFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord int_v0 = 0, int_v1 = 0;
    for (int i = 0; i < n; i++) int_v0 += wt[i] * v->val0[i];
    for (int i = 0; i < n; i++) int_v1 += wt[i] * v->val1[i];
    return const_coeff0 * int_v0 + const_coeff1 * int_v1;
  }

  VectorFormVol<Scalar>* DefaultVectorFormVol<Scalar>::clone() 
  {
    return new DefaultVectorFormVol<Scalar>(*this);
  }


  DefaultResidualVol::DefaultResidualVol(int i, std::string area, Scalar const_coeff,
    DefaultFunction<Scalar>* f_coeff,
    GeomType gt)
    : VectorFormVol<Scalar>(i, area),
    idx_i(i), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
    else error("Nonconstant functions in Hcurl forms not implemented yet.");
  }

  DefaultResidualVol::DefaultResidualVol(int i, Hermes::vector<std::string> areas, Scalar const_coeff,
    DefaultFunction<Scalar>* f_coeff,
    GeomType gt)
    : VectorFormVol<Scalar>(i, areas),
    idx_i(i), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
    else error("Nonconstant functions in Hcurl forms not implemented yet.");
  }

  DefaultResidualVol::~DefaultResidualVol() 
  {
    if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
  };

  Scalar DefaultResidualVol::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
    Geom<double> *e, ExtData<Scalar> *ext) const 
  {
    Scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * (u_ext[idx_i]->val0[i] * v->val0[i] +
          u_ext[idx_i]->val1[i] * v->val1[i]);
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return const_coeff * result;
  }

  Ord DefaultResidualVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  VectorFormVol<Scalar>* DefaultResidualVol::clone() 
  {
    return new DefaultResidualVol(*this);
  }


  DefaultResidualCurlCurl::DefaultResidualCurlCurl(int i, std::string area, Scalar const_coeff,
    CubicSpline* c_spline,
    GeomType gt)
    : VectorFormVol<Scalar>(i, area),
    idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
  {
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  };

  DefaultResidualCurlCurl::DefaultResidualCurlCurl(int i, Hermes::vector<std::string> areas, Scalar const_coeff,
    CubicSpline* c_spline, 
    GeomType gt)
    : VectorFormVol<Scalar>(i, areas),
    idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
  {
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultResidualCurlCurl::~DefaultResidualCurlCurl() 
  {
    if (spline_coeff != HERMES_DEFAULT_SPLINE) delete spline_coeff;
  };

  Scalar DefaultResidualCurlCurl::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
    Geom<double> *e, ExtData<Scalar> *ext) const
  {
    Func<Scalar>* u_prev = u_ext[idx_i];
    Scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        double mag0_i = std::abs(u_prev->val0[i]);
        double mag1_i = std::abs(u_prev->val1[i]);
        double mag_i = sqrt(sqr(mag0_i) + sqr(mag1_i));
        result += wt[i] * const_coeff*spline_coeff->get_value(mag_i) 
          * (u_prev->curl[i] * conj(v->curl[i]));
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  Ord DefaultResidualCurlCurl::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Func<Ord>* u_prev = u_ext[idx_i];
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        Ord mag0_i = u_prev->val0[i];
        Ord mag1_i = u_prev->val1[i];
        Ord mag_i = sqrt(sqr(mag0_i) + sqr(mag1_i));
        result += wt[i] * const_coeff*spline_coeff->get_value(mag_i) 
          * (u_prev->curl[i] * conj(v->curl[i]));
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  VectorFormVol<Scalar>* DefaultResidualCurlCurl::clone() 
  {
    return new DefaultResidualCurlCurl(*this);
  }


  DefaultMatrixFormSurf<Scalar>::DefaultMatrixFormSurf<Scalar>(int i, int j, std::string area,
    Scalar const_coeff, DefaultFunction<Scalar>* f_coeff,
    GeomType gt)
    : MatrixFormSurf<Scalar>(i, j, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
    else error("Nonconstant functions in Hcurl forms not implemented yet.");
  }

  DefaultMatrixFormSurf<Scalar>::DefaultMatrixFormSurf<Scalar>(int i, int j, Hermes::vector<std::string> areas,
    Scalar const_coeff, DefaultFunction<Scalar>* f_coeff,
    GeomType gt)
    : MatrixFormSurf<Scalar>(i, j, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
    else error("Nonconstant functions in Hcurl forms not implemented yet.");
  }

  DefaultMatrixFormSurf<Scalar>::~DefaultMatrixFormSurf<Scalar>() 
  {
    if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
  };

  Scalar DefaultMatrixFormSurf<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
    Geom<double> *e, ExtData<Scalar> *ext) const 
  {
    Scalar result = 0;
    if (gt == HERMES_PLANAR) {
      result = const_coeff * int_e_tau_f_tau<double, Scalar>(n, wt, u, v, e);
    }
    else error("Axisymmetric Hcurl forms not implemnted yet.");

    return result;
  }

  Ord DefaultMatrixFormSurf<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      result = const_coeff * int_e_tau_f_tau<Ord, Ord>(n, wt, u, v, e);
    }
    else error("Axisymmetric Hcurl forms not implemnted yet.");

    return result;
  }

  MatrixFormSurf<Scalar>* DefaultMatrixFormSurf<Scalar>::clone() 
  {
    return new DefaultMatrixFormSurf<Scalar>(*this);
  }


  DefaultResidualSurf::DefaultResidualSurf(int i, std::string area,
    Scalar const_coeff, DefaultFunction<Scalar>* f_coeff,
    GeomType gt)
    : VectorFormSurf<Scalar>(i, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
    else error("Nonconstant functions in Hcurl forms not implemented yet.");
  }

  DefaultResidualSurf::DefaultResidualSurf(int i, Hermes::vector<std::string> areas,
    Scalar const_coeff, DefaultFunction<Scalar>* f_coeff,
    GeomType gt)
    : VectorFormSurf<Scalar>(i, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
    else error("Nonconstant functions in Hcurl forms not implemented yet.");
  }

  DefaultResidualSurf::~DefaultResidualSurf()
  {
    if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
  };

  Scalar DefaultResidualSurf::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
    Geom<double> *e, ExtData<Scalar> *ext) const
  {
    Scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++)
        result += wt[i] * (    (u_ext[0]->val0[i] * e->tx[i] + u_ext[0]->val1[i] * e->ty[i]) *
        conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
      result *= const_coeff;
    }
    else error("Axisymmetric Hcurl forms not implemnted yet.");

    return result;
  }

  Ord DefaultResidualSurf::ord(int n, double *wt, Func<Ord> *u_ext[],
    Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++)
        result += wt[i] * (  (u_ext[0]->val0[i] * e->tx[i] + u_ext[0]->val1[i] * e->ty[i]) *
        conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
    }
    else error("Axisymmetric Hcurl forms not implemnted yet.");

    return result;
  }

  VectorFormSurf<Scalar>* DefaultResidualSurf::clone()
  {
    return new DefaultResidualSurf(*this);
  }


  DefaultVectorFormSurf<Scalar>::DefaultVectorFormSurf<Scalar>(int i, std::string area, Scalar const_coeff,
    DefaultFunction<Scalar>* f_coeff,
    GeomType gt)
    : VectorFormSurf<Scalar>(i, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultVectorFormSurf<Scalar>::DefaultVectorFormSurf<Scalar>(int i, Hermes::vector<std::string> areas, Scalar const_coeff,
    DefaultFunction<Scalar>* f_coeff,
    GeomType gt)
    : VectorFormSurf<Scalar>(i, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
    else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
  }

  DefaultVectorFormSurf<Scalar>::~DefaultVectorFormSurf<Scalar>() 
  {
    if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
  };

  Scalar DefaultVectorFormSurf<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
    Geom<double> *e, ExtData<Scalar> *ext) const 
  {
    Scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]);
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  Ord DefaultVectorFormSurf<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]);
      }
    }
    else error("Axisymmetric Hcurl forms not implemented yet.");

    return result;
  }

  VectorFormSurf<Scalar>* DefaultVectorFormSurf<Scalar>::clone() 
  {
    return new DefaultVectorFormSurf<Scalar>(*this);
  }


};
