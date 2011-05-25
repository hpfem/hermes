#include "../hermes2d.h"

namespace WeakFormsH1 
{
  DefaultMatrixFormVol::DefaultMatrixFormVol
    (int i, int j, std::string area, scalar const_coeff, 
    DefaultFunction* f_coeff, SymFlag sym, 
    GeomType gt)
    : WeakForm::MatrixFormVol(i, j, area, sym), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }

  DefaultMatrixFormVol::DefaultMatrixFormVol
    (int i, int j, Hermes::vector<std::string> areas,scalar const_coeff, 
    DefaultFunction* f_coeff, SymFlag sym, GeomType gt)
    : WeakForm::MatrixFormVol(i, j, areas, sym), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }

  DefaultMatrixFormVol::~DefaultMatrixFormVol() 
  {
    if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
  };

  scalar DefaultMatrixFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                                     Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
    }

    return const_coeff * result;
  }

  Ord DefaultMatrixFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
    }

    return result;
  }

  WeakForm::MatrixFormVol* DefaultMatrixFormVol::clone() 
  {
    return new DefaultMatrixFormVol(*this);
  }


  DefaultJacobianDiffusion::DefaultJacobianDiffusion(int i, int j, std::string area, scalar const_coeff,
                                                     CubicSpline* c_spline,
                                                     SymFlag sym, GeomType gt)
    : WeakForm::MatrixFormVol(i, j, area, sym), 
      idx_j(j), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
  {
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
  };

  DefaultJacobianDiffusion::DefaultJacobianDiffusion(int i, int j, Hermes::vector<std::string> areas, 
                                                     scalar const_coeff, CubicSpline* c_spline,
                                                     SymFlag sym, GeomType gt)
    : WeakForm::MatrixFormVol(i, j, areas, sym),
      idx_j(j), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
  {
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
  }

  DefaultJacobianDiffusion::~DefaultJacobianDiffusion() 
  {
    if (spline_coeff != HERMES_DEFAULT_SPLINE) delete spline_coeff;
  };

  scalar DefaultJacobianDiffusion::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * (const_coeff*spline_coeff->get_derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                  (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                  + const_coeff*spline_coeff->get_value(u_ext[idx_j]->val[i])
            * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * (const_coeff*spline_coeff->get_derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                    (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                    + const_coeff*spline_coeff->get_value(u_ext[idx_j]->val[i])
            * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * (const_coeff*spline_coeff->get_derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                    (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                    + const_coeff*spline_coeff->get_value(u_ext[idx_j]->val[i])
            * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
        }
      }
    }

    return result;
  }

  Ord DefaultJacobianDiffusion::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                    Geom<Ord> *e, ExtData<Ord> *ext) const
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * (const_coeff*spline_coeff->get_derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                  (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                  + const_coeff*spline_coeff->get_value(u_ext[idx_j]->val[i])
                                  * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
                result += wt[i] * e->y[i] * (const_coeff*spline_coeff->get_derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                          (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                          + const_coeff*spline_coeff->get_value(u_ext[idx_j]->val[i])
              * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * (const_coeff*spline_coeff->get_derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                    (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                    + const_coeff*spline_coeff->get_value(u_ext[idx_j]->val[i])
                      * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
        }
      }
    }

    return result;
  }

  WeakForm::MatrixFormVol* DefaultJacobianDiffusion::clone() 
  {
    return new DefaultJacobianDiffusion(*this);
  }
  

  DefaultJacobianAdvection::DefaultJacobianAdvection(int i, int j, std::string area, 
                                                     scalar const_coeff1, scalar const_coeff2,
                                                     CubicSpline* c_spline1,
                                                     CubicSpline* c_spline2, GeomType gt)
    : WeakForm::MatrixFormVol(i, j, area, HERMES_NONSYM),
      idx_j(j), const_coeff1(const_coeff1), const_coeff2(const_coeff2), 
      spline_coeff1(c_spline1), spline_coeff2(c_spline2), gt(gt)
  {
    if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

    // If spline1 == HERMES_DEFAULT_SPLINE or spline2 == HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline1 == HERMES_DEFAULT_SPLINE) this->spline_coeff1 = new CubicSpline(1.0);
    if (c_spline2 == HERMES_DEFAULT_SPLINE) this->spline_coeff2 = new CubicSpline(1.0);
  }

  DefaultJacobianAdvection::DefaultJacobianAdvection(int i, int j, Hermes::vector<std::string> areas, 
                                                     scalar const_coeff1, scalar const_coeff2,
                                                     CubicSpline* c_spline1,
                                                     CubicSpline* c_spline2,
                                                     GeomType gt)
    : WeakForm::MatrixFormVol(i, j, areas, HERMES_NONSYM),
      idx_j(j), const_coeff1(const_coeff1), const_coeff2(const_coeff2), 
      spline_coeff1(spline_coeff1), spline_coeff2(spline_coeff2), gt(gt)
  {
    if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

    // If spline1 == HERMES_DEFAULT_SPLINE or spline2 == HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline1 == HERMES_DEFAULT_SPLINE) this->spline_coeff1 = new CubicSpline(1.0);
    if (c_spline2 == HERMES_DEFAULT_SPLINE) this->spline_coeff2 = new CubicSpline(1.0);
  }

  DefaultJacobianAdvection::~DefaultJacobianAdvection() 
  {
    if (spline_coeff1 != HERMES_DEFAULT_SPLINE) delete spline_coeff1;
    if (spline_coeff2 != HERMES_DEFAULT_SPLINE) delete spline_coeff2;
  };

  template<typename Real, typename Scalar>
  Scalar DefaultJacobianAdvection::matrix_form(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                               Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
  {
    Scalar result = 0;
        for (int i = 0; i < n; i++) {
          result += wt[i] * (  const_coeff1*spline_coeff1->get_derivative(u_ext[idx_j]->val[i]) * u->val[i] * u_ext[idx_j]->dx[i] * v->val[i]
                             + const_coeff1*spline_coeff1->get_value(u_ext[idx_j]->val[i]) * u->dx[i] * v->val[i]
                             + const_coeff2*spline_coeff2->get_derivative(u_ext[idx_j]->val[i]) * u->val[i] * u_ext[idx_j]->dy[i] * v->val[i]
                             + const_coeff2*spline_coeff2->get_value(u_ext[idx_j]->val[i]) * u->dy[i] * v->val[i]);
        }
        return result;
      }

  scalar DefaultJacobianAdvection::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u,
                                         Func<double> *v, Geom<double> *e, ExtData<scalar> *ext) const 
  {
    return matrix_form<double, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord DefaultJacobianAdvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
                                    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    return matrix_form<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  // This is to make the form usable in rk_time_step().
  WeakForm::MatrixFormVol* DefaultJacobianAdvection::clone() 
  {
    return new DefaultJacobianAdvection(*this);
  }


  DefaultVectorFormVol::DefaultVectorFormVol(int i, std::string area, scalar const_coeff,
                                             DefaultFunction* f_coeff,
                                             GeomType gt)
    : WeakForm::VectorFormVol(i, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }

  DefaultVectorFormVol::DefaultVectorFormVol(int i, Hermes::vector<std::string> areas, scalar const_coeff,
                                             DefaultFunction* f_coeff,
                                             GeomType gt)
    : WeakForm::VectorFormVol(i, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }

  DefaultVectorFormVol::~DefaultVectorFormVol() 
  {
    if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
  };

  scalar DefaultVectorFormVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                     Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
            result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
    }
    return const_coeff * result;
  }

  Ord DefaultVectorFormVol::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
        }
      }
    }

    return result;
  }

  WeakForm::VectorFormVol* DefaultVectorFormVol::clone() 
  {
    return new DefaultVectorFormVol(*this);
  }


  DefaultResidualVol::DefaultResidualVol(int i, std::string area, scalar const_coeff,
                                         DefaultFunction* f_coeff,
                                         GeomType gt)
    : WeakForm::VectorFormVol(i, area),
      idx_i(i), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }

  DefaultResidualVol::DefaultResidualVol(int i, Hermes::vector<std::string> areas, scalar const_coeff,
                                         DefaultFunction* f_coeff,
                                         GeomType gt)
    : WeakForm::VectorFormVol(i, areas),
      idx_i(i), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }

  DefaultResidualVol::~DefaultResidualVol() 
  {
    if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
  };

  scalar DefaultResidualVol::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                   Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
    }
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
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
    }

    return result;
  }

  WeakForm::VectorFormVol* DefaultResidualVol::clone() 
  {
    return new DefaultResidualVol(*this);
  }

    
  DefaultResidualDiffusion::DefaultResidualDiffusion(int i, std::string area, scalar const_coeff,
                                                     CubicSpline* c_spline,
                                                     GeomType gt)
    : WeakForm::VectorFormVol(i, area),
      idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
  {
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
  };

  DefaultResidualDiffusion::DefaultResidualDiffusion(int i, Hermes::vector<std::string> areas, scalar const_coeff,
                                                     CubicSpline* c_spline, 
                                                     GeomType gt)
    : WeakForm::VectorFormVol(i, areas),
      idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
  {
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
  }

  DefaultResidualDiffusion::~DefaultResidualDiffusion() 
  {
    if (spline_coeff != HERMES_DEFAULT_SPLINE) delete spline_coeff;
  };

  scalar DefaultResidualDiffusion::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                         Geom<double> *e, ExtData<scalar> *ext) const
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * const_coeff * spline_coeff->get_value(u_ext[idx_i]->val[i])
                        * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * const_coeff * spline_coeff->get_value(u_ext[idx_i]->val[i])
                          * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * const_coeff * spline_coeff->get_value(u_ext[idx_i]->val[i])
                          * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
        }
      }
    }

    return result;
  }

  Ord DefaultResidualDiffusion::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    // Planar base.
    for (int i = 0; i < n; i++) {
      result += wt[i] * const_coeff * spline_coeff->get_value(u_ext[idx_i]->val[i])
                      * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
    }
    if (gt != HERMES_PLANAR) result = result * Ord(1);

    return result;
  }

  WeakForm::VectorFormVol* DefaultResidualDiffusion::clone() 
  {
    return new DefaultResidualDiffusion(*this);
  }
 
  
  DefaultResidualAdvection::DefaultResidualAdvection(int i, std::string area, 
                                                     scalar const_coeff1, scalar const_coeff2, 
                                                     CubicSpline* c_spline1,
                                                     CubicSpline* c_spline2,
                                                     GeomType gt)
    : WeakForm::VectorFormVol(i, area),
      idx_i(i), const_coeff1(const_coeff1), const_coeff2(const_coeff2), 
      spline_coeff1(c_spline1), spline_coeff2(c_spline2), gt(gt)
  {
    if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

    // If spline1 == HERMES_DEFAULT_SPLINE or spline2 == HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline1 == HERMES_DEFAULT_SPLINE) this->spline_coeff1 = new CubicSpline(1.0);
    if (c_spline2 == HERMES_DEFAULT_SPLINE) this->spline_coeff2 = new CubicSpline(1.0);
  }
  
  DefaultResidualAdvection::DefaultResidualAdvection(int i, Hermes::vector<std::string> areas,\
                                                     scalar const_coeff1, scalar const_coeff2,
                                                     CubicSpline* c_spline1,
                                                     CubicSpline* c_spline2, GeomType gt)
    : WeakForm::VectorFormVol(i, areas),
      idx_i(i), const_coeff1(const_coeff1), const_coeff2(const_coeff2), 
      spline_coeff1(spline_coeff1), spline_coeff2(spline_coeff2), gt(gt)
  {
    if (gt != HERMES_PLANAR) error("Axisymmetric advection forms not implemented yet.");

    // If spline1 == HERMES_DEFAULT_SPLINE or spline2 == HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline1 == HERMES_DEFAULT_SPLINE) this->spline_coeff1 = new CubicSpline(1.0);
    if (c_spline2 == HERMES_DEFAULT_SPLINE) this->spline_coeff2 = new CubicSpline(1.0);
  }

  DefaultResidualAdvection::~DefaultResidualAdvection() 
  {
    if (spline_coeff1 != HERMES_DEFAULT_SPLINE) delete spline_coeff1;
    if (spline_coeff2 != HERMES_DEFAULT_SPLINE) delete spline_coeff2;
  };

  template<typename Real, typename Scalar>
  Scalar DefaultResidualAdvection::vector_form(int n, double *wt, Func<Scalar> *u_ext[],
                                               Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
  {
    Scalar result = 0;
    Func<Scalar>* u_prev = u_ext[idx_i];
    for (int i = 0; i < n; i++) {
      result += wt[i] * (const_coeff1*spline_coeff1->get_value(u_prev->val[i]) * (u_prev->dx[i] * v->val[i])
                          + const_coeff2*spline_coeff2->get_value(u_prev->val[i]) * (u_prev->dy[i] * v->val[i]));
    }
    return result;
  }

  scalar DefaultResidualAdvection::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                         Geom<double> *e, ExtData<scalar> *ext) const 
  {
    return vector_form<double, scalar>(n, wt, u_ext, v, e, ext);
  }

  Ord DefaultResidualAdvection::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                    Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    return vector_form<Ord, Ord>(n, wt, u_ext, v, e, ext);
  }

  WeakForm::VectorFormVol* DefaultResidualAdvection::clone() 
  {
    return new DefaultResidualAdvection(*this);
  }
 
  
  DefaultMatrixFormSurf::DefaultMatrixFormSurf(int i, int j, std::string area,
                                               scalar const_coeff, DefaultFunction* f_coeff,
                                               GeomType gt)
    : WeakForm::MatrixFormSurf(i, j, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }
  
  DefaultMatrixFormSurf::DefaultMatrixFormSurf(int i, int j, Hermes::vector<std::string> areas,
                                               scalar const_coeff, DefaultFunction* f_coeff,
                                               GeomType gt)
    : WeakForm::MatrixFormSurf(i, j, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }

  DefaultMatrixFormSurf::~DefaultMatrixFormSurf() 
  {
    if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
  };

  scalar DefaultMatrixFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                                      Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
    }

    return const_coeff * result;
  }

  Ord DefaultMatrixFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                 Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * u->val[i] * v->val[i];
        }
      }
    }

    return result;
  }

  // This is to make the form usable in rk_time_step().
  WeakForm::MatrixFormSurf* DefaultMatrixFormSurf::clone() 
  {
    return new DefaultMatrixFormSurf(*this);
  }
 
  
  DefaultJacobianFormSurf::DefaultJacobianFormSurf(int i, int j, std::string area, scalar const_coeff,
                                                   CubicSpline* c_spline,
                                                   GeomType gt)
    : WeakForm::MatrixFormSurf(i, j, area),
      idx_j(j), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
  {
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
  }
  
  DefaultJacobianFormSurf::DefaultJacobianFormSurf(int i, int j, Hermes::vector<std::string> areas, scalar const_coeff,
                                                   CubicSpline* c_spline,
                                                   GeomType gt)
    : WeakForm::MatrixFormSurf(i, j, areas), const_coeff(const_coeff), spline_coeff(spline_coeff), gt(gt)
  {
    // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
    if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
  }

  DefaultJacobianFormSurf::~DefaultJacobianFormSurf() 
  {
    if (spline_coeff != HERMES_DEFAULT_SPLINE) delete spline_coeff;
  };

  template<typename Real, typename Scalar>
  Scalar DefaultJacobianFormSurf::matrix_form_surf(int n, double *wt, Func<Scalar> *u_ext[], Func<Real> *u,
                                                   Func<Real> *v, Geom<Real> *e, ExtData<Scalar> *ext) const 
  {
    Scalar result = 0;
    for (int i = 0; i < n; i++) {
      result += wt[i] * (const_coeff*spline_coeff->get_derivative(u_ext[idx_j]->val[i]) * u_ext[idx_j]->val[i]
                          + const_coeff*spline_coeff->get_value(u_ext[idx_j]->val[i]))
                * u->val[i] * v->val[i];
    }
    return result;
  }

  scalar DefaultJacobianFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *u, Func<double> *v,
                                        Geom<double> *e, ExtData<scalar> *ext) const 
  {
    return matrix_form_surf<double, scalar>(n, wt, u_ext, u, v, e, ext);
  }

  Ord DefaultJacobianFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
                                   Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    return matrix_form_surf<Ord, Ord>(n, wt, u_ext, u, v, e, ext);
  }

  // This is to make the form usable in rk_time_step().
  WeakForm::MatrixFormSurf* DefaultJacobianFormSurf::clone() 
  {
    return new DefaultJacobianFormSurf(*this);
  }


  DefaultVectorFormSurf::DefaultVectorFormSurf(int i, std::string area, scalar const_coeff,
                                               DefaultFunction* f_coeff,
                                               GeomType gt)
    : WeakForm::VectorFormSurf(i, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }
  
  DefaultVectorFormSurf::DefaultVectorFormSurf(int i, Hermes::vector<std::string> areas, scalar const_coeff,
                                               DefaultFunction* f_coeff,
                                               GeomType gt)
    : WeakForm::VectorFormSurf(i, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }

  DefaultVectorFormSurf::~DefaultVectorFormSurf() 
  {
    if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
  };

  scalar DefaultVectorFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                      Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * v->val[i];
        }
      }
    }

    return const_coeff * result;
  }

  Ord DefaultVectorFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                 Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * v->val[i];
        }
      }
    }

    return result;
  }

  WeakForm::VectorFormSurf* DefaultVectorFormSurf::clone() 
  {
    return new DefaultVectorFormSurf(*this);
  }
 
  
  DefaultMultiComponentVectorFormSurf::DefaultMultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates,
                                                                           std::string area,
                                                                           Hermes::vector<scalar> coeffs,
                                                                           GeomType gt)
  : WeakForm::MultiComponentVectorFormSurf(coordinates, area), coeffs(coeffs), gt(gt) 
  { 
  }
  
  DefaultMultiComponentVectorFormSurf::DefaultMultiComponentVectorFormSurf(Hermes::vector<unsigned int> coordinates,
                                                                           Hermes::vector<std::string> areas,
                                                                           Hermes::vector<scalar> coeffs, GeomType gt)
  : WeakForm::MultiComponentVectorFormSurf(coordinates, areas), coeffs(coeffs), gt(gt) 
  { 
  }

  void DefaultMultiComponentVectorFormSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                                  Geom<double> *e, ExtData<scalar> *ext, 
                                                  Hermes::vector<scalar>& result) const 
  {
    scalar result_base = 0;
    if (gt == HERMES_PLANAR)
      result_base = int_v<double>(n, wt, v);
    else
      if (gt == HERMES_AXISYM_X)
        result_base = int_y_v<double>(n, wt, v, e);
      else
        result_base = int_x_v<double>(n, wt, v, e);

    for(unsigned int result_i = 0; result_i < this->coordinates.size(); result_i++)
      result.push_back(result_base * coeffs[result_i]);
  }

  Ord DefaultMultiComponentVectorFormSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                                               Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    if (gt == HERMES_PLANAR)
      return int_v<Ord>(n, wt, v);
    else
      if (gt == HERMES_AXISYM_X)
        return int_y_v<Ord>(n, wt, v, e);
      else
        return int_x_v<Ord>(n, wt, v, e);
  }

  /*
  WeakForm::VectorFormSurf* DefaultMultiComponentVectorFormSurf::clone() {
    return new DefaultMultiComponentVectorFormSurf(*this);
  }
  */

  DefaultResidualSurf::DefaultResidualSurf(int i, std::string area, scalar const_coeff,
                                           DefaultFunction* f_coeff,
                                           GeomType gt)
    : WeakForm::VectorFormSurf(i, area),
      idx_i(i), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }
  
  DefaultResidualSurf::DefaultResidualSurf(int i, Hermes::vector<std::string> areas, scalar const_coeff,
                                           DefaultFunction* f_coeff,
                                           GeomType gt)
    : WeakForm::VectorFormSurf(i, areas),
      idx_i(i), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
  {
    // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
    if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new DefaultFunction(1.0);
  }

  DefaultResidualSurf::~DefaultResidualSurf() 
  {
    if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
  };

  scalar DefaultResidualSurf::value(int n, double *wt, Func<scalar> *u_ext[], Func<double> *v,
                                    Geom<double> *e, ExtData<scalar> *ext) const 
  {
    scalar result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
    }

    return const_coeff * result;
  }

  Ord DefaultResidualSurf::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
                               Geom<Ord> *e, ExtData<Ord> *ext) const 
  {
    Ord result = 0;
    if (gt == HERMES_PLANAR) {
      for (int i = 0; i < n; i++) {
        result += wt[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
      }
    }
    else {
      if (gt == HERMES_AXISYM_X) {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->y[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
      else {
        for (int i = 0; i < n; i++) {
          result += wt[i] * e->x[i] * function_coeff->ord(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
        }
      }
    }

    return result;
  }

  // This is to make the form usable in rk_time_step().
  WeakForm::VectorFormSurf* DefaultResidualSurf::clone() 
  {
    return new DefaultResidualSurf(*this);
  }

  
  DefaultWeakFormLaplace::DefaultWeakFormLaplace(std::string area, scalar const_coeff,
                                                 CubicSpline* spline_coeff,
                                                 GeomType gt) : WeakForm()
  {
    // Jacobian.
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, area, const_coeff,
                                                  spline_coeff, HERMES_SYM, gt));
    // Residual.
    add_vector_form(new DefaultResidualDiffusion(0, area, const_coeff,
                                                  spline_coeff, gt));
  };
  
  
  DefaultWeakFormPoisson::DefaultWeakFormPoisson(DefaultFunction* rhs,
                                                 std::string area, scalar const_coeff,
                                                 CubicSpline* spline_coeff,
                                                 GeomType gt) : WeakForm()
  {
    // Jacobian.
    add_matrix_form(new DefaultJacobianDiffusion(0, 0, area, const_coeff,
                                                 spline_coeff, HERMES_NONSYM, gt));
    // Residual.
    add_vector_form(new DefaultResidualDiffusion(0, area, const_coeff,
                                                 spline_coeff, gt));
    add_vector_form(new DefaultVectorFormVol(0, area, -1.0, rhs, gt));
  };
};
