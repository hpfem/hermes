// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "weakforms_hcurl.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace WeakFormsHcurl 
    {
      template<typename Scalar>
      DefaultMatrixFormVol<Scalar>::DefaultMatrixFormVol
        (int i, int j, std::string area, Scalar const_coeff, 
        HermesFunction<Scalar>* f_coeff, SymFlag sym, 
        GeomType gt)
        : MatrixFormVol<Scalar>(i, j, area, sym), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultMatrixFormVol<Scalar>::DefaultMatrixFormVol
        (int i, int j, Hermes::vector<std::string> areas,Scalar const_coeff, 
        HermesFunction<Scalar>* f_coeff, SymFlag sym, GeomType gt)
        : MatrixFormVol<Scalar>(i, j, areas, sym), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultMatrixFormVol<Scalar>::~DefaultMatrixFormVol() 
      {
        if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
      };

      template<typename Scalar>
      Scalar DefaultMatrixFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const 
      {
        Scalar result = 0;
        if (gt == HERMES_PLANAR) 
        {
          result = const_coeff * int_e_f<double, Scalar>(n, wt, u, v);
        }
        else error("Axisymmetric Hcurl forms not implemented yet.");

        return result;
      }

      template<typename Scalar>
      Ord DefaultMatrixFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        Ord result = 0;
        if (gt == HERMES_PLANAR) 
        {
          result = const_coeff * int_e_f<Ord, Ord>(n, wt, u, v);
        }
        else error("Axisymmetric Hcurl forms not implemented yet.");

        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultMatrixFormVol<Scalar>::clone() 
      {
        return new DefaultMatrixFormVol<Scalar>(*this);
      }


      template<typename Scalar>
      DefaultJacobianCurlCurl<Scalar>::DefaultJacobianCurlCurl(int i, int j, std::string area, Scalar const_coeff,
        CubicSpline* c_spline,
        SymFlag sym, GeomType gt)
        : MatrixFormVol<Scalar>(i, j, area, sym), 
        idx_j(j), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
      };

      template<typename Scalar>
      DefaultJacobianCurlCurl<Scalar>::DefaultJacobianCurlCurl(int i, int j, Hermes::vector<std::string> areas, 
        Scalar const_coeff, CubicSpline* c_spline,
        SymFlag sym, GeomType gt)
        : MatrixFormVol<Scalar>(i, j, areas, sym),
        idx_j(j), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultJacobianCurlCurl<Scalar>::~DefaultJacobianCurlCurl() 
      {
        if (spline_coeff != HERMES_DEFAULT_SPLINE) delete spline_coeff;
      };

      template<typename Scalar>
      Scalar DefaultJacobianCurlCurl<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, ExtData<Scalar> *ext) const 
      {
        Scalar result = 0;
        if (gt == HERMES_PLANAR) 
        {
          result = const_coeff * int_curl_e_curl_f<double, Scalar>(n, wt, u, v);
        }
        else error("Axisymmetric Hcurl forms not implemented yet.");

        return result;
      }

      template<typename Scalar>
      Ord DefaultJacobianCurlCurl<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        Ord result = 0;
        if (gt == HERMES_PLANAR) 
        {
          result = const_coeff * int_curl_e_curl_f<Ord, Ord>(n, wt, u, v);
        }
        else error("Axisymmetric Hcurl forms not implemented yet.");

        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultJacobianCurlCurl<Scalar>::clone() 
      {
        return new DefaultJacobianCurlCurl(*this);
      }


      template<typename Scalar>
      DefaultVectorFormVol<Scalar>::DefaultVectorFormVol(int i, std::string area, Scalar const_coeff0, Scalar const_coeff1,
        HermesFunction<Scalar>* f_coeff0, HermesFunction<Scalar>* f_coeff1,
        GeomType gt)
        : VectorFormVol<Scalar>(i, area), const_coeff0(const_coeff0), const_coeff1(const_coeff1),
        function_coeff0(f_coeff0), function_coeff1(f_coeff1), gt(gt)
      {
        // If f_coeff0 is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff0 == HERMES_DEFAULT_FUNCTION) this->function_coeff0 = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
        if (f_coeff1 == HERMES_DEFAULT_FUNCTION) this->function_coeff1 = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultVectorFormVol<Scalar>::DefaultVectorFormVol(int i, Hermes::vector<std::string> areas, 
        Scalar const_coeff0, Scalar const_coeff1,
        HermesFunction<Scalar>* f_coeff0, HermesFunction<Scalar>* f_coeff1,
        GeomType gt)
        : VectorFormVol<Scalar>(i, areas), const_coeff0(const_coeff0), const_coeff1(const_coeff1),
        function_coeff0(f_coeff0), function_coeff1(f_coeff1), gt(gt)
      {
        // If f_coeff0 is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff0 == HERMES_DEFAULT_FUNCTION) this->function_coeff0 = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
        if (f_coeff1 == HERMES_DEFAULT_FUNCTION) this->function_coeff1 = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultVectorFormVol<Scalar>::~DefaultVectorFormVol() 
      {
        if (function_coeff0 != HERMES_DEFAULT_FUNCTION) delete function_coeff0;
        if (function_coeff1 != HERMES_DEFAULT_FUNCTION) delete function_coeff1;
      };

      template<typename Scalar>
      Scalar DefaultVectorFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const 
      {
        Scalar int_v0 = 0, int_v1 = 0;
        for (int i = 0; i < n; i++) int_v0 += wt[i] * v->val0[i];
        for (int i = 0; i < n; i++) int_v1 += wt[i] * v->val1[i];
        return const_coeff0 * int_v0 + const_coeff1 * int_v1;
      }

      template<typename Scalar>
      Ord DefaultVectorFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        Ord int_v0 = 0, int_v1 = 0;
        for (int i = 0; i < n; i++) int_v0 += wt[i] * v->val0[i];
        for (int i = 0; i < n; i++) int_v1 += wt[i] * v->val1[i];
        return const_coeff0 * int_v0 + const_coeff1 * int_v1;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultVectorFormVol<Scalar>::clone() 
      {
        return new DefaultVectorFormVol<Scalar>(*this);
      }


      template<typename Scalar>
      DefaultResidualVol<Scalar>::DefaultResidualVol(int i, std::string area, Scalar const_coeff,
        HermesFunction<Scalar>* f_coeff,
        GeomType gt)
        : VectorFormVol<Scalar>(i, area),
        idx_i(i), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant functions in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultResidualVol<Scalar>::DefaultResidualVol(int i, Hermes::vector<std::string> areas, Scalar const_coeff,
        HermesFunction<Scalar>* f_coeff,
        GeomType gt)
        : VectorFormVol<Scalar>(i, areas),
        idx_i(i), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant functions in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultResidualVol<Scalar>::~DefaultResidualVol() 
      {
        if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
      };

      template<typename Scalar>
      Scalar DefaultResidualVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const 
      {
        Scalar result = 0;
        if (gt == HERMES_PLANAR) 
        {
          for (int i = 0; i < n; i++) 
          {
            result += wt[i] * function_coeff->value(e->x[i], e->y[i]) * (u_ext[idx_i]->val0[i] * v->val0[i] +
              u_ext[idx_i]->val1[i] * v->val1[i]);
          }
        }
        else error("Axisymmetric Hcurl forms not implemented yet.");

        return const_coeff * result;
      }

      template<typename Scalar>
      Ord DefaultResidualVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        Ord result = 0;
        if (gt == HERMES_PLANAR) 
        {
          for (int i = 0; i < n; i++) 
          {
            result += wt[i] * function_coeff->value_ord(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
          }
        }
        else error("Axisymmetric Hcurl forms not implemented yet.");

        return result;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultResidualVol<Scalar>::clone() 
      {
        return new DefaultResidualVol(*this);
      }


      template<typename Scalar>
      DefaultResidualCurlCurl<Scalar>::DefaultResidualCurlCurl(int i, std::string area, Scalar const_coeff,
        CubicSpline* c_spline,
        GeomType gt)
        : VectorFormVol<Scalar>(i, area),
        idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
      };

      template<typename Scalar>
      DefaultResidualCurlCurl<Scalar>::DefaultResidualCurlCurl(int i, Hermes::vector<std::string> areas, Scalar const_coeff,
        CubicSpline* c_spline, 
        GeomType gt)
        : VectorFormVol<Scalar>(i, areas),
        idx_i(i), const_coeff(const_coeff), spline_coeff(c_spline), gt(gt)
      {
        // If spline is HERMES_DEFAULT_SPLINE, initialize it to be constant 1.0.
        if (c_spline == HERMES_DEFAULT_SPLINE) this->spline_coeff = new CubicSpline(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultResidualCurlCurl<Scalar>::~DefaultResidualCurlCurl() 
      {
        if (spline_coeff != HERMES_DEFAULT_SPLINE) delete spline_coeff;
      };

      template<typename Scalar>
      Scalar DefaultResidualCurlCurl<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const
      {
        Func<Scalar>* u_prev = u_ext[idx_i];
        Scalar result = 0;
        if (gt == HERMES_PLANAR) 
        {
          for (int i = 0; i < n; i++) 
          {
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

      template<typename Scalar>
      Ord DefaultResidualCurlCurl<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        Func<Ord>* u_prev = u_ext[idx_i];
        Ord result = 0;
        if (gt == HERMES_PLANAR) 
        {
          for (int i = 0; i < n; i++) 
          {
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

      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultResidualCurlCurl<Scalar>::clone() 
      {
        return new DefaultResidualCurlCurl(*this);
      }


      template<typename Scalar>
      DefaultMatrixFormSurf<Scalar>::DefaultMatrixFormSurf(int i, int j, std::string area,
        Scalar const_coeff, HermesFunction<Scalar>* f_coeff,
        GeomType gt)
        : MatrixFormSurf<Scalar>(i, j, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant functions in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultMatrixFormSurf<Scalar>::DefaultMatrixFormSurf(int i, int j, Hermes::vector<std::string> areas,
        Scalar const_coeff, HermesFunction<Scalar>* f_coeff,
        GeomType gt)
        : MatrixFormSurf<Scalar>(i, j, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant functions in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultMatrixFormSurf<Scalar>::~DefaultMatrixFormSurf() 
      {
        if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
      };

      template<typename Scalar>
      Scalar DefaultMatrixFormSurf<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const 
      {
        Scalar result = 0;
        if (gt == HERMES_PLANAR) 
        {
          result = const_coeff * int_e_tau_f_tau<double, Scalar>(n, wt, u, v, e);
        }
        else error("Axisymmetric Hcurl forms not implemnted yet.");

        return result;
      }

      template<typename Scalar>
      Ord DefaultMatrixFormSurf<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        Ord result = 0;
        if (gt == HERMES_PLANAR) 
        {
          result = const_coeff * int_e_tau_f_tau<Ord, Ord>(n, wt, u, v, e);
        }
        else error("Axisymmetric Hcurl forms not implemnted yet.");

        return result;
      }

      template<typename Scalar>
      MatrixFormSurf<Scalar>* DefaultMatrixFormSurf<Scalar>::clone() 
      {
        return new DefaultMatrixFormSurf<Scalar>(*this);
      }


      template<typename Scalar>
      DefaultResidualSurf<Scalar>::DefaultResidualSurf(int i, std::string area,
        Scalar const_coeff, HermesFunction<Scalar>* f_coeff,
        GeomType gt)
        : VectorFormSurf<Scalar>(i, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant functions in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultResidualSurf<Scalar>::DefaultResidualSurf(int i, Hermes::vector<std::string> areas,
        Scalar const_coeff, HermesFunction<Scalar>* f_coeff,
        GeomType gt)
        : VectorFormSurf<Scalar>(i, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant functions in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultResidualSurf<Scalar>::~DefaultResidualSurf()
      {
        if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
      };

      template<typename Scalar>
      Scalar DefaultResidualSurf<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const
      {
        Scalar result = 0;
        if (gt == HERMES_PLANAR) 
        {
          for (int i = 0; i < n; i++)
            result += wt[i] * (    (u_ext[0]->val0[i] * e->tx[i] + u_ext[0]->val1[i] * e->ty[i]) *
            conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
          result *= const_coeff;
        }
        else error("Axisymmetric Hcurl forms not implemnted yet.");

        return result;
      }

      template<typename Scalar>
      Ord DefaultResidualSurf<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[],
        Func<Ord> *v, Geom<Ord> *e, ExtData<Ord> *ext) const
      {
        Ord result = 0;
        if (gt == HERMES_PLANAR) 
        {
          for (int i = 0; i < n; i++)
            result += wt[i] * (  (u_ext[0]->val0[i] * e->tx[i] + u_ext[0]->val1[i] * e->ty[i]) *
            conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]));
        }
        else error("Axisymmetric Hcurl forms not implemnted yet.");

        return result;
      }

      template<typename Scalar>
      VectorFormSurf<Scalar>* DefaultResidualSurf<Scalar>::clone()
      {
        return new DefaultResidualSurf(*this);
      }


      template<typename Scalar>
      DefaultVectorFormSurf<Scalar>::DefaultVectorFormSurf(int i, std::string area, Scalar const_coeff,
        HermesFunction<Scalar>* f_coeff,
        GeomType gt)
        : VectorFormSurf<Scalar>(i, area), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultVectorFormSurf<Scalar>::DefaultVectorFormSurf(int i, Hermes::vector<std::string> areas, Scalar const_coeff,
        HermesFunction<Scalar>* f_coeff,
        GeomType gt)
        : VectorFormSurf<Scalar>(i, areas), const_coeff(const_coeff), function_coeff(f_coeff), gt(gt)
      {
        // If f_coeff is HERMES_DEFAULT_FUNCTION, initialize it to be constant 1.0.
        if (f_coeff == HERMES_DEFAULT_FUNCTION) this->function_coeff = new HermesFunction<Scalar>(1.0);
        else error("Nonconstant coefficients in Hcurl forms not implemented yet.");
      }

      template<typename Scalar>
      DefaultVectorFormSurf<Scalar>::~DefaultVectorFormSurf() 
      {
        if (function_coeff != HERMES_DEFAULT_FUNCTION) delete function_coeff;
      };

      template<typename Scalar>
      Scalar DefaultVectorFormSurf<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, ExtData<Scalar> *ext) const 
      {
        Scalar result = 0;
        if (gt == HERMES_PLANAR) 
        {
          for (int i = 0; i < n; i++) 
          {
            result += wt[i] * conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]);
          }
        }
        else error("Axisymmetric Hcurl forms not implemented yet.");

        return result;
      }

      template<typename Scalar>
      Ord DefaultVectorFormSurf<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, ExtData<Ord> *ext) const 
      {
        Ord result = 0;
        if (gt == HERMES_PLANAR) 
        {
          for (int i = 0; i < n; i++) 
          {
            result += wt[i] * conj(v->val0[i] * e->tx[i] + v->val1[i] * e->ty[i]);
          }
        }
        else error("Axisymmetric Hcurl forms not implemented yet.");

        return result;
      }

      template<typename Scalar>
      VectorFormSurf<Scalar>* DefaultVectorFormSurf<Scalar>::clone() 
      {
        return new DefaultVectorFormSurf<Scalar>(*this);
      }

      template class HERMES_API DefaultJacobianCurlCurl<double>;
      template class HERMES_API DefaultJacobianCurlCurl<std::complex<double> >;
      template class HERMES_API DefaultMatrixFormVol<double>;
      template class HERMES_API DefaultMatrixFormVol<std::complex<double> >;
      template class HERMES_API DefaultVectorFormVol<double>;
      template class HERMES_API DefaultVectorFormVol<std::complex<double> >;
      template class HERMES_API DefaultResidualVol<double>;
      template class HERMES_API DefaultResidualVol<std::complex<double> >;
      template class HERMES_API DefaultResidualCurlCurl<double>;
      template class HERMES_API DefaultResidualCurlCurl<std::complex<double> >;
      template class HERMES_API DefaultMatrixFormSurf<double>;
      template class HERMES_API DefaultMatrixFormSurf<std::complex<double> >;
      template class HERMES_API DefaultVectorFormSurf<double>;
      template class HERMES_API DefaultVectorFormSurf<std::complex<double> >;
      template class HERMES_API DefaultResidualSurf<double>;
      template class HERMES_API DefaultResidualSurf<std::complex<double> >;
    };
  }
}