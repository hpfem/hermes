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

#include "weakforms_h1.h"
#include "api2d.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace WeakFormsH1
    {
      template<>
      DefaultMatrixFormVol<double>::DefaultMatrixFormVol
        (int i, int j, std::string area, Hermes2DFunction<double>* coeff, SymFlag sym, GeomType gt)
        : MatrixFormVol<double>(i, j), coeff(coeff), gt(gt)
      {
        this->set_area(area);
        this->setSymFlag(sym);

        if(coeff == HERMES_ONE)
          this->coeff = new Hermes2DFunction<double>(1.0);
      }
      template<>
      DefaultMatrixFormVol<std::complex<double> >::DefaultMatrixFormVol
        (int i, int j, std::string area, Hermes2DFunction<std::complex<double> >* coeff, SymFlag sym, GeomType gt)
        : MatrixFormVol<std::complex<double> >(i, j), coeff(coeff), gt(gt)
      {
        this->set_area(area);
        this->setSymFlag(sym);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE)
          this->coeff = new Hermes2DFunction<std::complex<double> >(std::complex<double>(1.0, 1.0));
      }

      template<>
      DefaultMatrixFormVol<double>::DefaultMatrixFormVol
        (int i, int j, Hermes::vector<std::string> areas,
        Hermes2DFunction<double>* coeff, SymFlag sym, GeomType gt)
        : MatrixFormVol<double>(i, j), coeff(coeff), gt(gt)
      {
        this->set_areas(areas);
        this->setSymFlag(sym);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE)
          this->coeff = new Hermes2DFunction<double>(1.0);
      }
      template<>
      DefaultMatrixFormVol<std::complex<double> >::DefaultMatrixFormVol
        (int i, int j, Hermes::vector<std::string> areas,
        Hermes2DFunction<std::complex<double> >* coeff, SymFlag sym, GeomType gt)
        : MatrixFormVol<std::complex<double> >(i, j), coeff(coeff), gt(gt)
      {
        this->set_areas(areas);
        this->setSymFlag(sym);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE)
          this->coeff = new Hermes2DFunction<std::complex<double> >(std::complex<double>(1.0, 1.0));
      }

      template<typename Scalar>
      DefaultMatrixFormVol<Scalar>::~DefaultMatrixFormVol()
      {
        
        if(coeff == HERMES_ONE) delete coeff;
      };

      template<typename Scalar>
      Scalar DefaultMatrixFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      Ord DefaultMatrixFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultMatrixFormVol<Scalar>::clone() const
      {
        return new DefaultMatrixFormVol<Scalar>(*this);
      }

      template<typename Scalar>
      DefaultJacobianDiffusion<Scalar>::DefaultJacobianDiffusion(int i, int j, std::string area,
        Hermes1DFunction<Scalar>* coeff,
        SymFlag sym, GeomType gt)
        : MatrixFormVol<Scalar>(i, j), idx_j(j), coeff(coeff), gt(gt)
      {
        this->set_area(area);
        this->setSymFlag(sym);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE)
          this->coeff = new Hermes1DFunction<Scalar>(1.0);
      };

      template<typename Scalar>
      DefaultJacobianDiffusion<Scalar>::DefaultJacobianDiffusion(int i, int j, Hermes::vector<std::string> areas,
        Hermes1DFunction<Scalar>* coeff, SymFlag sym, GeomType gt)
        : MatrixFormVol<Scalar>(i, j), idx_j(j), coeff(coeff), gt(gt)
      {
        this->set_areas(areas);
        this->setSymFlag(sym);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes1DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultJacobianDiffusion<Scalar>::~DefaultJacobianDiffusion()
      {
        if(coeff==HERMES_ONE) delete coeff;
      };

      template<typename Scalar>
      Scalar DefaultJacobianDiffusion<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
              (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
              + coeff->value(u_ext[idx_j]->val[i])
              * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                + coeff->value(u_ext[idx_j]->val[i])
                * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                + coeff->value(u_ext[idx_j]->val[i])
                * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      Ord DefaultJacobianDiffusion<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
              (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
              + coeff->value(u_ext[idx_j]->val[i])
              * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                + coeff->value(u_ext[idx_j]->val[i])
                * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u->val[i] *
                (u_ext[idx_j]->dx[i] * v->dx[i] + u_ext[idx_j]->dy[i] * v->dy[i])
                + coeff->value(u_ext[idx_j]->val[i])
                * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]));
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultJacobianDiffusion<Scalar>::clone() const
      {
        return new DefaultJacobianDiffusion<Scalar>(*this);
      }

      template<typename Scalar>
      DefaultMatrixFormDiffusion<Scalar>::DefaultMatrixFormDiffusion(int i, int j, std::string area,
        Hermes1DFunction<Scalar>* coeff,
        SymFlag sym, GeomType gt)
        : MatrixFormVol<Scalar>(i, j), idx_j(j), coeff(coeff), gt(gt)
      {
        this->set_area(area);
        this->setSymFlag(sym);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes1DFunction<Scalar>(1.0);
      };

      template<typename Scalar>
      DefaultMatrixFormDiffusion<Scalar>::DefaultMatrixFormDiffusion(int i, int j, Hermes::vector<std::string> areas,
        Hermes1DFunction<Scalar>* coeff, SymFlag sym, GeomType gt)
        : MatrixFormVol<Scalar>(i, j), idx_j(j), coeff(coeff), gt(gt)
      {
        this->set_areas(areas);
        this->setSymFlag(sym);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes1DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultMatrixFormDiffusion<Scalar>::~DefaultMatrixFormDiffusion()
      {
        if(coeff == HERMES_ONE) delete coeff;
      };

      template<typename Scalar>
      Scalar DefaultMatrixFormDiffusion<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      Ord DefaultMatrixFormDiffusion<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * (u->dx[i] * v->dx[i] + u->dy[i] * v->dy[i]);
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultMatrixFormDiffusion<Scalar>::clone() const
      {
        return new DefaultMatrixFormDiffusion<Scalar>(*this);
      }

      template<typename Scalar>
      DefaultJacobianAdvection<Scalar>::DefaultJacobianAdvection(int i, int j, std::string area,
        Hermes1DFunction<Scalar>* coeff1,
        Hermes1DFunction<Scalar>* coeff2,
        GeomType gt)
        : MatrixFormVol<Scalar>(i, j),
        idx_j(j), coeff1(coeff1), coeff2(coeff2), gt(gt)
      {
        this->set_area(area);
        
        if(gt != HERMES_PLANAR) throw Hermes::Exceptions::Exception("Axisymmetric advection forms not implemented yet.");

        // If coeff1 == HERMES_ONE or coeff22 == HERMES_ONE, initialize it to be constant 1.0.
        if(coeff1 == HERMES_ONE) this->coeff1 = new Hermes1DFunction<Scalar>(1.0);
        if(coeff2 == HERMES_ONE) this->coeff2 = new Hermes1DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultJacobianAdvection<Scalar>::DefaultJacobianAdvection(int i, int j, Hermes::vector<std::string> areas,
        Hermes1DFunction<Scalar>* coeff1,
        Hermes1DFunction<Scalar>* coeff2,
        GeomType gt)
        : MatrixFormVol<Scalar>(i, j),
        idx_j(j), coeff1(coeff1), coeff2(coeff2), gt(gt)
      {
        this->set_areas(areas);
        
        if(gt != HERMES_PLANAR) throw Hermes::Exceptions::Exception("Axisymmetric advection forms not implemented yet.");

        // If coeff1 == HERMES_ONE or coeff22 == HERMES_ONE, initialize it to be constant 1.0.
        if(coeff1 == HERMES_ONE) this->coeff1 = new Hermes1DFunction<Scalar>(1.0);
        if(coeff2 == HERMES_ONE) this->coeff2 = new Hermes1DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultJacobianAdvection<Scalar>::~DefaultJacobianAdvection()
      {
        
        //if(coeff1 != HERMES_ONE) delete coeff1;
        //if(coeff2 != HERMES_ONE) delete coeff2;
      };

      template<typename Scalar>
      Scalar DefaultJacobianAdvection<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u,
        Func<double> *v, Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++) {
          result += wt[i] * (  coeff1->derivative(u_ext[idx_j]->val[i]) * u->val[i] * u_ext[idx_j]->dx[i] * v->val[i]
          + coeff1->value(u_ext[idx_j]->val[i]) * u->dx[i] * v->val[i]
          + coeff2->derivative(u_ext[idx_j]->val[i]) * u->val[i] * u_ext[idx_j]->dy[i] * v->val[i]
          + coeff2->value(u_ext[idx_j]->val[i]) * u->dy[i] * v->val[i]);
        }
        return result;
      }

      template<typename Scalar>
      Ord DefaultJacobianAdvection<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u, Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++) {
          result += wt[i] * (  coeff1->derivative(u_ext[idx_j]->val[i]) * u->val[i] * u_ext[idx_j]->dx[i] * v->val[i]
          + coeff1->value(u_ext[idx_j]->val[i]) * u->dx[i] * v->val[i]
          + coeff2->derivative(u_ext[idx_j]->val[i]) * u->val[i] * u_ext[idx_j]->dy[i] * v->val[i]
          + coeff2->value(u_ext[idx_j]->val[i]) * u->dy[i] * v->val[i]);
        }
        return result;
      }

      // This is to make the form usable in rk_time_step_newton().
      template<typename Scalar>
      MatrixFormVol<Scalar>* DefaultJacobianAdvection<Scalar>::clone() const
      {
        return new DefaultJacobianAdvection<Scalar>(*this);
      }

      template<typename Scalar>
      DefaultVectorFormVol<Scalar>::DefaultVectorFormVol(int i, std::string area,
        Hermes2DFunction<Scalar>* coeff,
        GeomType gt)
        : VectorFormVol<Scalar>(i), coeff(coeff), gt(gt)
      {
        this->set_area(area);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes2DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultVectorFormVol<Scalar>::DefaultVectorFormVol(int i, Hermes::vector<std::string> areas,
        Hermes2DFunction<Scalar>* coeff,
        GeomType gt)
        : VectorFormVol<Scalar>(i), coeff(coeff), gt(gt)
      {
        this->set_areas(areas);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes2DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultVectorFormVol<Scalar>::~DefaultVectorFormVol()
      {
        
        if(coeff == HERMES_ONE) delete coeff;
      };

      template<typename Scalar>
      Scalar DefaultVectorFormVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
            }
          }
        }
        return result;
      }

      template<typename Scalar>
      Ord DefaultVectorFormVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultVectorFormVol<Scalar>::clone() const
      {
        return new DefaultVectorFormVol<Scalar>(*this);
      }

      template<typename Scalar>
      DefaultResidualVol<Scalar>::DefaultResidualVol(int i, std::string area,
        Hermes2DFunction<Scalar>* coeff,
        GeomType gt)
        : VectorFormVol<Scalar>(i), idx_i(i), coeff(coeff), gt(gt)
      {
        this->set_area(area);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes2DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultResidualVol<Scalar>::DefaultResidualVol(int i, Hermes::vector<std::string> areas,
        Hermes2DFunction<Scalar>* coeff,
        GeomType gt)
        : VectorFormVol<Scalar>(i), idx_i(i), coeff(coeff), gt(gt)
      {
        this->set_areas(areas);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes2DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultResidualVol<Scalar>::~DefaultResidualVol()
      {
        
        if(coeff == HERMES_ONE) delete coeff;
      };

      template<typename Scalar>
      Scalar DefaultResidualVol<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
            }
          }
        }
        return result;
      }

      template<typename Scalar>
      Ord DefaultResidualVol<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultResidualVol<Scalar>::clone() const
      {
        return new DefaultResidualVol(*this);
      }

      template<typename Scalar>
      DefaultResidualDiffusion<Scalar>::DefaultResidualDiffusion(int i, std::string area,
        Hermes1DFunction<Scalar>* coeff, GeomType gt)
        : VectorFormVol<Scalar>(i), idx_i(i), coeff(coeff), gt(gt)
      {
        this->set_area(area);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes1DFunction<Scalar>(1.0);
      };

      template<typename Scalar>
      DefaultResidualDiffusion<Scalar>::DefaultResidualDiffusion(int i, Hermes::vector<std::string> areas,
        Hermes1DFunction<Scalar>* coeff, GeomType gt)
        : VectorFormVol<Scalar>(i), idx_i(i), coeff(coeff), gt(gt)
      {
        this->set_areas(areas);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes1DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultResidualDiffusion<Scalar>::~DefaultResidualDiffusion()
      {
        
        if(coeff == HERMES_ONE) delete coeff;
      };

      template<typename Scalar>
      Scalar DefaultResidualDiffusion<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(u_ext[idx_i]->val[i])
              * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(u_ext[idx_i]->val[i])
                * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(u_ext[idx_i]->val[i])
                * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      Ord DefaultResidualDiffusion<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++) {
          result += wt[i] * coeff->value(u_ext[idx_i]->val[i])
            * (u_ext[idx_i]->dx[i] * v->dx[i] + u_ext[idx_i]->dy[i] * v->dy[i]);
        }
        if(gt != HERMES_PLANAR) result = result * Ord(1);

        return result;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultResidualDiffusion<Scalar>::clone() const
      {
        return new DefaultResidualDiffusion<Scalar>(*this);
      }

      template<typename Scalar>
      DefaultResidualAdvection<Scalar>::DefaultResidualAdvection(int i, std::string area,
        Hermes1DFunction<Scalar>* coeff1,
        Hermes1DFunction<Scalar>* coeff2,
        GeomType gt)
        : VectorFormVol<Scalar>(i), idx_i(i), coeff1(coeff1), coeff2(coeff2), gt(gt)
      {
        this->set_area(area);
        
        if(gt != HERMES_PLANAR) throw Hermes::Exceptions::Exception("Axisymmetric advection forms not implemented yet.");

        // If coeff1 == HERMES_ONE or coeff22 == HERMES_ONE, initialize it to be constant 1.0.
        if(coeff1 == HERMES_ONE) this->coeff1 = new Hermes1DFunction<Scalar>(1.0);
        if(coeff2 == HERMES_ONE) this->coeff2 = new Hermes1DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultResidualAdvection<Scalar>::DefaultResidualAdvection(int i, Hermes::vector<std::string> areas, \
        Hermes1DFunction<Scalar>* coeff1,
        Hermes1DFunction<Scalar>* coeff2,
        GeomType gt)
        : VectorFormVol<Scalar>(i),
        idx_i(i), coeff1(coeff1), coeff2(coeff2), gt(gt)
      {
        this->set_areas(areas);
        
        if(gt != HERMES_PLANAR) throw Hermes::Exceptions::Exception("Axisymmetric advection forms not implemented yet.");

        // If coeff1 == HERMES_ONE or coeff22 == HERMES_ONE, initialize it to be constant 1.0.
        if(coeff1 == HERMES_ONE) this->coeff1 = new Hermes1DFunction<Scalar>(1.0);
        if(coeff2 == HERMES_ONE) this->coeff2 = new Hermes1DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultResidualAdvection<Scalar>::~DefaultResidualAdvection()
      {
        
        //if(coeff1 != HERMES_ONE) delete coeff1;
        //if(coeff2 != HERMES_ONE) delete coeff2;
      };

      template<typename Scalar>
      Scalar DefaultResidualAdvection<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        Func<Scalar>* u_prev = u_ext[idx_i];
        for (int i = 0; i < n; i++) {
          result += wt[i] * (coeff1->value(u_prev->val[i]) * (u_prev->dx[i] * v->val[i])
            + coeff2->value(u_prev->val[i]) * (u_prev->dy[i] * v->val[i]));
        }
        return result;
      }

      template<typename Scalar>
      Ord DefaultResidualAdvection<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        Func<Ord>* u_prev = u_ext[idx_i];
        for (int i = 0; i < n; i++) {
          result += wt[i] * (coeff1->value(u_prev->val[i]) * (u_prev->dx[i] * v->val[i])
            + coeff2->value(u_prev->val[i]) * (u_prev->dy[i] * v->val[i]));
        }
        return result;
      }

      template<typename Scalar>
      VectorFormVol<Scalar>* DefaultResidualAdvection<Scalar>::clone() const
      {
        return new DefaultResidualAdvection<Scalar>(*this);
      }

      template<typename Scalar>
      DefaultMatrixFormSurf<Scalar>::DefaultMatrixFormSurf(int i, int j, std::string area,
        Hermes2DFunction<Scalar>* coeff,
        GeomType gt)
        : MatrixFormSurf<Scalar>(i, j), coeff(coeff), gt(gt)
      {
        this->set_area(area);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes2DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultMatrixFormSurf<Scalar>::DefaultMatrixFormSurf(int i, int j, Hermes::vector<std::string> areas,
        Hermes2DFunction<Scalar>* coeff,
        GeomType gt)
        : MatrixFormSurf<Scalar>(i, j), coeff(coeff), gt(gt)
      {
        this->set_areas(areas);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes2DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultMatrixFormSurf<Scalar>::~DefaultMatrixFormSurf()
      {
        
        if(coeff == HERMES_ONE) delete coeff;
      };

      template<typename Scalar>
      Scalar DefaultMatrixFormSurf<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      Ord DefaultMatrixFormSurf<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u->val[i] * v->val[i];
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      MatrixFormSurf<Scalar>* DefaultMatrixFormSurf<Scalar>::clone() const
      {
        return new DefaultMatrixFormSurf<Scalar>(*this);
      }

      template<typename Scalar>
      DefaultJacobianFormSurf<Scalar>::DefaultJacobianFormSurf(int i, int j, std::string area,
        Hermes1DFunction<Scalar>* coeff,
        GeomType gt)
        : MatrixFormSurf<Scalar>(i, j),
        idx_j(j), coeff(coeff), gt(gt)
      {
        this->set_area(area);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes1DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultJacobianFormSurf<Scalar>::DefaultJacobianFormSurf(int i, int j, Hermes::vector<std::string> areas,
        Hermes1DFunction<Scalar>* coeff,
        GeomType gt)
        : MatrixFormSurf<Scalar>(i, j), coeff(coeff), gt(gt)
      {
        this->set_areas(areas);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes1DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultJacobianFormSurf<Scalar>::~DefaultJacobianFormSurf()
      {
        
        if(coeff == HERMES_ONE) delete coeff;
      };

      template<typename Scalar>
      Scalar DefaultJacobianFormSurf<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *u, Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        for (int i = 0; i < n; i++) {
          result += wt[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u_ext[idx_j]->val[i]
          + coeff->value(u_ext[idx_j]->val[i]))
            * u->val[i] * v->val[i];
        }
        return result;
      }

      template<typename Scalar>
      Ord DefaultJacobianFormSurf<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *u,
        Func<Ord> *v, Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        for (int i = 0; i < n; i++) {
          result += wt[i] * (coeff->derivative(u_ext[idx_j]->val[i]) * u_ext[idx_j]->val[i]
          + coeff->value(u_ext[idx_j]->val[i]))
            * u->val[i] * v->val[i];
        }
        return result;
      }

      template<typename Scalar>
      MatrixFormSurf<Scalar>* DefaultJacobianFormSurf<Scalar>::clone() const
      {
        return new DefaultJacobianFormSurf<Scalar>(*this);
      }

      template<typename Scalar>
      DefaultVectorFormSurf<Scalar>::DefaultVectorFormSurf(int i, std::string area,
        Hermes2DFunction<Scalar>* coeff,
        GeomType gt)
        : VectorFormSurf<Scalar>(i), coeff(coeff), gt(gt)
      {
        this->set_area(area);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes2DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultVectorFormSurf<Scalar>::DefaultVectorFormSurf(int i, Hermes::vector<std::string> areas,
        Hermes2DFunction<Scalar>* coeff,
        GeomType gt)
        : VectorFormSurf<Scalar>(i), coeff(coeff), gt(gt)
      {
        this->set_areas(areas);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes2DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultVectorFormSurf<Scalar>::~DefaultVectorFormSurf()
      {
        
        if(coeff == HERMES_ONE) delete coeff;
      };

      template<typename Scalar>
      Scalar DefaultVectorFormSurf<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      Ord DefaultVectorFormSurf<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * v->val[i];
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      VectorFormSurf<Scalar>* DefaultVectorFormSurf<Scalar>::clone() const
      {
        return new DefaultVectorFormSurf<Scalar>(*this);
      }

      template<typename Scalar>
      DefaultResidualSurf<Scalar>::DefaultResidualSurf(int i, std::string area,
        Hermes2DFunction<Scalar>* coeff,
        GeomType gt)
        : VectorFormSurf<Scalar>(i), idx_i(i), coeff(coeff), gt(gt)
      {
        this->set_area(area);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes2DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultResidualSurf<Scalar>::DefaultResidualSurf(int i, Hermes::vector<std::string> areas,
        Hermes2DFunction<Scalar>* coeff,
        GeomType gt)
        : VectorFormSurf<Scalar>(i), idx_i(i), coeff(coeff), gt(gt)
      {
        this->set_areas(areas);
        
        // If coeff is HERMES_ONE, initialize it to be constant 1.0.
        if(coeff == HERMES_ONE) this->coeff = new Hermes2DFunction<Scalar>(1.0);
      }

      template<typename Scalar>
      DefaultResidualSurf<Scalar>::~DefaultResidualSurf()
      {
        
        if(coeff == HERMES_ONE) delete coeff;
      };

      template<typename Scalar>
      Scalar DefaultResidualSurf<Scalar>::value(int n, double *wt, Func<Scalar> *u_ext[], Func<double> *v,
        Geom<double> *e, Func<Scalar> **ext) const
      {
        Scalar result = 0;
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      Ord DefaultResidualSurf<Scalar>::ord(int n, double *wt, Func<Ord> *u_ext[], Func<Ord> *v,
        Geom<Ord> *e, Func<Ord> **ext) const
      {
        Ord result = Ord(0);
        if(gt == HERMES_PLANAR) {
          for (int i = 0; i < n; i++) {
            result += wt[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
          }
        }
        else {
          if(gt == HERMES_AXISYM_X) {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->y[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
            }
          }
          else {
            for (int i = 0; i < n; i++) {
              result += wt[i] * e->x[i] * coeff->value(e->x[i], e->y[i]) * u_ext[idx_i]->val[i] * v->val[i];
            }
          }
        }

        return result;
      }

      template<typename Scalar>
      VectorFormSurf<Scalar>* DefaultResidualSurf<Scalar>::clone() const
      {
        return new DefaultResidualSurf(*this);
      }

      template<typename Scalar>
      DefaultWeakFormLaplace<Scalar>::DefaultWeakFormLaplace(std::string area,
        Hermes1DFunction<Scalar>* coeff,
        GeomType gt) : WeakForm<Scalar>()
      {
        // Jacobian.
        this->add_matrix_form(new DefaultJacobianDiffusion<Scalar>(0, 0, area, coeff, HERMES_NONSYM, gt));

        // Residual.
        this->add_vector_form(new DefaultResidualDiffusion<Scalar>(0, area, coeff, gt));
      };

      template<typename Scalar>
      DefaultWeakFormPoisson<Scalar>::DefaultWeakFormPoisson() : WeakForm<Scalar>()
      {
      };

      template<typename Scalar>
      DefaultWeakFormPoisson<Scalar>::DefaultWeakFormPoisson(std::string area,
        Hermes1DFunction<Scalar>* coeff,
        Hermes2DFunction<Scalar>* f,
        GeomType gt) : WeakForm<Scalar>()
      {
        // Jacobian.
        // NOTE: The flag HERMES_NONSYM is important here.
        this->add_matrix_form(new DefaultJacobianDiffusion<Scalar>(0, 0, area, coeff, HERMES_NONSYM, gt));

        // Residual.
        this->add_vector_form(new DefaultResidualDiffusion<Scalar>(0, area, coeff, gt));
        this->add_vector_form(new DefaultVectorFormVol<Scalar>(0, area, f, gt));
      };

      template class HERMES_API DefaultMatrixFormVol<double>;
      template class HERMES_API DefaultMatrixFormVol<std::complex<double> >;
      template class HERMES_API DefaultJacobianDiffusion<double>;
      template class HERMES_API DefaultJacobianDiffusion<std::complex<double> >;
      template class HERMES_API DefaultMatrixFormDiffusion<double>;
      template class HERMES_API DefaultMatrixFormDiffusion<std::complex<double> >;
      template class HERMES_API DefaultResidualAdvection<double>;
      template class HERMES_API DefaultResidualAdvection<std::complex<double> >;
      template class HERMES_API DefaultJacobianAdvection<double>;
      template class HERMES_API DefaultJacobianAdvection<std::complex<double> >;
      template class HERMES_API DefaultResidualDiffusion<double>;
      template class HERMES_API DefaultResidualDiffusion<std::complex<double> >;
      template class HERMES_API DefaultMatrixFormSurf<double>;
      template class HERMES_API DefaultMatrixFormSurf<std::complex<double> >;
      template class HERMES_API DefaultVectorFormSurf<double>;
      template class HERMES_API DefaultVectorFormSurf<std::complex<double> >;
      template class HERMES_API DefaultVectorFormVol<double>;
      template class HERMES_API DefaultVectorFormVol<std::complex<double> >;
      template class HERMES_API DefaultJacobianFormSurf<double>;
      template class HERMES_API DefaultJacobianFormSurf<std::complex<double> >;
      template class HERMES_API DefaultWeakFormLaplace<double>;
      template class HERMES_API DefaultWeakFormLaplace<std::complex<double> >;
      template class HERMES_API DefaultWeakFormPoisson<double>;
      template class HERMES_API DefaultWeakFormPoisson<std::complex<double> >;
      template class HERMES_API DefaultResidualSurf<double>;
      template class HERMES_API DefaultResidualSurf<std::complex<double> >;
      template class HERMES_API DefaultResidualVol<double>;
      template class HERMES_API DefaultResidualVol<std::complex<double> >;
    };
  }
}
