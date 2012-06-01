// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "global.h"
#include "quadrature/quad_all.h"
#include "mesh.h"
#include "traverse.h"
#include "refmap.h"
#include "function/solution.h"
#include "quadrature/limit_order.h"
#include "integrals/h1.h"
#include "discrete_problem.h"

namespace Hermes
{
  namespace Hermes2D
  {
    class Transformable;

    Hermes::Hermes2D::Api HermesApi;

    Api::Parameter::Parameter(int defaultVal)
    {
      this->defaultVal = defaultVal;
      this->userSet = false;
    }

    Api::Api()
    {
      this->parameters.insert(std::pair<std::string, Parameter*> ("num_threads",new Parameter(NUM_THREADS)));
    }

    Api::~Api()
    {
      this->parameters.clear();
    }

    int Api::getParamValue(std::string param)
    {
      if(this->parameters.find(param) == parameters.end())
        throw new Hermes::Exceptions::ValueException("parameter name", param);
      if(this->parameters.find(param)->second->userSet)
        return this->parameters.find(param)->second->userVal;
      else
        return this->parameters.find(param)->second->defaultVal;
    }

    void Api::setParamValue(std::string param, int value)
    {
      if(this->parameters.find(param) == parameters.end())
        throw new Hermes::Exceptions::ValueException("parameter name", param);
      this->parameters.find(param)->second->userSet = true;
      this->parameters.find(param)->second->userVal = value;
    }

    template<typename Scalar>
    double Global<Scalar>::get_l2_norm(Vector<Scalar>* vec)
    {
      _F_;
      Scalar val = 0;
      for (unsigned int i = 0; i < vec->length(); i++)
      {
        Scalar inc = vec->get(i);
        val = val + inc*conj(inc);
      }
      return sqrt(std::abs(val));
    }

    template<typename Scalar>
    double Global<Scalar>::calc_abs_error(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, int norm_type)
    {
      // sanity checks
      if (sln1 == NULL) error("sln1 is NULL in calc_abs_error().");
      if (sln2 == NULL) error("sln2 is NULL in calc_abs_error().");

      Quad2D* quad = &g_quad_2d_std;
      sln1->set_quad_2d(quad);
      sln2->set_quad_2d(quad);

      Mesh* meshes[2] = { sln1->get_mesh(), sln2->get_mesh() };
      Transformable* tr[2] = { sln1, sln2 };
      Traverse trav(true);
      trav.begin(2, meshes, tr);

      double error = 0.0;
      Traverse::State* ee;
      while ((ee = trav.get_next_state()) != NULL)
      {
        update_limit_table(ee->e[0]->get_mode());

        RefMap* ru = sln1->get_refmap();
        RefMap* rv = sln2->get_refmap();
        switch (norm_type)
        {
        case HERMES_L2_NORM:
          error += error_fn_l2(sln1, sln2, ru, rv);
          break;
        case HERMES_H1_NORM:
          error += error_fn_h1(sln1, sln2, ru, rv);
          break;
        case HERMES_HCURL_NORM:
          error += error_fn_hc(sln1, sln2, ru, rv);
          break;
        case HERMES_HDIV_NORM:
          error += error_fn_hdiv(sln1, sln2, ru, rv);
          break;
        default: error("Unknown norm in calc_error().");
        }
      }
      trav.finish();
      return sqrt(error);
    }

    template<typename Scalar>
    double Global<Scalar>::calc_rel_error(MeshFunction<Scalar>* sln, MeshFunction<Scalar>* ref_sln, int norm_type)
    {
      double error = calc_abs_error(sln, ref_sln, norm_type);
      double norm = calc_norm(ref_sln, norm_type);

      return error/norm;
    }

    template<typename Scalar>
    double Global<Scalar>::calc_norm(MeshFunction<Scalar>* sln, int norm_type)
    {
      Quad2D* quad = &g_quad_2d_std;
      sln->set_quad_2d(quad);

      double norm = 0.0;
      Element* e;
      Mesh* mesh = sln->get_mesh();

      for_all_active_elements(e, mesh)
      {
        // set maximum integration order for use in integrals, see limit_order()
        update_limit_table(e->get_mode());

        sln->set_active_element(e);
        RefMap* ru = sln->get_refmap();

        switch (norm_type)
        {
        case HERMES_L2_NORM:
          norm += norm_fn_l2(sln, ru);
          break;
        case HERMES_H1_NORM:
          norm += norm_fn_h1(sln, ru);
          break;
        case HERMES_HCURL_NORM:
          norm += norm_fn_hc(sln, ru);
          break;
        case HERMES_HDIV_NORM:
          norm += norm_fn_hdiv(sln, ru);
          break;
        default: error("Unknown norm in calc_norm().");
        }
      }
      return sqrt(norm);
    }

    template<typename Scalar>
    double Global<Scalar>::calc_norms(Hermes::vector<Solution<Scalar>*> slns)
    {
      // Calculate norms for all solutions.
      Hermes::vector<double> norms;
      int n = slns.size();
      for (int i = 0; i<n; i++)
      {
        switch (slns[i]->get_space_type())
        {
        case HERMES_H1_SPACE: norms.push_back(calc_norm(slns[i], HERMES_H1_NORM)); break;
        case HERMES_HCURL_SPACE: norms.push_back(calc_norm(slns[i], HERMES_HCURL_NORM)); break;
        case HERMES_HDIV_SPACE: norms.push_back(calc_norm(slns[i], HERMES_HDIV_NORM)); break;
        case HERMES_L2_SPACE: norms.push_back(calc_norm(slns[i], HERMES_L2_NORM)); break;
        default: error("Internal in calc_norms(): unknown space type.");
        }
      }
      // Calculate the resulting norm.
      double result = 0;
      for (int i = 0; i < n; i++)
        result += norms[i] * norms[i];
      return sqrt(result);
    }

    template<typename Scalar>
    double Global<Scalar>::calc_abs_errors(Hermes::vector<Solution<Scalar>*> slns1, Hermes::vector<Solution<Scalar>*> slns2)
    {
      // Calculate errors for all solutions.
      Hermes::vector<double> errors;
      int n = slns1.size();
      for (int i = 0; i < n; i++)
      {
        switch (slns1[i]->get_space_type())
        {
        case HERMES_H1_SPACE: errors.push_back(calc_abs_error(slns1[i], slns2[i], HERMES_H1_NORM)); break;
        case HERMES_HCURL_SPACE: errors.push_back(calc_abs_error(slns1[i], slns2[i], HERMES_HCURL_NORM)); break;
        case HERMES_HDIV_SPACE: errors.push_back(calc_abs_error(slns1[i], slns2[i], HERMES_HDIV_NORM)); break;
        case HERMES_L2_SPACE: errors.push_back(calc_abs_error(slns1[i], slns2[i], HERMES_L2_NORM)); break;
        default: error("Internal in calc_norms(): unknown space type.");
        }
      }
      // Calculate the resulting error.
      double result = 0;
      for (int i = 0; i < n; i++)
        result += errors[i] * errors[i];
      return sqrt(result);
    }

    template<typename Scalar>
    double Global<Scalar>::calc_rel_errors(Hermes::vector<Solution<Scalar>*> slns1, Hermes::vector<Solution<Scalar>*> slns2)
    {
      return calc_abs_errors(slns1, slns2) / calc_norms(slns2);
    }

    template<typename Scalar>
    double Global<Scalar>::error_fn_h1(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv)
    {
      Quad2D* quad = sln1->get_quad_2d();

      int o = 2*std::max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
      limit_order_nowarn(o, ru->get_active_element()->get_mode());

      sln1->set_quad_order(o);
      sln2->set_quad_order(o);

      Scalar *uval, *vval, *dudx, *dudy, *dvdx, *dvdy;
      uval = sln1->get_fn_values();
      vval = sln2->get_fn_values();
      sln1->get_dx_dy_values(dudx, dudy);
      sln2->get_dx_dy_values(dvdx, dvdy);

      double result = 0.0;
      h1_integrate_expression(Hermes::sqr(uval[i] - vval[i]) +
        sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]));
      return result;
    }

    template<typename Scalar>
    double Global<Scalar>::norm_fn_h1(MeshFunction<Scalar>* sln, RefMap* ru)
    {
      Quad2D* quad = sln->get_quad_2d();

      int o = 2 * sln->get_fn_order() + ru->get_inv_ref_order();
      limit_order_nowarn(o, ru->get_active_element()->get_mode());

      sln->set_quad_order(o);

      Scalar *uval, *dudx, *dudy;
      uval = sln->get_fn_values();
      sln->get_dx_dy_values(dudx, dudy);

      double result = 0.0;
      h1_integrate_expression(sqr(uval[i]) + sqr(dudx[i]) + sqr(dudy[i]));
      return result;
    }

    template<typename Scalar>
    double Global<Scalar>::error_fn_l2(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv)
    {
      Quad2D* quad = sln1->get_quad_2d();

      int o = 2*std::max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
      limit_order_nowarn(o, ru->get_active_element()->get_mode());

      sln1->set_quad_order(o, H2D_FN_VAL);
      sln2->set_quad_order(o, H2D_FN_VAL);

      Scalar *uval, *vval;
      uval = sln1->get_fn_values();
      vval = sln2->get_fn_values();

      double result = 0.0;
      h1_integrate_expression(Hermes::sqr(uval[i] - vval[i]));
      return result;
    }

    template<typename Scalar>
    double Global<Scalar>::norm_fn_l2(MeshFunction<Scalar>* sln, RefMap* ru)
    {
      Quad2D* quad = sln->get_quad_2d();

      int o = 2 *sln->get_fn_order() + ru->get_inv_ref_order();
      limit_order_nowarn(o, ru->get_active_element()->get_mode());

      sln->set_quad_order(o, H2D_FN_VAL);

      Scalar* uval = sln->get_fn_values();

      double result = 0.0;
      h1_integrate_expression(Hermes::sqr(uval[i]));
      return result;
    }

    template<typename Scalar>
    double Global<Scalar>::error_fn_hc(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv)
    {
      Quad2D* quad = sln1->get_quad_2d();

      int o = 2 * std::max(sln1->get_fn_order(), sln2->get_fn_order()) + 2 + ru->get_inv_ref_order();
      limit_order_nowarn(o, ru->get_active_element()->get_mode());

      sln1->set_quad_order(o);
      sln2->set_quad_order(o);


      Scalar *uval0 = sln1->get_fn_values(0), *uval1 = sln1->get_fn_values(1);
      Scalar *udx1  = sln1->get_dx_values(1), *udy0  = sln1->get_dy_values(0);
      Scalar *vval0 = sln2->get_fn_values(0), *vval1 = sln2->get_fn_values(1);
      Scalar *vdx1  = sln2->get_dx_values(1), *vdy0  = sln2->get_dy_values(0);

      double result = 0.0;
      h1_integrate_expression(Hermes::sqr(uval0[i] - vval0[i]) + Hermes::sqr(uval1[i] - vval1[i]) +
        sqr((udx1[i] - udy0[i]) - (vdx1[i] - vdy0[i])));
      return result;
    }

    template<typename Scalar>
    double Global<Scalar>::norm_fn_hc(MeshFunction<Scalar>* sln, RefMap* ru)
    {
      Quad2D* quad = sln->get_quad_2d();

      int o = 2 * sln->get_fn_order() + 2 + ru->get_inv_ref_order();
      limit_order_nowarn(o, ru->get_active_element()->get_mode());

      sln->set_quad_order(o);

      Scalar *uval0 = sln->get_fn_values(0), *uval1 = sln->get_fn_values(1);
      Scalar *udx1  = sln->get_dx_values(1), *udy0  = sln->get_dy_values(0);

      double result = 0.0;
      h1_integrate_expression(Hermes::sqr(uval0[i]) + Hermes::sqr(uval1[i]) + Hermes::sqr(udx1[i] - udy0[i]));
      return result;
    }

    template<typename Scalar>
    double Global<Scalar>::error_fn_hcl2(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv)
    {
      Quad2D* quad = sln1->get_quad_2d();

      int o = 2 * std::max(sln1->get_fn_order(), sln2->get_fn_order()) + 2 + ru->get_inv_ref_order();
      limit_order_nowarn(o, ru->get_active_element()->get_mode());

      sln1->set_quad_order(o);
      sln2->set_quad_order(o);


      Scalar *uval0 = sln1->get_fn_values(0), *uval1 = sln1->get_fn_values(1);
      Scalar *vval0 = sln2->get_fn_values(0), *vval1 = sln2->get_fn_values(1);

      double result = 0.0;
      h1_integrate_expression(Hermes::sqr(uval0[i] - vval0[i]) + Hermes::sqr(uval1[i] - vval1[i]));
      return result;
    }

    template<typename Scalar>
    double Global<Scalar>::norm_fn_hcl2(MeshFunction<Scalar>* sln, RefMap* ru)
    {
      Quad2D* quad = sln->get_quad_2d();

      int o = 2 * sln->get_fn_order() + 2 + ru->get_inv_ref_order();
      limit_order_nowarn(o, ru->get_active_element()->get_mode());

      sln->set_quad_order(o);

      Scalar *uval0 = sln->get_fn_values(0), *uval1 = sln->get_fn_values(1);
      Scalar *udx1  = sln->get_dx_values(1), *udy0  = sln->get_dy_values(0);

      double result = 0.0;
      h1_integrate_expression(Hermes::sqr(uval0[i]) + Hermes::sqr(uval1[i]));
      return result;
    }

    template<typename Scalar>
    double Global<Scalar>::error_fn_hdiv(MeshFunction<Scalar>* sln1, MeshFunction<Scalar>* sln2, RefMap* ru, RefMap* rv)
    {
      error("error_fn_hdiv() not implemented yet.");

      // Hcurl code
      Quad2D* quad = sln1->get_quad_2d();

      int o = 2 * std::max(sln1->get_fn_order(), sln2->get_fn_order()) + 2 + ru->get_inv_ref_order();
      limit_order_nowarn(o, ru->get_active_element()->get_mode());

      sln1->set_quad_order(o);
      sln2->set_quad_order(o);


      Scalar *uval0 = sln1->get_fn_values(0), *uval1 = sln1->get_fn_values(1);
      Scalar *udx1  = sln1->get_dx_values(1), *udy0  = sln1->get_dy_values(0);
      Scalar *vval0 = sln2->get_fn_values(0), *vval1 = sln2->get_fn_values(1);
      Scalar *vdx1  = sln2->get_dx_values(1), *vdy0  = sln2->get_dy_values(0);

      double result = 0.0;
      h1_integrate_expression(Hermes::sqr(uval0[i] - vval0[i]) + Hermes::sqr(uval1[i] - vval1[i]) +
        Hermes::sqr((udx1[i] - udy0[i]) - (vdx1[i] - vdy0[i])));
      return result;
    }

    template<typename Scalar>
    double Global<Scalar>::norm_fn_hdiv(MeshFunction<Scalar>* sln, RefMap* ru)
    {
      error("norm_fn_hdiv() not implemented yet.");

      // Hcurl code
      Quad2D* quad = sln->get_quad_2d();

      int o = 2 * sln->get_fn_order() + 2 + ru->get_inv_ref_order();
      limit_order_nowarn(o, ru->get_active_element()->get_mode());

      sln->set_quad_order(o);

      Scalar *uval0 = sln->get_fn_values(0), *uval1 = sln->get_fn_values(1);
      Scalar *udx1  = sln->get_dx_values(1), *udy0  = sln->get_dy_values(0);

      double result = 0.0;
      h1_integrate_expression(Hermes::sqr(uval0[i]) + Hermes::sqr(uval1[i]) + Hermes::sqr(udx1[i] - udy0[i]));
      return result;
    }

    HERMES_API void throw_exception(char *text)
    {
      throw std::runtime_error(text);
    }

    template class HERMES_API Global<double>;
    template class HERMES_API Global<std::complex<double> >;
  }
}