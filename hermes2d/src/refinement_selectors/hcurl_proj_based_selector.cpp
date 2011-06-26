#include "hermes2d_common_defs.h"
#include "matrix.h"
#include "solution.h"
#include "shapeset/shapeset_hc_all.h"
#include "element_to_refine.h"
#include "hcurl_proj_based_selector.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors 
    {
      HcurlShapeset HcurlProjBasedSelector::default_shapeset;

      const int HcurlProjBasedSelector::H2DRS_MAX_HCURL_ORDER = 6;

      HcurlProjBasedSelector::HcurlProjBasedSelector(CandList cand_list, double conv_exp, int max_order, HcurlShapeset* user_shapeset)
        : ProjBasedSelector<std::complex<double> >(cand_list, conv_exp, max_order, user_shapeset == NULL ? &default_shapeset : user_shapeset, Range<int>(), Range<int>(0, H2DRS_MAX_HCURL_ORDER))
        , precalc_rvals_curl(NULL) {}

      HcurlProjBasedSelector::~HcurlProjBasedSelector() 
      {
        delete[] precalc_rvals_curl;
      }

      void HcurlProjBasedSelector::set_current_order_range(Element* element) 
      {
        current_max_order = this->max_order;
        if (current_max_order == H2DRS_DEFAULT_ORDER)
          current_max_order = std::min(H2DRS_MAX_HCURL_ORDER, (20 - element->iro_cache)/2 - 1); // default
        else
          current_max_order = std::min(max_order, (20 - element->iro_cache)/2 - 1); // user specified
        current_min_order = 0;
      }

      void HcurlProjBasedSelector::precalc_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const Hermes::vector<ShapeInx>& shapes, const int max_shape_inx, TrfShape& svals) 
      {
        //for all transformations
        bool done = false;
        int inx_trf = 0;
        while (!done && inx_trf < H2D_TRF_NUM) 
        {
          //prepare data for processing
          const Trf& trf = trfs[inx_trf];
          Hermes::vector<TrfShapeExp>& trf_svals = svals[inx_trf];

          //allocate
          trf_svals.resize(max_shape_inx + 1);

          //for all shapes
          const int num_shapes = (int)shapes.size();
          for(int i = 0; i < num_shapes; i++) 
          {
            int inx_shape = shapes[i].inx;
            TrfShapeExp& shape_exp = trf_svals[inx_shape];

            //allocate
            shape_exp.allocate(H2D_HCFE_NUM, num_gip_points);

            //for all GIP points
            for(int k = 0; k < num_gip_points; k++) 
            {
              //transform coordinates
              double ref_x = gip_points[k][H2D_GIP2D_X] * trf.m[0] + trf.t[0];
              double ref_y = gip_points[k][H2D_GIP2D_Y] * trf.m[1] + trf.t[1];

              //for all expansions: retrieve values
              shape_exp[H2D_HCFE_VALUE0][k] = shapeset->get_fn_value(inx_shape, ref_x, ref_y, 0);
              shape_exp[H2D_HCFE_VALUE1][k] = shapeset->get_fn_value(inx_shape, ref_x, ref_y, 1);
              shape_exp[H2D_HCFE_CURL][k] = shapeset->get_dx_value(inx_shape, ref_x, ref_y, 1) - shapeset->get_dy_value(inx_shape, ref_x, ref_y, 0);
            }
          }

          //move to the next transformation
          if (inx_trf == H2D_TRF_IDENTITY)
            done = true;
          else 
          {
            inx_trf++;
            if (inx_trf >= num_noni_trfs) //if all transformations were processed, move to the identity transformation
              inx_trf = H2D_TRF_IDENTITY;
          }
        }
        error_if(!done, "All transformation processed but identity transformation not found."); //identity transformation has to be the last transformation
      }

      void HcurlProjBasedSelector::precalc_ortho_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const Hermes::vector<ShapeInx>& shapes, const int max_shape_inx, TrfShape& svals) 
      {
        //calculate values
        precalc_shapes(gip_points, num_gip_points, trfs, num_noni_trfs, shapes, max_shape_inx, svals);

        //calculate orthonormal basis
        const int num_shapes = (int)shapes.size();
        for(int i = 0; i < num_shapes; i++) 
        {
          const int inx_shape_i = shapes[i].inx;

          //orthogonalize
          for(int j = 0; j < i; j++) 
          {
            const int inx_shape_j = shapes[j].inx;

            //calculate product of non-transformed functions
            double product = 0.0;
            for(int k = 0; k < num_gip_points; k++) 
            {
              double sum = 0.0;
              sum += svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_HCFE_VALUE0][k] * svals[H2D_TRF_IDENTITY][inx_shape_j][H2D_HCFE_VALUE0][k];
              sum += svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_HCFE_VALUE1][k] * svals[H2D_TRF_IDENTITY][inx_shape_j][H2D_HCFE_VALUE1][k];
              sum += svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_HCFE_CURL][k] * svals[H2D_TRF_IDENTITY][inx_shape_j][H2D_HCFE_CURL][k];
              product += gip_points[k][H2D_GIP2D_W] * sum;
            }

            //for all transformations
            int inx_trf = 0;
            bool done = false;
            while (!done && inx_trf < H2D_TRF_NUM) 
            {
              //for all integration points
              for(int k = 0; k < num_gip_points; k++) 
              {
                svals[inx_trf][inx_shape_i][H2D_HCFE_VALUE0][k] -= product * svals[inx_trf][inx_shape_j][H2D_HCFE_VALUE0][k];
                svals[inx_trf][inx_shape_i][H2D_HCFE_VALUE1][k] -= product * svals[inx_trf][inx_shape_j][H2D_HCFE_VALUE1][k];
                svals[inx_trf][inx_shape_i][H2D_HCFE_CURL][k] -= product * svals[inx_trf][inx_shape_j][H2D_HCFE_CURL][k];
              }

              //move to the next transformation
              if (inx_trf == H2D_TRF_IDENTITY)
                done = true;
              else 
              {
                inx_trf++;
                if (inx_trf >= num_noni_trfs) //if all transformations were processed, move to the identity transformation
                  inx_trf = H2D_TRF_IDENTITY;
              }
            }
            error_if(!done, "All transformation processed but identity transformation not found."); //identity transformation has to be the last transformation
          }

          //normalize
          //calculate norm
          double norm_squared = 0.0;
          for(int k = 0; k < num_gip_points; k++) 
          {
            double sum = 0.0;
            sum += sqr(svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_HCFE_VALUE0][k]);
            sum += sqr(svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_HCFE_VALUE1][k]);
            sum += sqr(svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_HCFE_CURL][k]);
            norm_squared += gip_points[k][H2D_GIP2D_W] * sum;
          }
          double norm = sqrt(norm_squared);
          assert_msg(finite(1/norm), "Norm (%g) is almost zero.", norm);

          //for all transformations: normalize
          int inx_trf = 0;
          bool done = false;
          while (!done && inx_trf < H2D_TRF_NUM) 
          {
            //for all integration points
            for(int k = 0; k < num_gip_points; k++) 
            {
              svals[inx_trf][inx_shape_i][H2D_HCFE_VALUE0][k] /= norm;
              svals[inx_trf][inx_shape_i][H2D_HCFE_VALUE1][k] /= norm;
              svals[inx_trf][inx_shape_i][H2D_HCFE_CURL][k] /= norm;
            }

            //move to the next transformation
            if (inx_trf == H2D_TRF_IDENTITY)
              done = true;
            else 
            {
              inx_trf++;
              if (inx_trf >= num_noni_trfs) //if all transformations were processed, move to the identity transformation
                inx_trf = H2D_TRF_IDENTITY;
            }
          }
          error_if(!done, "All transformation processed but identity transformation not found."); //identity transformation has to be the last transformation
        }
      }

      std::complex<double>** HcurlProjBasedSelector::precalc_ref_solution(int inx_son, Solution<std::complex<double> >* rsln, Element* element, int intr_gip_order) 
      {
        //set element and integration order
        rsln->set_active_element(element);
        rsln->set_quad_order(intr_gip_order);
        const int num_gip = rsln->get_quad_2d()->get_num_points(intr_gip_order);

        //allocate space for Curl
        if (precalc_rvals_curl == NULL)
          precalc_rvals_curl = new_matrix<std::complex<double> >(H2D_MAX_ELEMENT_SONS, num_gip);

        //prepre for curl
        std::complex<double>* curl = precalc_rvals_curl[inx_son];
        std::complex<double>* d1dx = rsln->get_dx_values(1);
        std::complex<double>* d0dy = rsln->get_dy_values(0);
        for(int i = 0; i < num_gip; i++)
          curl[i] = d1dx[i] - d0dy[i];

        //fill with values
        std::complex<double>** rvals_son = precalc_rvals[inx_son];
        rvals_son[H2D_HCFE_VALUE0] = rsln->get_fn_values(0);
        rvals_son[H2D_HCFE_VALUE1] = rsln->get_fn_values(1);
        rvals_son[H2D_HCFE_CURL] = curl;

        return rvals_son;
      }

      double** HcurlProjBasedSelector::build_projection_matrix(double3* gip_points, int num_gip_points,
        const int* shape_inx, const int num_shapes) 
      {
        //allocate
        double** matrix = new_matrix<double>(num_shapes, num_shapes);

        //calculate products
        int inx_row = 0;
        for(int i = 0; i < num_shapes; i++, inx_row += num_shapes) 
        {
          double* matrix_row = matrix[i];
          int shape0_inx = shape_inx[i];
          for(int k = 0; k < num_shapes; k++) 
          {
            int shape1_inx = shape_inx[k];

            double value = 0.0;
            for(int j = 0; j < num_gip_points; j++) 
            {
              double gip_x = gip_points[j][H2D_GIP2D_X], gip_y = gip_points[j][H2D_GIP2D_Y];
              double value0[2] = { shapeset->get_value(H2D_FEI_VALUE, shape0_inx, gip_x, gip_y, 0), shapeset->get_value(H2D_FEI_VALUE, shape0_inx, gip_x, gip_y, 1) };
              double value1[2] = { shapeset->get_value(H2D_FEI_VALUE, shape1_inx, gip_x, gip_y, 0), shapeset->get_value(H2D_FEI_VALUE, shape1_inx, gip_x, gip_y, 1) };
              double d1dx0 = shapeset->get_value(H2D_FEI_DX, shape0_inx, gip_x, gip_y, 1);
              double d1dx1 = shapeset->get_value(H2D_FEI_DX, shape1_inx, gip_x, gip_y, 1);
              double d0dy0 = shapeset->get_value(H2D_FEI_DY, shape0_inx, gip_x, gip_y, 0);
              double d0dy1 = shapeset->get_value(H2D_FEI_DY, shape1_inx, gip_x, gip_y, 0);
              double curl0 = d1dx0 - d0dy0;
              double curl1 = d1dx1 - d0dy1;

              value += gip_points[j][H2D_GIP2D_W] * (value0[0]*value1[0] + value0[1]*value1[1] + curl0*curl1);
            }

            matrix_row[k] = value;
          }
        }

        return matrix;
      }

      std::complex<double> HcurlProjBasedSelector::evaluate_rhs_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemSubShapeFunc& sub_shape) 
      {
        double coef_curl = std::abs(sub_trf.coef_mx * sub_trf.coef_my);
        std::complex<double> total_value = 0;
        for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++) 
        {
          //get location and transform it
          double3 &gip_pt = sub_gip.gip_points[gip_inx];

          //get value of a shape function
          std::complex<double> shape_value0 = sub_shape.svals[H2D_HCFE_VALUE0][gip_inx];
          std::complex<double> shape_value1 = sub_shape.svals[H2D_HCFE_VALUE1][gip_inx];
          std::complex<double> shape_curl = sub_shape.svals[H2D_HCFE_CURL][gip_inx];

          ////DEBUG-BEGIN
          //double ref_x = gip_pt[H2D_GIP2D_X] * sub_trf.trf->m[0] + sub_trf.trf->t[0];
          //double ref_y = gip_pt[H2D_GIP2D_Y] * sub_trf.trf->m[1] + sub_trf.trf->t[1];
          //std::complex<double> shape_value0A = shapeset->get_fn_value(sub_shape.inx, ref_x, ref_y, 0);
          //std::complex<double> shape_value1A = shapeset->get_fn_value(sub_shape.inx, ref_x, ref_y, 1);
          //std::complex<double> shape_curlA = shapeset->get_dx_value(sub_shape.inx, ref_x, ref_y, 1) - shapeset->get_dy_value(sub_shape.inx, ref_x, ref_y, 0);
          //error_if(std::abs(shape_value0 - shape_value0A) > 1E-15
          //  || std::abs(shape_value1 - shape_value1A) > 1E-15
          //  || std::abs(shape_curl - shape_curlA) > 1E-15, "A1");
          ////DEBUG-END

          //get value of ref. solution
          std::complex<double> ref_value0 = sub_trf.coef_mx * sub_gip.rvals[H2D_HCFE_VALUE0][gip_inx];
          std::complex<double> ref_value1 = sub_trf.coef_my * sub_gip.rvals[H2D_HCFE_VALUE1][gip_inx];
          std::complex<double> ref_curl = coef_curl * sub_gip.rvals[H2D_HCFE_CURL][gip_inx]; //coef_curl * curl

          //evaluate a right-hand value
          std::complex<double> value = (shape_value0 * ref_value0)
            + (shape_value1 * ref_value1)
            + (shape_curl * ref_curl);

          total_value += gip_pt[H2D_GIP2D_W] * value;
        }
        return total_value;
      }

      double HcurlProjBasedSelector::evaluate_error_squared_subdomain(Element* sub_elem, const ElemGIP& sub_gip, const ElemSubTrf& sub_trf, const ElemProj& elem_proj) 
      {
        double total_error_squared = 0;
        double coef_curl = std::abs(sub_trf.coef_mx * sub_trf.coef_my);
        for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++) 
        {
          //get location and transform it
          double3 &gip_pt = sub_gip.gip_points[gip_inx];

          //calculate value of projected solution
          std::complex<double> proj_value0 = 0, proj_value1 = 0, proj_curl = 0;
          for(int i = 0; i < elem_proj.num_shapes; i++) 
          {
            int shape_inx = elem_proj.shape_inxs[i];
            proj_value0 += elem_proj.shape_coefs[i] * elem_proj.svals[shape_inx][H2D_HCFE_VALUE0][gip_inx];
            proj_value1 += elem_proj.shape_coefs[i] * elem_proj.svals[shape_inx][H2D_HCFE_VALUE1][gip_inx];
            proj_curl += elem_proj.shape_coefs[i] * elem_proj.svals[shape_inx][H2D_HCFE_CURL][gip_inx];
          }

          ////DEBUG-BEGIN
          //double ref_x = gip_pt[H2D_GIP2D_X] * sub_trf.trf->m[0] + sub_trf.trf->t[0];
          //double ref_y = gip_pt[H2D_GIP2D_Y] * sub_trf.trf->m[1] + sub_trf.trf->t[1];
          //std::complex<double> proj_value0A = 0, proj_value1A = 0, proj_curlA = 0;
          //for(int i = 0; i < elem_proj.num_shapes; i++) 
          {
            //  int shape_inx = elem_proj.shape_inxs[i];
            //  proj_value0A += elem_proj.shape_coefs[i] * shapeset->get_fn_value(shape_inx, ref_x, ref_y, 0);
            //  proj_value1A += elem_proj.shape_coefs[i] * shapeset->get_fn_value(shape_inx, ref_x, ref_y, 1);
            //  proj_curlA += elem_proj.shape_coefs[i] * (shapeset->get_dx_value(shape_inx, ref_x, ref_y, 1) - shapeset->get_dy_value(shape_inx, ref_x, ref_y, 0));
            //}
            //error_if(std::abs(proj_value0 - proj_value0A) > 1E-15
            //  || std::abs(proj_value1 - proj_value1A) > 1E-15
            //  || std::abs(proj_curl - proj_curlA) > 1E-15, "A1");
            ////DEBUG-END

            //get value of ref. solution
            std::complex<double> ref_value0 = sub_trf.coef_mx * sub_gip.rvals[H2D_HCFE_VALUE0][gip_inx];
            std::complex<double> ref_value1 = sub_trf.coef_my * sub_gip.rvals[H2D_HCFE_VALUE1][gip_inx];
            std::complex<double> ref_curl = coef_curl * sub_gip.rvals[H2D_HCFE_CURL][gip_inx]; //coef_curl * curl

            //evaluate error
            double error_squared = sqr(proj_value0 - ref_value0)
              + sqr(proj_value1 - ref_value1)
              + sqr(proj_curl - ref_curl);

            total_error_squared += gip_pt[H2D_GIP2D_W] * error_squared;
          }
        }
        return total_error_squared;
      }
    }
  }
}