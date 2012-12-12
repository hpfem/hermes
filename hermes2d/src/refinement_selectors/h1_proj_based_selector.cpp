#include "global.h"
#include "matrix.h"
#include "solution.h"
#include "shapeset/shapeset_h1_all.h"
#include "element_to_refine.h"
#include "h1_proj_based_selector.h"
namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors
    {
      template<typename Scalar>
      const int H1ProjBasedSelector<Scalar>::H2DRS_MAX_H1_ORDER = H2DRS_MAX_ORDER;

      template<typename Scalar>
      H1ProjBasedSelector<Scalar>::H1ProjBasedSelector(CandList cand_list, double conv_exp, int max_order, H1Shapeset* user_shapeset)
        : ProjBasedSelector<Scalar>(cand_list, conv_exp, max_order, user_shapeset == NULL ? new H1Shapeset() : user_shapeset, typename OptimumSelector<Scalar>::Range(1, 1), typename OptimumSelector<Scalar>::Range(2, H2DRS_MAX_H1_ORDER)), user_shapeset(user_shapeset == NULL ? false : true)
      {
        if(user_shapeset != NULL)
        {
          this->warn("Warn: The user shapeset provided for the selector has to have a correct copy constructor implemented.");
          this->warn("Warning: The functionality for cloning user shapeset is to be implemented yet.");
        }
      }

      template<typename Scalar>
      H1ProjBasedSelector<Scalar>::~H1ProjBasedSelector()
      {
        if(!this->user_shapeset)
          delete this->shapeset;
      }
      
      template<typename Scalar>
      Selector<Scalar>* H1ProjBasedSelector<Scalar>::clone()
      {
        H1ProjBasedSelector<Scalar>* newSelector = new H1ProjBasedSelector(this->cand_list, this->conv_exp, this->max_order, (H1Shapeset*)this->shapeset);
        newSelector->set_error_weights(this->error_weight_h, this->error_weight_p, this->error_weight_aniso);
        newSelector->isAClone = true;
        return newSelector;
      }

      template<typename Scalar>
      void H1ProjBasedSelector<Scalar>::set_current_order_range(Element* element)
      {
        this->current_max_order = this->max_order;
        int max_element_order = (20 - element->iro_cache)/2 - 1;
        if(this->current_max_order == H2DRS_DEFAULT_ORDER)
          this->current_max_order = max_element_order; // default
        else
          this->current_max_order = std::min(this->current_max_order, max_element_order); // user specified
        this->current_min_order = 1;
      }

      template<typename Scalar>
      void H1ProjBasedSelector<Scalar>::precalc_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const Hermes::vector<typename OptimumSelector<Scalar>::ShapeInx>& shapes, const int max_shape_inx, typename ProjBasedSelector<Scalar>::TrfShape& svals, ElementMode2D mode)
      {
        //for all transformations
        bool done = false;
        int inx_trf = 0;
        while (!done && inx_trf < H2D_TRF_NUM)
        {
          //prepare data for processing
          const Trf& trf = trfs[inx_trf];
          Hermes::vector<typename ProjBasedSelector<Scalar>::TrfShapeExp>& trf_svals = svals[inx_trf];

          //allocate
          trf_svals.resize(max_shape_inx + 1);

          //for all shapes
          const int num_shapes = (int)shapes.size();
          for(int i = 0; i < num_shapes; i++)
          {
            int inx_shape = shapes[i].inx;
            typename ProjBasedSelector<Scalar>::TrfShapeExp& shape_exp = trf_svals[inx_shape];

            //allocate
            shape_exp.allocate(H2D_H1FE_NUM, num_gip_points);

            //for all GIP points
            for(int k = 0; k < num_gip_points; k++)
            {
              //transform coordinates
              double ref_x = gip_points[k][H2D_GIP2D_X] * trf.m[0] + trf.t[0];
              double ref_y = gip_points[k][H2D_GIP2D_Y] * trf.m[1] + trf.t[1];

              //for all expansions: retrieve values
              shape_exp[H2D_H1FE_VALUE][k] = this->shapeset->get_fn_value(inx_shape, ref_x, ref_y, 0, mode);
              shape_exp[H2D_H1FE_DX][k] = this->shapeset->get_dx_value(inx_shape, ref_x, ref_y, 0, mode);
              shape_exp[H2D_H1FE_DY][k] = this->shapeset->get_dy_value(inx_shape, ref_x, ref_y, 0, mode);
            }
          }

          //move to the next transformation
          if(inx_trf == H2D_TRF_IDENTITY)
            done = true;
          else
          {
            inx_trf++;
            if(inx_trf >= num_noni_trfs) //if all transformations were processed, move to the identity transformation
              inx_trf = H2D_TRF_IDENTITY;
          }
        }
        if(!done)
          throw Exceptions::Exception("All transformation processed but identity transformation not found."); //identity transformation has to be the last transformation
      }

      template<typename Scalar>
      void H1ProjBasedSelector<Scalar>::precalc_ortho_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const Hermes::vector<typename OptimumSelector<Scalar>::ShapeInx>& shapes, const int max_shape_inx, typename ProjBasedSelector<Scalar>::TrfShape& svals, ElementMode2D mode)
      {
        //calculate values
        precalc_shapes(gip_points, num_gip_points, trfs, num_noni_trfs, shapes, max_shape_inx, svals, mode);

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
              sum += svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_H1FE_VALUE][k] * svals[H2D_TRF_IDENTITY][inx_shape_j][H2D_H1FE_VALUE][k];
              sum += svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_H1FE_DX][k] * svals[H2D_TRF_IDENTITY][inx_shape_j][H2D_H1FE_DX][k];
              sum += svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_H1FE_DY][k] * svals[H2D_TRF_IDENTITY][inx_shape_j][H2D_H1FE_DY][k];
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
                svals[inx_trf][inx_shape_i][H2D_H1FE_VALUE][k] -= product * svals[inx_trf][inx_shape_j][H2D_H1FE_VALUE][k];
                svals[inx_trf][inx_shape_i][H2D_H1FE_DX][k] -= product * svals[inx_trf][inx_shape_j][H2D_H1FE_DX][k];
                svals[inx_trf][inx_shape_i][H2D_H1FE_DY][k] -= product * svals[inx_trf][inx_shape_j][H2D_H1FE_DY][k];
              }

              //move to the next transformation
              if(inx_trf == H2D_TRF_IDENTITY)
                done = true;
              else
              {
                inx_trf++;
                if(inx_trf >= num_noni_trfs) //if all transformations were processed, move to the identity transformation
                  inx_trf = H2D_TRF_IDENTITY;
              }
            }
            if(!done)
              throw Exceptions::Exception("All transformation processed but identity transformation not found."); //identity transformation has to be the last transformation
          }

          //normalize
          //calculate norm
          double norm_squared = 0.0;
          for(int k = 0; k < num_gip_points; k++)
          {
            double sum = 0.0;
            sum += sqr(svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_H1FE_VALUE][k]);
            sum += sqr(svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_H1FE_DX][k]);
            sum += sqr(svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_H1FE_DY][k]);
            norm_squared += gip_points[k][H2D_GIP2D_W] * sum;
          }
          double norm = sqrt(norm_squared);
          if(!finite(1/norm))
            throw Exceptions::Exception("Norm (%g) is almost zero.", norm);

          //for all transformations: normalize
          int inx_trf = 0;
          bool done = false;
          while (!done && inx_trf < H2D_TRF_NUM)
          {
            //for all integration points
            for(int k = 0; k < num_gip_points; k++)
            {
              svals[inx_trf][inx_shape_i][H2D_H1FE_VALUE][k] /= norm;
              svals[inx_trf][inx_shape_i][H2D_H1FE_DX][k] /= norm;
              svals[inx_trf][inx_shape_i][H2D_H1FE_DY][k] /= norm;
            }

            //move to the next transformation
            if(inx_trf == H2D_TRF_IDENTITY)
              done = true;
            else
            {
              inx_trf++;
              if(inx_trf >= num_noni_trfs) //if all transformations were processed, move to the identity transformation
                inx_trf = H2D_TRF_IDENTITY;
            }
          }
          if(!done)
            throw Exceptions::Exception("All transformation processed but identity transformation not found."); //identity transformation has to be the last transformation
        }
      }

      template<typename Scalar>
      Scalar** H1ProjBasedSelector<Scalar>::precalc_ref_solution(int inx_son, Solution<Scalar>* rsln, Element* element, int intr_gip_order)
      {
        //fill with values
        Scalar** rvals_son = precalc_rvals[inx_son];
        rvals_son[H2D_H1FE_VALUE] = rsln->get_fn_values(0);
        rvals_son[H2D_H1FE_DX] = rsln->get_dx_values(0);
        rvals_son[H2D_H1FE_DY] = rsln->get_dy_values(0);

        return rvals_son;
      }

      template<typename Scalar>
      double** H1ProjBasedSelector<Scalar>::build_projection_matrix(double3* gip_points, int num_gip_points,
        const int* shape_inx, const int num_shapes, ElementMode2D mode)
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
              double value0 = this->shapeset->get_value(H2D_FEI_VALUE, shape0_inx, gip_x, gip_y, 0, mode);
              double value1 = this->shapeset->get_value(H2D_FEI_VALUE, shape1_inx, gip_x, gip_y, 0, mode);
              double dx0 = this->shapeset->get_value(H2D_FEI_DX, shape0_inx, gip_x, gip_y, 0, mode);
              double dx1 = this->shapeset->get_value(H2D_FEI_DX, shape1_inx, gip_x, gip_y, 0, mode);
              double dy0 = this->shapeset->get_value(H2D_FEI_DY, shape0_inx, gip_x, gip_y, 0, mode);
              double dy1 = this->shapeset->get_value(H2D_FEI_DY, shape1_inx, gip_x, gip_y, 0, mode);

              value += gip_points[j][H2D_GIP2D_W] * (value0*value1 + dx0*dx1 + dy0*dy1);
            }

            matrix_row[k] = value;
          }
        }

        return matrix;
      }

      template<typename Scalar>
      Scalar H1ProjBasedSelector<Scalar>::evaluate_rhs_subdomain(Element* sub_elem, const typename ProjBasedSelector<Scalar>::ElemGIP& sub_gip, const typename ProjBasedSelector<Scalar>::ElemSubTrf& sub_trf, const typename ProjBasedSelector<Scalar>::ElemSubShapeFunc& sub_shape)
      {
        Scalar total_value = 0;
        for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++)
        {
          double3 &gip_pt = sub_gip.gip_points[gip_inx];

          //get value of a shape function
          double shape_value[H2D_H1FE_NUM] = {0, 0, 0};
          shape_value[H2D_H1FE_VALUE] = sub_shape.svals[H2D_H1FE_VALUE][gip_inx];
          shape_value[H2D_H1FE_DX] = sub_shape.svals[H2D_H1FE_DX][gip_inx];
          shape_value[H2D_H1FE_DY] = sub_shape.svals[H2D_H1FE_DY][gip_inx];

          ////DEBUG-BEGIN
          //double ref_x = gip_pt[H2D_GIP2D_X] * sub_trf.trf->m[0] + sub_trf.trf->t[0];
          //double ref_y = gip_pt[H2D_GIP2D_Y] * sub_trf.trf->m[1] + sub_trf.trf->t[1];
          //Scalar shape_value[H2D_H1FE_NUM] = {0, 0, 0};
          //shape_value[H2D_H1FE_VALUE] = shapeset->get_fn_value(sub_shape.inx, ref_x, ref_y, 0);
          //shape_value[H2D_H1FE_DX] = shapeset->get_dx_value(sub_shape.inx, ref_x, ref_y, 0);
          //shape_value[H2D_H1FE_DY] = shapeset->get_dy_value(sub_shape.inx, ref_x, ref_y, 0);
          //error_if(std::abs(shape_value[H2D_H1FE_VALUE] - shape_valueA[H2D_H1FE_VALUE]) > 1E-15
          //  || std::abs(shape_value[H2D_H1FE_DX] - shape_valueA[H2D_H1FE_DX]) > 1E-15
          //  || std::abs(shape_value[H2D_H1FE_DY] - shape_valueA[H2D_H1FE_DY]) > 1E-15, "A1");
          ////DEBUG-END

          //get value of ref. solution
          Scalar ref_value[H2D_H1FE_NUM];
          ref_value[H2D_H1FE_VALUE] = sub_gip.rvals[H2D_H1FE_VALUE][gip_inx];
          ref_value[H2D_H1FE_DX] = sub_trf.coef_mx * sub_gip.rvals[H2D_H1FE_DX][gip_inx];
          ref_value[H2D_H1FE_DY] = sub_trf.coef_my * sub_gip.rvals[H2D_H1FE_DY][gip_inx];

          //evaluate a right-hand value
          Scalar value = (shape_value[H2D_H1FE_VALUE] * ref_value[H2D_H1FE_VALUE])
            + (shape_value[H2D_H1FE_DX] * ref_value[H2D_H1FE_DX])
            + (shape_value[H2D_H1FE_DY] * ref_value[H2D_H1FE_DY]);

          total_value += gip_pt[H2D_GIP2D_W] * value;
        }
        return total_value;
      }

      template<typename Scalar>
      double H1ProjBasedSelector<Scalar>::evaluate_error_squared_subdomain(Element* sub_elem, const typename ProjBasedSelector<Scalar>::ElemGIP& sub_gip, const typename ProjBasedSelector<Scalar>::ElemSubTrf& sub_trf, const typename ProjBasedSelector<Scalar>::ElemProj& elem_proj)
      {
        double total_error_squared = 0;
        for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++)
        {
          double3 &gip_pt = sub_gip.gip_points[gip_inx];

          //calculate value of projected solution
          Scalar proj_value[H2D_H1FE_NUM] = {0, 0, 0};
          for(int i = 0; i < elem_proj.num_shapes; i++)
          {
            int shape_inx = elem_proj.shape_inxs[i];
            proj_value[H2D_H1FE_VALUE] += elem_proj.shape_coeffs[i] * elem_proj.svals[shape_inx][H2D_H1FE_VALUE][gip_inx];
            proj_value[H2D_H1FE_DX] += elem_proj.shape_coeffs[i] * elem_proj.svals[shape_inx][H2D_H1FE_DX][gip_inx];
            proj_value[H2D_H1FE_DY] += elem_proj.shape_coeffs[i] * elem_proj.svals[shape_inx][H2D_H1FE_DY][gip_inx];
          }

          ////DEBUG-BEGIN
          //double ref_x = gip_pt[H2D_GIP2D_X] * sub_trf.trf->m[0] + sub_trf.trf->t[0];
          //double ref_y = gip_pt[H2D_GIP2D_Y] * sub_trf.trf->m[1] + sub_trf.trf->t[1];
          //Scalar proj_valueA[H2D_H1FE_NUM] = {0, 0, 0};
          //for(int i = 0; i < elem_proj.num_shapes; i++)
          {
            //  int shape_inx = elem_proj.shape_inxs[i];
            //  proj_valueA[H2D_H1FE_VALUE] += elem_proj.shape_coeffs[i] * shapeset->get_fn_value(shape_inx, ref_x, ref_y, 0);
            //  proj_valueA[H2D_H1FE_DX] += elem_proj.shape_coeffs[i] * shapeset->get_dx_value(shape_inx, ref_x, ref_y, 0);
            //  proj_valueA[H2D_H1FE_DY] += elem_proj.shape_coeffs[i] * shapeset->get_dy_value(shape_inx, ref_x, ref_y, 0);
            //}
            //error_if(std::abs(proj_value[H2D_H1FE_VALUE] - proj_valueA[H2D_H1FE_VALUE]) > 1E-15
            //  || std::abs(proj_value[H2D_H1FE_DX] - proj_valueA[H2D_H1FE_DX]) > 1E-15
            //  || std::abs(proj_value[H2D_H1FE_DY] - proj_valueA[H2D_H1FE_DY]) > 1E-15, "A2");
            ////DEBUG-END

            //get value of ref. solution
            Scalar ref_value[3];
            ref_value[H2D_H1FE_VALUE] = sub_gip.rvals[H2D_H1FE_VALUE][gip_inx];
            ref_value[H2D_H1FE_DX] = sub_trf.coef_mx * sub_gip.rvals[H2D_H1FE_DX][gip_inx];
            ref_value[H2D_H1FE_DY] = sub_trf.coef_my * sub_gip.rvals[H2D_H1FE_DY][gip_inx];

            //evaluate error
            double error_squared = sqr(proj_value[H2D_H1FE_VALUE] - ref_value[H2D_H1FE_VALUE])
              + sqr(proj_value[H2D_H1FE_DX] - ref_value[H2D_H1FE_DX])
              + sqr(proj_value[H2D_H1FE_DY] - ref_value[H2D_H1FE_DY]);

            total_error_squared += gip_pt[H2D_GIP2D_W] * error_squared;
          }
        }
        return total_error_squared;
      }
      template class HERMES_API H1ProjBasedSelector<double>;
      template class HERMES_API H1ProjBasedSelector<std::complex<double> >;
    }
  }
}