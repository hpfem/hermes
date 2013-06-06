#include "matrix.h"
#include "l2_proj_based_selector.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors
    {
      template<typename Scalar>
      L2ProjBasedSelector<Scalar>::L2ProjBasedSelector(CandList cand_list, int max_order, L2Shapeset* user_shapeset)
        : ProjBasedSelector<Scalar>(cand_list, max_order, user_shapeset == NULL ? new L2Shapeset() : user_shapeset, Range(1, 1), Range(0, H2DRS_MAX_ORDER)), user_shapeset(user_shapeset == NULL ? false : true)
      {
        if(user_shapeset != NULL)
        {
          this->warn("Warning: The user shapeset provided for the selector has to have a correct copy constructor implemented.");
          this->warn("Warning: The functionality for cloning user shapeset is to be implemented yet.");
        }
      }

      template<typename Scalar>
      L2ProjBasedSelector<Scalar>::~L2ProjBasedSelector()
      {
        if(!this->user_shapeset)
          delete this->shapeset;
      }

      template<typename Scalar>
      void L2ProjBasedSelector<Scalar>::get_current_order_range(Element* element, int& min_order_, int& max_order_)
      {
        int max_element_order = (20 - element->iro_cache)/2 - 1;
        if(this->max_order == H2DRS_DEFAULT_ORDER)
          max_order_ = max_element_order;
        else
          max_order_ = std::min(this->max_order, max_element_order);
        min_order_ = 0;
      }

      template<typename Scalar>
      void L2ProjBasedSelector<Scalar>::precalc_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const Hermes::vector<typename OptimumSelector<Scalar>::ShapeInx>& shapes, const int max_shape_inx, typename ProjBasedSelector<Scalar>::TrfShape& svals, ElementMode2D mode)
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
            shape_exp.allocate(H2D_L2FE_NUM, num_gip_points);

            //for all GIP points
            for(int k = 0; k < num_gip_points; k++)
            {
              //transform coordinates
              double ref_x = gip_points[k][H2D_GIP2D_X] * trf.m[0] + trf.t[0];
              double ref_y = gip_points[k][H2D_GIP2D_Y] * trf.m[1] + trf.t[1];

              //for all expansions: retrieve values
              shape_exp[H2D_L2FE_VALUE][k] = this->shapeset->get_fn_value(inx_shape, ref_x, ref_y, 0, mode);
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
          throw Exceptions::Exception("All transformation processed but identity transformation not found.");
      }

      template<typename Scalar>
      void L2ProjBasedSelector<Scalar>::precalc_ortho_shapes(const double3* gip_points, const int num_gip_points, const Trf* trfs, const int num_noni_trfs, const Hermes::vector<typename OptimumSelector<Scalar>::ShapeInx>& shapes, const int max_shape_inx, typename ProjBasedSelector<Scalar>::TrfShape& svals, ElementMode2D mode)
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
              sum += svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_L2FE_VALUE][k] * svals[H2D_TRF_IDENTITY][inx_shape_j][H2D_L2FE_VALUE][k];
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
                svals[inx_trf][inx_shape_i][H2D_L2FE_VALUE][k] -= product * svals[inx_trf][inx_shape_j][H2D_L2FE_VALUE][k];
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
            //identity transformation has to be the last transformation
            if(!done)
              throw Exceptions::Exception("All transformation processed but identity transformation not found.");
          }

          //normalize
          //calculate norm
          double norm_squared = 0.0;
          for(int k = 0; k < num_gip_points; k++)
          {
            double sum = 0.0;
            sum += sqr(svals[H2D_TRF_IDENTITY][inx_shape_i][H2D_L2FE_VALUE][k]);
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
              svals[inx_trf][inx_shape_i][H2D_L2FE_VALUE][k] /= norm;
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
          //identity transformation has to be the last transformation
          if(!done)
            throw Exceptions::Exception("All transformation processed but identity transformation not found.");
        }
      }

      template<typename Scalar>
      Scalar** L2ProjBasedSelector<Scalar>::precalc_ref_solution(int inx_son, MeshFunction<Scalar>* rsln, Element* element, int intr_gip_order)
      {
        // fill with values
        Scalar** rvals_son = new Scalar*[H2D_L2FE_NUM];
        rvals_son[H2D_L2FE_VALUE] = rsln->get_fn_values(0);
        return rvals_son;
      }

      template<typename Scalar>
      double** L2ProjBasedSelector<Scalar>::build_projection_matrix(double3* gip_points, int num_gip_points,
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

              value += gip_points[j][H2D_GIP2D_W] * (value0*value1);
            }

            matrix_row[k] = value;
          }
        }

        return matrix;
      }

      template<typename Scalar>
      Scalar L2ProjBasedSelector<Scalar>::evaluate_rhs_subdomain(Element* sub_elem, const typename ProjBasedSelector<Scalar>::ElemGIP& sub_gip, const typename ProjBasedSelector<Scalar>::ElemSubTrf& sub_trf, const typename ProjBasedSelector<Scalar>::ElemSubShapeFunc& sub_shape)
      {
        Scalar total_value = 0;
        for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++)
        {
          //get location and transform it
          double3 &gip_pt = sub_gip.gip_points[gip_inx];

          //get value of a shape function
          double shape_value = sub_shape.svals[H2D_L2FE_VALUE][gip_inx];

          //get value of ref. solution
          Scalar ref_value;
          ref_value = sub_gip.rvals[H2D_L2FE_VALUE][gip_inx];

          //evaluate a right-hand value
          Scalar value = (shape_value * ref_value);

          total_value += gip_pt[H2D_GIP2D_W] * value;
        }
        return total_value;
      }

      template<typename Scalar>
      double L2ProjBasedSelector<Scalar>::evaluate_error_squared_subdomain(Element* sub_elem, const typename ProjBasedSelector<Scalar>::ElemGIP& sub_gip, const typename ProjBasedSelector<Scalar>::ElemSubTrf& sub_trf, const typename ProjBasedSelector<Scalar>::ElemProj& elem_proj)
      {
        double total_error_squared = 0;
        for(int gip_inx = 0; gip_inx < sub_gip.num_gip_points; gip_inx++)
        {
          //get location and transform it
          double3 &gip_pt = sub_gip.gip_points[gip_inx];

          //calculate value of projected solution
          Scalar proj_value = 0;
          for(int i = 0; i < elem_proj.num_shapes; i++)
          {
            int shape_inx = elem_proj.shape_inxs[i];
            proj_value += elem_proj.shape_coeffs[i] * elem_proj.svals[shape_inx][H2D_L2FE_VALUE][gip_inx];
          }

          //get value of ref. solution
          Scalar ref_value = sub_gip.rvals[H2D_L2FE_VALUE][gip_inx];

          //evaluate error
          double error_squared = sqr(proj_value - ref_value);

          total_error_squared += gip_pt[H2D_GIP2D_W] * error_squared;
        }

        return total_error_squared;
      }
      template class HERMES_API L2ProjBasedSelector<double>;
      template class HERMES_API L2ProjBasedSelector<std::complex<double> >;
    }
  }
}