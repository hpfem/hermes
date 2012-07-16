#include <iostream>
#include "global.h"
#include "solution.h"
#include "discrete_problem.h"
#include "quad_all.h"
#include "element_to_refine.h"
#include "optimum_selector.h"

#define H2DST_ANY H2DST_VERTEX | H2DST_HORIZ_EDGE | H2DST_VERT_EDGE | H2DST_TRI_EDGE | H2DST_BUBBLE ///< Any type of shape. Used just for masky

namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors
    {
      HERMES_API const char* get_cand_list_str(const CandList cand_list)
      {
        switch(cand_list)
        {
        case H2D_P_ISO: return "P_ISO";
        case H2D_P_ANISO: return "P_ANISO";
        case H2D_H_ISO: return "H_ISO";
        case H2D_H_ANISO: return "H_ANISO";
        case H2D_HP_ISO: return "HP_ISO";
        case H2D_HP_ANISO_H: return "HP_ANISO_H";
        case H2D_HP_ANISO_P: return "HP_ANISO_P";
        case H2D_HP_ANISO: return "HP_ANISO";
        default:
          throw Hermes::Exceptions::Exception("Invalid adapt type %d.", cand_list);
          return NULL;
        }
      }

      HERMES_API bool is_hp(const CandList cand_list)
      {
        switch(cand_list)
        {
        case H2D_P_ISO:
        case H2D_P_ANISO:
        case H2D_H_ISO:
        case H2D_H_ANISO: return false; break;
        case H2D_HP_ISO:
        case H2D_HP_ANISO_H:
        case H2D_HP_ANISO_P:
        case H2D_HP_ANISO: return true; break;
        default: throw Hermes::Exceptions::Exception("Invalid adapt type %d.", cand_list); return false;
        }
      }

      HERMES_API bool is_p_aniso(const CandList cand_list)
      {
        switch(cand_list)
        {
        case H2D_P_ISO: return false;
        case H2D_P_ANISO: return true;
        case H2D_H_ISO: return false;
        case H2D_H_ANISO: return false;
        case H2D_HP_ISO: return false;
        case H2D_HP_ANISO_H: return false;
        case H2D_HP_ANISO_P: return true;
        case H2D_HP_ANISO: return true;
        default: throw Hermes::Exceptions::Exception("Invalid adapt type %d.", cand_list); return false;
        }
      }

      template<typename Scalar>
      OptimumSelector<Scalar>::OptimumSelector(CandList cand_list, double conv_exp, int
        max_order, Shapeset* shapeset, const Range& vertex_order, const
        Range& edge_bubble_order) :
      Selector<Scalar>(max_order),
        opt_symmetric_mesh(true),
        opt_apply_exp_dof(false),
        cand_list(cand_list),
        conv_exp(conv_exp),
        shapeset(shapeset)
      {
        if(shapeset == NULL)
          throw Exceptions::NullException(3);

        num_shapes = new int***[2];
        for(int i = 0; i < 2; i++)
        {
          num_shapes[i] = new int**[H2DRS_MAX_ORDER + 2];
          for(int j = 0; j < H2DRS_MAX_ORDER + 2; j++)
          {
            num_shapes[i][j] = new int*[H2DRS_MAX_ORDER + 2];
            for(int k = 0; k < H2DRS_MAX_ORDER + 2; k++)
            {
              num_shapes[i][j][k] = new int[6];
              memset(num_shapes[i][j][k], 0, 6*sizeof(int));
            }
          }
        }

        //build shape indices
        build_shape_indices(HERMES_MODE_TRIANGLE, vertex_order, edge_bubble_order);
        build_shape_indices(HERMES_MODE_QUAD, vertex_order, edge_bubble_order);
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::add_bubble_shape_index(int order_h, int order_v, std::map<int, bool>& used_shape_index, Hermes::vector<ShapeInx>& indices, ElementMode2D mode)
      {
        int quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v);
        const int num_bubbles = shapeset->get_num_bubbles(quad_order, mode);
        int* bubble_inxs = shapeset->get_bubble_indices(quad_order, mode);
        for(int j = 0; j < num_bubbles; j++)
        {
          int inx_bubble = bubble_inxs[j];
          if(used_shape_index.find(inx_bubble) == used_shape_index.end())
          {
            used_shape_index[inx_bubble] = true;
            indices.push_back(ShapeInx(order_h, order_v, inx_bubble, H2DST_BUBBLE));
            for(int order_h_i = order_h+1; order_h_i < H2DRS_MAX_ORDER + 2; order_h_i++)
              for(int order_v_i = order_v+1; order_v_i < H2DRS_MAX_ORDER + 2; order_v_i++)
              {
                num_shapes[mode][order_h_i][order_v_i][H2DSI_BUBBLE]++;
                num_shapes[mode][order_h_i][order_v_i][H2DSI_ANY]++;
              }
            num_shapes[mode][0][0][H2DSI_BUBBLE]++;
            num_shapes[mode][0][0][H2DSI_ANY]++;
          }
        }
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::build_shape_indices(const ElementMode2D mode, const Range& vertex_order, const Range& edge_bubble_order)
      {
        Hermes::vector<ShapeInx> &indices = shape_indices[mode];
        int* next_order = this->next_order_shape[mode];
        int& max_shape_inx = this->max_shape_inx[mode];
        int num_edges = (mode == HERMES_MODE_QUAD) ? 4 : 3;
        bool &has_vertex = has_vertex_shape[mode];
        bool &has_edge = has_edge_shape[mode];
        bool &has_bubble = has_bubble_shape[mode];

        //cleanup
        indices.clear();
        indices.reserve((H2DRS_MAX_ORDER + 1) * (H2DRS_MAX_ORDER + 1));
        has_vertex = has_edge = has_bubble = false;

        //get total range of orders
        Range order_range = Range::make_envelope(vertex_order, edge_bubble_order);

        //prepare array of used shape indices
        std::map<int, bool> used_shape_index;

        //for all orders
        max_shape_inx = 0;
        int examined_shape = 0;
        for(int i = order_range.lower(); i <= order_range.upper(); i++)
        {
          //vertex functions
          if(vertex_order.is_in_closed(i))
          {
            for (int j = 0; j < num_edges; j++)
            {
              int inx = shapeset->get_vertex_index(j, mode);
              if(inx >= 0)
              {
                used_shape_index[inx] = true;
                indices.push_back(ShapeInx(1, 1, inx, H2DST_VERTEX));
                if(mode == HERMES_MODE_QUAD)
                {
                  for(int order_h_i = 2; order_h_i < H2DRS_MAX_ORDER + 2; order_h_i++)
                    for(int order_v_i = 2; order_v_i < H2DRS_MAX_ORDER + 2; order_v_i++)
                    {
                      num_shapes[mode][order_h_i][order_v_i][H2DSI_VERTEX]++;
                      num_shapes[mode][order_h_i][order_v_i][H2DSI_ANY]++;
                    }
                  num_shapes[mode][0][0][H2DSI_VERTEX]++;
                  num_shapes[mode][0][0][H2DSI_ANY]++;
                }
                else
                {
                  for(int order_h_i = 1; order_h_i < H2DRS_MAX_ORDER + 1; order_h_i++)
                  {
                    num_shapes[mode][order_h_i][0][H2DSI_VERTEX]++;
                    num_shapes[mode][order_h_i][0][H2DSI_ANY]++;
                  }
                  num_shapes[mode][0][0][H2DSI_VERTEX]++;
                  num_shapes[mode][0][0][H2DSI_ANY]++;
                }
                has_vertex = true;
              }
            }
          }

          //edge functions
          if(edge_bubble_order.is_in_closed(i))
          {
            //edge functions
            if(mode == HERMES_MODE_QUAD)
            {
              for (int j = 0; j < num_edges; j++)
              {
                int inx = shapeset->get_edge_index(j, 0, i, HERMES_MODE_QUAD);
                if(inx >= 0)
                {
                  used_shape_index[inx] = true;
                  if((j&1) == 0)//horizontal edge
                  {
                    indices.push_back(ShapeInx(i, 0, inx, H2DST_HORIZ_EDGE));
                    for(int order_h_i = i+1; order_h_i < H2DRS_MAX_ORDER + 2; order_h_i++)
                      for(int order_v_i = 1; order_v_i < H2DRS_MAX_ORDER + 2; order_v_i++)
                      {
                        num_shapes[mode][order_h_i][order_v_i][H2DST_HORIZ_EDGE]++;
                        num_shapes[mode][order_h_i][order_v_i][H2DSI_ANY]++;
                      }
                    num_shapes[mode][0][0][H2DST_HORIZ_EDGE]++;
                    num_shapes[mode][0][0][H2DSI_ANY]++;
                  }
                  else  //vertical edge
                  {
                    indices.push_back(ShapeInx(0, i, inx, H2DST_VERT_EDGE));
                    for(int order_h_i = 1; order_h_i < H2DRS_MAX_ORDER + 2; order_h_i++)
                      for(int order_v_i = i+1; order_v_i < H2DRS_MAX_ORDER + 2; order_v_i++)
                      {
                        num_shapes[mode][order_h_i][order_v_i][H2DSI_VERT_EDGE]++;
                        num_shapes[mode][order_h_i][order_v_i][H2DSI_ANY]++;
                      }
                    num_shapes[mode][0][0][H2DSI_VERT_EDGE]++;
                    num_shapes[mode][0][0][H2DSI_ANY]++;
                  }
                  has_edge = true;
                }
              }
            }
            else
            {
              for (int j = 0; j < num_edges; j++)
              {
                int inx = shapeset->get_edge_index(j, 0, i, HERMES_MODE_TRIANGLE);
                if(inx >= 0)
                {
                  used_shape_index[inx] = true;
                  indices.push_back(ShapeInx(i, i, inx, H2DST_TRI_EDGE));
                  for(int order_h_i = i+1; order_h_i < H2DRS_MAX_ORDER + 2; order_h_i++)
                  {
                    num_shapes[mode][order_h_i][0][H2DSI_TRI_EDGE]++;
                    num_shapes[mode][order_h_i][0][H2DSI_ANY]++;
                  }
                  num_shapes[mode][0][0][H2DSI_TRI_EDGE]++;
                  num_shapes[mode][0][0][H2DSI_ANY]++;
                  has_edge = true;
                }
              }
            }

            //bubble functions
            if(mode == HERMES_MODE_QUAD)
            {
              //NOTE: shapeset returns a set of all possible bubble functions and it is not possible to identify which is the smallest
              // order of an element which contains this function, e.g., in a case of a Hcurl and an element of an order 1/1, it returns
              // a list that contains a function of poly-order 2/0. Also, order of indices is not given.
              int order = i;
              unsigned num_indices_prev = indices.size();
              for(int order_other = edge_bubble_order.lower(); order_other <= order; order_other++)
              {
                add_bubble_shape_index(order, order_other, used_shape_index, indices, mode);
                add_bubble_shape_index(order_other, order, used_shape_index, indices, mode);
              }

              //check if indices were added
              if(num_indices_prev < indices.size())
                has_bubble = true;
            }
            else { //triangles
              int order = i;
              int num_bubbles = shapeset->get_num_bubbles(order, mode);
              int* bubble_inxs = shapeset->get_bubble_indices(order, mode);
              for(int j = 0; j < num_bubbles; j++)
              {
                int inx_bubble = bubble_inxs[j];
                if(used_shape_index.find(inx_bubble) == used_shape_index.end())
                {
                  used_shape_index[inx_bubble] = true;
                  indices.push_back(ShapeInx(order, order, inx_bubble, H2DST_BUBBLE));
                  for(int order_h_i = order+1; order_h_i < H2DRS_MAX_ORDER + 2; order_h_i++)
                  {
                    num_shapes[mode][order_h_i][0][H2DSI_BUBBLE]++;
                    num_shapes[mode][order_h_i][0][H2DSI_ANY]++;
                  }
                  num_shapes[mode][0][0][H2DSI_BUBBLE]++;
                  num_shapes[mode][0][0][H2DSI_ANY]++;
                  has_bubble = true;
                }
              }
            }

            //store index of the next order
            next_order[i] = (int)indices.size();

            //update maximum
            while(examined_shape < next_order[i])
            {
              max_shape_inx = std::max(max_shape_inx, indices[examined_shape].inx);
              examined_shape++;
            }
          }
          else
            next_order[i] = (int)indices.size();
        }
      }

      template<typename Scalar>
      int OptimumSelector<Scalar>::calc_num_shapes(int mode, int order_h, int order_v, int allowed_type_mask)
      {
        //test whether the evaluation is necessary
        bool full_eval = false;
        if((allowed_type_mask & H2DST_VERTEX) != 0)
          full_eval |= has_vertex_shape[mode];
        if((allowed_type_mask & (H2DST_HORIZ_EDGE | H2DST_VERT_EDGE | H2DST_TRI_EDGE)) != 0)
          full_eval |= has_edge_shape[mode];
        if((allowed_type_mask & H2DST_BUBBLE) != 0)
          full_eval |= has_bubble_shape[mode];

        //evaluate
        if(full_eval)
        {
          Hermes::vector<ShapeInx>& shapes = shape_indices[mode];
          int num = 0;
          typename Hermes::vector<ShapeInx>::const_iterator shape = shapes.begin();
          while (shape != shapes.end())
          {
            if(((int)shape->type & allowed_type_mask) != 0)
            {
              if((order_h == H2DRS_ORDER_ANY || shape->order_h <= order_h) && (order_v == H2DRS_ORDER_ANY || shape->order_v <= order_v))
                num++;
            }
            shape++;
          }
          return num;
        }
        else
          return 0;
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::append_candidates_split(const int start_quad_order, const int last_quad_order, const int split, bool iso_p)
      {
        //check whether end orders are not lower than start orders
        if(last_quad_order < 0 || start_quad_order < 0)
          return;
        if(H2D_GET_H_ORDER(start_quad_order) > H2D_GET_H_ORDER(last_quad_order) || H2D_GET_V_ORDER(start_quad_order) > H2D_GET_V_ORDER(last_quad_order))
          return;

        //get number of sons
        const int num_sons = get_refin_sons(split);

        //initialize orders
        int quad_orders[H2D_MAX_ELEMENT_SONS];
        OrderPermutator quad_perms[H2D_MAX_ELEMENT_SONS];
        for(int i = 0; i < num_sons; i++)
        {
          quad_orders[i] = start_quad_order;
          quad_perms[i] = OrderPermutator(start_quad_order, last_quad_order, iso_p, &quad_orders[i]);
        }
        for(int i = num_sons; i < H2D_MAX_ELEMENT_SONS; i++)
          quad_orders[i] = 0;

        //generate permutations of orders
        bool quit = false;
        while(!quit)
        {
          do { //create permutation of son 0
            candidates.push_back(Cand(split, quad_orders));
          } while (quad_perms[0].next());

          //reset son 0
          quad_perms[0].reset();

          //increment orders of other sons
          int inx_son = 1;
          while (inx_son < num_sons && !quad_perms[inx_son].next())
          {
            quad_perms[inx_son].reset(); //reset order of the son
            inx_son++;
          }
          if(inx_son >= num_sons)
            quit = true;
        }
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::create_candidates(Element* e, int quad_order, int max_ha_quad_order, int max_p_quad_order)
      {
        int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);
        int max_p_order_h = H2D_GET_H_ORDER(max_p_quad_order), max_p_order_v = H2D_GET_V_ORDER(max_p_quad_order);
        int max_ha_order_h = H2D_GET_H_ORDER(max_ha_quad_order), max_ha_order_v = H2D_GET_V_ORDER(max_ha_quad_order);
        bool tri = e->is_triangle();

        //clear list of candidates
        candidates.clear();
        if(candidates.capacity() < H2DRS_ASSUMED_MAX_CANDS)
          candidates.reserve(H2DRS_ASSUMED_MAX_CANDS);

        //generate all P-candidates (start from intention of generating all possible candidates
        //and restrict it according to the given adapt-type)
        bool iso_p = false;
        int start_quad_order = quad_order;
        int last_quad_order = H2D_MAKE_QUAD_ORDER(std::min(max_p_order_h, order_h + H2DRS_MAX_ORDER_INC), std::min(max_p_order_v, order_v + H2DRS_MAX_ORDER_INC));
        switch(cand_list)
        {
        case H2D_H_ISO:
        case H2D_H_ANISO: last_quad_order = start_quad_order; break; //no P-candidates except the original candidate
        case H2D_P_ISO:
        case H2D_HP_ISO:
        case H2D_HP_ANISO_H: iso_p = true; break; //iso change of orders
        }
        append_candidates_split(quad_order, last_quad_order, H2D_REFINEMENT_P, tri || iso_p);

        //generate all H-candidates
        iso_p = false;
        int start_order_h = std::max(current_min_order, (order_h + 1) / 2), start_order_v = std::max(current_min_order, (order_v + 1) / 2);
        start_quad_order = H2D_MAKE_QUAD_ORDER(start_order_h, start_order_v);
        last_quad_order = H2D_MAKE_QUAD_ORDER(std::min(max_ha_order_h, std::min(start_order_h + H2DRS_MAX_ORDER_INC, order_h)), std::min(max_ha_order_v, std::min(start_order_v + H2DRS_MAX_ORDER_INC, order_v)));
        switch(cand_list)
        {
        case H2D_H_ISO:
        case H2D_H_ANISO:
          last_quad_order = start_quad_order = quad_order; break; //no only one candidate will be created
        case H2D_P_ISO:
        case H2D_P_ANISO: last_quad_order = -1; break; //no H-candidate will be generated
        case H2D_HP_ISO:
        case H2D_HP_ANISO_H: iso_p = true; break; //iso change of orders
        }
        append_candidates_split(start_quad_order, last_quad_order, H2D_REFINEMENT_H, tri || iso_p);

        //generate all ANISO-candidates
        if(!tri && e->iro_cache < 8 /** \todo Find and why is iro_cache compared with the number 8. What does the number 8 mean? */
          && (cand_list == H2D_H_ANISO || cand_list == H2D_HP_ANISO_H || cand_list == H2D_HP_ANISO))
        {
          iso_p = false;
          int start_quad_order_hz = H2D_MAKE_QUAD_ORDER(order_h, std::max(current_min_order, (order_v + 1) / 2));
          int last_quad_order_hz = H2D_MAKE_QUAD_ORDER(std::min(max_ha_order_h, order_h + H2DRS_MAX_ORDER_INC), std::min(order_v, H2D_GET_V_ORDER(start_quad_order) + H2DRS_MAX_ORDER_INC));
          int start_quad_order_vt = H2D_MAKE_QUAD_ORDER(std::max(current_min_order, (order_h + 1) / 2), order_v);
          int last_quad_order_vt = H2D_MAKE_QUAD_ORDER(std::min(order_h, H2D_GET_H_ORDER(start_quad_order) + H2DRS_MAX_ORDER_INC), std::min(max_ha_order_v, order_v + H2DRS_MAX_ORDER_INC));
          switch(cand_list)
          {
          case H2D_H_ANISO:
            last_quad_order_hz = start_quad_order_hz = quad_order;
            last_quad_order_vt = start_quad_order_vt = quad_order;
            break; //only one candidate will be created
          case H2D_HP_ANISO_H: iso_p = true; break; //iso change of orders
          }
          if(iso_p) { //make orders uniform: take mininmum order since nonuniformity is caused by different handling of orders along directions
            int order = std::min(H2D_GET_H_ORDER(start_quad_order_hz), H2D_GET_V_ORDER(start_quad_order_hz));
            start_quad_order_hz = H2D_MAKE_QUAD_ORDER(order, order);
            order = std::min(H2D_GET_H_ORDER(start_quad_order_vt), H2D_GET_V_ORDER(start_quad_order_vt));
            start_quad_order_vt = H2D_MAKE_QUAD_ORDER(order, order);

            order = std::min(H2D_GET_H_ORDER(last_quad_order_hz), H2D_GET_V_ORDER(last_quad_order_hz));
            last_quad_order_hz = H2D_MAKE_QUAD_ORDER(order, order);
            order = std::min(H2D_GET_H_ORDER(last_quad_order_vt), H2D_GET_V_ORDER(last_quad_order_vt));
            last_quad_order_vt = H2D_MAKE_QUAD_ORDER(order, order);
          }
          append_candidates_split(start_quad_order_hz, last_quad_order_hz, H2D_REFINEMENT_ANISO_H, iso_p);
          append_candidates_split(start_quad_order_vt, last_quad_order_vt, H2D_REFINEMENT_ANISO_V, iso_p);
        }
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::update_cands_info(CandsInfo& info_h, CandsInfo& info_p, CandsInfo& info_aniso) const
      {
        typename Hermes::vector<Cand>::const_iterator cand = candidates.begin();
        while (cand != candidates.end())
        {
          CandsInfo* info = NULL;
          if(cand->split == H2D_REFINEMENT_H) info = &info_h;
          else if(cand->split == H2D_REFINEMENT_P) info = &info_p;
          else if(cand->split == H2D_REFINEMENT_ANISO_H || cand->split == H2D_REFINEMENT_ANISO_V) info = &info_aniso;
          else { throw Hermes::Exceptions::Exception("Invalid candidate type: %d.", cand->split); };

          //evaluate elements of candidates
          const int num_elems = cand->get_num_elems();
          for(int i = 0; i < num_elems; i++)
          {
            int elem_order_h = H2D_GET_H_ORDER(cand->p[i]), elem_order_v = H2D_GET_V_ORDER(cand->p[i]);
            if(elem_order_h != elem_order_v)
              info->uniform_orders = false;
            if(info->min_quad_order < 0 || info->max_quad_order < 0)
              info->min_quad_order = info->max_quad_order = H2D_MAKE_QUAD_ORDER(elem_order_h, elem_order_v);
            else
            {
              info->min_quad_order = H2D_MAKE_QUAD_ORDER(std::min(H2D_GET_H_ORDER(info->min_quad_order), elem_order_h), std::min(H2D_GET_V_ORDER(info->min_quad_order), elem_order_v));
              info->max_quad_order = H2D_MAKE_QUAD_ORDER(std::max(H2D_GET_H_ORDER(info->max_quad_order), elem_order_h), std::max(H2D_GET_V_ORDER(info->max_quad_order), elem_order_v));
            }
          }

          //next candidate
          cand++;
        }
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::evaluate_cands_dof(Element* e, Solution<Scalar>* rsln)
      {
        bool tri = e->is_triangle();

        for (unsigned i = 0; i < candidates.size(); i++)
        {
          Cand& c = candidates[i];
          if(tri) { //triangle
            switch(c.split)
            {
            case H2D_REFINEMENT_H:
              {
                const int central = 3; //central triangle
                c.dofs = 0;
                for(int j = 0; j < H2D_MAX_ELEMENT_SONS; j++)
                {
                  c.dofs += num_shapes[HERMES_MODE_TRIANGLE][H2D_GET_H_ORDER(c.p[j]) + 1][H2DRS_ORDER_ANY + 1][H2DSI_ANY];
                  if(j != central)
                  {
                    c.dofs -= num_shapes[HERMES_MODE_TRIANGLE][std::min(H2D_GET_H_ORDER(c.p[j]), H2D_GET_H_ORDER(c.p[central])) + 1][H2DRS_ORDER_ANY + 1][H2DSI_TRI_EDGE] / 3; //shared edge: since triangle has three edges which are identified by a single order this will find 3 x different edge of a given order
                  }
                }
                if(has_vertex_shape[HERMES_MODE_TRIANGLE])
                  c.dofs -= 2*3; // Every vertex function belonging to vertices of the middle triangle is added 3-times, so it has to be deducted 2 times.
              }
              break;

            case H2D_REFINEMENT_P:
              c.dofs = num_shapes[HERMES_MODE_TRIANGLE][H2D_GET_H_ORDER(c.p[0]) + 1][H2DRS_ORDER_ANY + 1][H2DSI_ANY];
              break;

            default:
              throw Hermes::Exceptions::Exception("Unknown split type \"%d\" at candidate %d (element #%d)", c.split, i, e->id);
            }
          }
          else { //quad
            switch(c.split)
            {
            case H2D_REFINEMENT_H:
              c.dofs = 0;
              for(int j = 0; j < H2D_MAX_ELEMENT_SONS; j++)
                c.dofs += num_shapes[HERMES_MODE_QUAD][H2D_GET_H_ORDER(c.p[j]) + 1][H2D_GET_V_ORDER(c.p[j]) + 1][H2DSI_ANY];
              for(int j = 0; j < 2; j++) { //shared edge functions
                c.dofs -= num_shapes[HERMES_MODE_QUAD][H2DRS_ORDER_ANY + 1][std::min(H2D_GET_V_ORDER(c.p[2*j]), H2D_GET_V_ORDER(c.p[2*j + 1])) + 1][H2DSI_VERT_EDGE] / 2; //shared vertical edge functions: every edge is twice there
                c.dofs -= num_shapes[HERMES_MODE_QUAD][std::min(H2D_GET_H_ORDER(c.p[j]), H2D_GET_H_ORDER(c.p[j^3])) + 1][H2DRS_ORDER_ANY + 1][H2DSI_HORIZ_EDGE] / 2; //shared horizontal edge functions: every edge is twice there
              }
              if(has_vertex_shape[HERMES_MODE_QUAD])
                c.dofs -= 4 + 3; //edge vertex + central vertex
              break;

            case H2D_REFINEMENT_ANISO_H:
              c.dofs = num_shapes[HERMES_MODE_QUAD][H2D_GET_H_ORDER(c.p[0]) + 1][H2D_GET_V_ORDER(c.p[0]) + 1][H2DSI_ANY];
              c.dofs += num_shapes[HERMES_MODE_QUAD][H2D_GET_H_ORDER(c.p[1]) + 1][H2D_GET_V_ORDER(c.p[1]) + 1][H2DSI_ANY];
              c.dofs -= num_shapes[HERMES_MODE_QUAD][std::min(H2D_GET_H_ORDER(c.p[0]), H2D_GET_H_ORDER(c.p[1])) + 1][H2DRS_ORDER_ANY + 1][H2DSI_HORIZ_EDGE] / 2; //shared edge functions
              if(has_vertex_shape[HERMES_MODE_QUAD])
                c.dofs -= 2; //shared vertex functions
              break;

            case H2D_REFINEMENT_ANISO_V:
              c.dofs = num_shapes[HERMES_MODE_QUAD][H2D_GET_H_ORDER(c.p[0]) + 1][H2D_GET_V_ORDER(c.p[0]) + 1][H2DSI_ANY];
              c.dofs += num_shapes[HERMES_MODE_QUAD][H2D_GET_H_ORDER(c.p[1]) + 1][H2D_GET_V_ORDER(c.p[1]) + 1][H2DSI_ANY];
              c.dofs -= num_shapes[HERMES_MODE_QUAD][H2DRS_ORDER_ANY + 1][std::min(H2D_GET_V_ORDER(c.p[0]), H2D_GET_V_ORDER(c.p[1])) + 1][H2DSI_VERT_EDGE] / 2; //shared edge functions
              if(has_vertex_shape[HERMES_MODE_QUAD])
                c.dofs -= 2; //shared vertex functions
              break;

            case H2D_REFINEMENT_P:
              c.dofs = num_shapes[HERMES_MODE_QUAD][H2D_GET_H_ORDER(c.p[0]) + 1][H2D_GET_V_ORDER(c.p[0]) + 1][H2DSI_ANY];
              break;

            default:
              throw Hermes::Exceptions::Exception("Unknown split type \"%d\" at candidate %d", c.split, i);
            }
          }
        }
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::evaluate_candidates(Element* e, Solution<Scalar>* rsln, double* avg_error, double* dev_error)
      {
        evaluate_cands_error(e, rsln, avg_error, dev_error);

        evaluate_cands_dof(e, rsln);

        evaluate_cands_score(e);
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::evaluate_cands_score(Element* e)
      {
        //calculate score of candidates
        Cand& unrefined = candidates[0];
        const int num_cands = (int)candidates.size();
        unrefined.score = 0;
        const double unrefined_dofs_exp = std::pow(unrefined.dofs, conv_exp);
        for (int i = 1; i < num_cands; i++)
        {
          Cand& cand = candidates[i];
          if(cand.error < unrefined.error && cand.dofs > unrefined.dofs)
          {
            double delta_dof_exp = std::pow(cand.dofs - unrefined.dofs, conv_exp);
            if(opt_apply_exp_dof)
              delta_dof_exp = std::pow(cand.dofs, conv_exp) - unrefined_dofs_exp;
            candidates[i].score = (log10(unrefined.error) - log10(cand.error)) / delta_dof_exp;
          }
          else
            candidates[i].score = 0;
        }
      }

      template<typename Scalar>
      bool OptimumSelector<Scalar>::compare_cand_score(const Cand& a, const Cand& b)
      {
        return a.score > b.score;
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::select_best_candidate(Element* e, const double avg_error, const double dev_error, int* selected_cand, int* selected_h_cand)
      {
        // select an above-average candidate with the steepest error decrease

        //sort according to the score
        const int num_cands = (int)candidates.size();
        if(num_cands > 2)
          std::sort(candidates.begin() + 1, candidates.end(), compare_cand_score);

        //select best candidate
        int imax = 1, h_imax = 1;
        if(opt_symmetric_mesh) { //prefer symmetric mesh
          //find first valid score that diffres from the next scores
          while ((imax + 1) < num_cands && std::abs(candidates[imax].score - candidates[imax + 1].score) < H2DRS_SCORE_DIFF_ZERO)
          {
            //find the first candidate with a different score
            Cand& cand_current = candidates[imax];
            int imax_end = imax + 2;
            while (imax_end < num_cands && std::abs(cand_current.score - candidates[imax_end].score) < H2DRS_SCORE_DIFF_ZERO)
              imax_end++;

            imax = imax_end;
          }

          //find valid H-refinement candidate
          h_imax = imax;
          while (h_imax < num_cands && candidates[h_imax].split != H2D_REFINEMENT_H)
            h_imax++;
        }

        //make sure the result is valid: index is valid, selected candidate has a valid score
        if(imax >= num_cands || candidates[imax].score == 0)
          imax = 0;
        if(h_imax >= num_cands || candidates[h_imax].score == 0)
          h_imax = 0;

        //report result
        *selected_cand = imax;
        *selected_h_cand = h_imax;
      }

      template<typename Scalar>
      bool OptimumSelector<Scalar>::select_refinement(Element* element, int quad_order, Solution<Scalar>* rsln, ElementToRefine& refinement)
      {
        //make an uniform order in a case of a triangle
        int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);
        if(element->is_triangle())
        {
          order_v = order_h;
          quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v); //in a case of a triangle, order_v is zero. Set it to order_h in order to simplify the routines.
        }

        //set orders
        set_current_order_range(element);

        // To generate always at least the unchanged candidate.
        if(current_max_order < std::max(order_h, order_v))
          current_max_order = std::max(order_h, order_v);

        //build candidates
        int inx_cand, inx_h_cand;
        create_candidates(element, quad_order
          , H2D_MAKE_QUAD_ORDER(current_max_order, current_max_order)
          , H2D_MAKE_QUAD_ORDER(current_max_order, current_max_order));

        if(candidates.size() > 1) { //there are candidates to choose from
          // evaluate candidates (sum partial projection errors, calculate dofs)
          double avg_error, dev_error;
          evaluate_candidates(element, rsln, &avg_error, &dev_error);

          //select candidate
          select_best_candidate(element, avg_error, dev_error, &inx_cand, &inx_h_cand);
        }
        else { //there is not candidate to choose from, select the original candidate
          inx_cand = 0;
          inx_h_cand = 0;
        }

        //copy result to output
        Cand& cand = candidates[inx_cand];
        Cand& cand_h = candidates[inx_h_cand];
        refinement.split = cand.split;
        ElementToRefine::copy_orders(refinement.p, cand.p);
        if(candidates[inx_h_cand].split == H2D_REFINEMENT_H) { //inx_h_cand points to a candidate which is a H-candidate: copy orders
          ElementToRefine::copy_orders(refinement.q, cand_h.p);
        }
        else { //the index is not H-candidate because not candidate was generate: so, fake orders
          int h_cand_orders[H2D_MAX_ELEMENT_SONS] = { cand_h.p[0], cand_h.p[0], cand_h.p[0], cand_h.p[0] };
          ElementToRefine::copy_orders(refinement.q, h_cand_orders);
        }

        //modify orders in a case of a triangle such that order_v is zero
        if(element->is_triangle())
        {
          for(int i = 0; i < H2D_MAX_ELEMENT_SONS; i++)
          {
            if(!(H2D_GET_V_ORDER(refinement.p[i]) == 0 || H2D_GET_H_ORDER(refinement.p[i]) == H2D_GET_V_ORDER(refinement.p[i])))
              throw Exceptions::Exception("Triangle processed but the resulting order (%d, %d) of son %d is not uniform", H2D_GET_H_ORDER(refinement.p[i]), H2D_GET_V_ORDER(refinement.p[i]), i);
            refinement.p[i] = H2D_MAKE_QUAD_ORDER(H2D_GET_H_ORDER(refinement.p[i]), 0);
            if(!(H2D_GET_V_ORDER(refinement.q[i]) == 0 || H2D_GET_H_ORDER(refinement.q[i]) == H2D_GET_V_ORDER(refinement.q[i])))
              throw Exceptions::Exception("Triangle processed but the resulting q-order (%d, %d) of son %d is not uniform", H2D_GET_H_ORDER(refinement.q[i]), H2D_GET_V_ORDER(refinement.q[i]), i);
            refinement.q[i] = H2D_MAKE_QUAD_ORDER(H2D_GET_H_ORDER(refinement.q[i]), 0);
          }
        }

        if(inx_cand == 0)
          return false;
        else
          return true;
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::generate_shared_mesh_orders(const Element* element, const int orig_quad_order, const int refinement, int tgt_quad_orders[H2D_MAX_ELEMENT_SONS], const int* suggested_quad_orders)
      {
        if(refinement == H2D_REFINEMENT_P)
          throw Exceptions::Exception("P-candidate not supported for updating shared orders");
        const int num_sons = get_refin_sons(refinement);
        if(suggested_quad_orders != NULL)
        {
          for(int i = 0; i < num_sons; i++)
            tgt_quad_orders[i] = suggested_quad_orders[i];
        }
        else
        {
          //calculate new quad_orders
          int quad_order = orig_quad_order;
          if(cand_list != H2D_H_ISO && cand_list != H2D_H_ANISO) { //H_ISO and H_ANISO has to keep given order
            int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);
            switch(refinement)
            {
            case H2D_REFINEMENT_H:
              order_h = std::max(1, (order_h + 1)/2);
              order_v = std::max(1, (order_v + 1)/2);
              break;
            case H2D_REFINEMENT_ANISO_H:
              order_v = std::max(1, 2*(order_v + 1)/3);
              break;
            case H2D_REFINEMENT_ANISO_V:
              order_h = std::max(1, 2*(order_h + 1)/3);
              break;
            }
            if(element->is_triangle())
              quad_order = H2D_MAKE_QUAD_ORDER(order_h, 0);
            else
              quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v);
          }
          for(int i = 0; i < num_sons; i++)
            tgt_quad_orders[i] = quad_order;
        }
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::set_option(const SelOption option, bool enable)
      {
        switch(option)
        {
        case H2D_PREFER_SYMMETRIC_MESH: opt_symmetric_mesh = enable; break;
        case H2D_APPLY_CONV_EXP_DOF: opt_apply_exp_dof = enable; break;
        default: throw Hermes::Exceptions::Exception("Unknown option %d.", (int)option);
        }
      }

      template<typename Scalar>
      OptimumSelector<Scalar>::Range::Range() : empty_range(true) {}

      template<typename Scalar>
      OptimumSelector<Scalar>::Range::Range(const int& lower_bound, const int& upper_bound) : lower_bound(lower_bound), upper_bound(upper_bound), empty_range(lower_bound > upper_bound)
      {
      }

      template<typename Scalar>
      bool OptimumSelector<Scalar>::Range::empty() const
      {
        return empty_range;
      }

      template<typename Scalar>
      const int& OptimumSelector<Scalar>::Range::lower() const
      {
        return lower_bound;
      }

      template<typename Scalar>
      const int& OptimumSelector<Scalar>::Range::upper() const
      {
        return upper_bound;
      }

      template<typename Scalar>
      bool OptimumSelector<Scalar>::Range::is_in_closed(const typename OptimumSelector<Scalar>::Range& range) const
      {
        return (range.lower_bound >= lower_bound && range.upper_bound <= upper_bound);
      }

      template<typename Scalar>
      bool OptimumSelector<Scalar>::Range::is_in_closed(const int& value) const
      {
        return (value >= lower_bound && value <= upper_bound);
      }

      template<typename Scalar>
      bool OptimumSelector<Scalar>::Range::is_in_open(const int& value) const
      {
        return (value > lower_bound && value < upper_bound);
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::Range::enlarge_to_include(const int& value)
      {
        if(empty_range) {
          lower_bound = upper_bound = value;
          empty_range = false;
        }
        else {
          if(lower_bound > value)
            lower_bound = value;
          if(upper_bound < value)
            upper_bound = value;
        }
      }

      template<typename Scalar>
      typename OptimumSelector<Scalar>::Range OptimumSelector<Scalar>::Range::make_envelope(const typename OptimumSelector<Scalar>::Range& a, const typename OptimumSelector<Scalar>::Range& b)
      {
        if(a.empty())
          return b;
        else if(b.empty())
          return a;
        else
          return Range(std::min(a.lower(), b.lower()), std::max(a.upper(), b.upper()));
      }

      template class HERMES_API OptimumSelector<double>;
      template class HERMES_API OptimumSelector<std::complex<double> >;
    }
  }
}