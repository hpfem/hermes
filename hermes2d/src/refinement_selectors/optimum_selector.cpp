#include <iostream>
#include "element_to_refine.h"
#include "optimum_selector.h"
#include "order_permutator.h"

#define H2DST_ANY H2DST_VERTEX | H2DST_HORIZ_EDGE | H2DST_VERT_EDGE | H2DST_TRI_EDGE | H2DST_BUBBLE ///< Any type of shape. Used just for masky

namespace Hermes
{
  namespace Hermes2D
  {
    namespace RefinementSelectors
    {
      template<typename Scalar>
      OptimumSelector<Scalar>::OptimumSelector(CandList cand_list, int
        max_order, Shapeset* shapeset, const Range& vertex_order, const
        Range& edge_bubble_order) : Selector<Scalar>(shapeset->get_min_order(), max_order), cand_list(cand_list), shapeset(shapeset), dof_score_exponent(1.0)
      {
        if(shapeset == NULL)
          throw Exceptions::NullException(3);

        num_shapes = new int***[2];
        for(int i = 0; i < 2; i++)
        {
          /// Why H2DRS_MAX_ORDER + 3?
          /// Because H2DRS_MAX_ORDER ( = 9 ) is the max order refinement selectors support (\todo verify the purpose of this).
          /// But we need to know number of shapes for polynomial degrees 0, ..., 10 ( = 11 ) + 1 special (any order - at 0th position).
          num_shapes[i] = new int**[H2D_NUM_SHAPES_SIZE];
          for(int j = 0; j < H2D_NUM_SHAPES_SIZE; j++)
          {
            num_shapes[i][j] = new int*[H2D_NUM_SHAPES_SIZE];
            for(int k = 0; k < H2D_NUM_SHAPES_SIZE; k++)
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
      OptimumSelector<Scalar>::~OptimumSelector()
      {
        for(int i = 0; i < 2; i++)
        {
          for(int j = 0; j < H2D_NUM_SHAPES_SIZE; j++)
          {
            for(int k = 0; k < H2D_NUM_SHAPES_SIZE; k++)
            {
              delete [] num_shapes[i][j][k];
            }
            delete [] num_shapes[i][j];
          }
          delete [] num_shapes[i];
        }
        delete [] num_shapes;
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
            for(int order_h_i = order_h+1; order_h_i < H2D_NUM_SHAPES_SIZE; order_h_i++)
              for(int order_v_i = order_v+1; order_v_i < H2D_NUM_SHAPES_SIZE; order_v_i++)
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
        indices.reserve((H2DRS_MAX_ORDER + 2) * (H2DRS_MAX_ORDER + 2));
        has_vertex = has_edge = has_bubble = false;

        //get total range of orders
        Range order_range = Range::make_envelope(vertex_order, edge_bubble_order);

        //prepare array of used shape indices
        std::map<int, bool> used_shape_index;

        //for all orders
        max_shape_inx = 0;
        int examined_shape = 0;
        for(int i = order_range.lower(); i <= 10; i++)
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
                  for(int order_h_i = 2; order_h_i < H2D_NUM_SHAPES_SIZE; order_h_i++)
                    for(int order_v_i = 2; order_v_i < H2D_NUM_SHAPES_SIZE; order_v_i++)
                    {
                      num_shapes[mode][order_h_i][order_v_i][H2DSI_VERTEX]++;
                      num_shapes[mode][order_h_i][order_v_i][H2DSI_ANY]++;
                    }
                  num_shapes[mode][0][0][H2DSI_VERTEX]++;
                  num_shapes[mode][0][0][H2DSI_ANY]++;
                }
                else
                {
                  for(int order_h_i = 2; order_h_i < H2D_NUM_SHAPES_SIZE; order_h_i++)
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
                    for(int order_h_i = i+1; order_h_i < H2D_NUM_SHAPES_SIZE; order_h_i++)
                      for(int order_v_i = 1; order_v_i < H2D_NUM_SHAPES_SIZE; order_v_i++)
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
                    for(int order_h_i = 1; order_h_i < H2D_NUM_SHAPES_SIZE; order_h_i++)
                      for(int order_v_i = i+1; order_v_i < H2D_NUM_SHAPES_SIZE; order_v_i++)
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
                  for(int order_h_i = i+1; order_h_i < H2D_NUM_SHAPES_SIZE; order_h_i++)
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
                  for(int order_h_i = order+1; order_h_i < H2D_NUM_SHAPES_SIZE; order_h_i++)
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
      void OptimumSelector<Scalar>::append_candidates_split(Hermes::vector<Cand>& candidates, const int start_order, const int last_order, const int split, bool iso_p)
      {
        //check whether end orders are not lower than start orders
        if(last_order < 0 || start_order < 0)
          return;
        if(H2D_GET_H_ORDER(start_order) > H2D_GET_H_ORDER(last_order) || H2D_GET_V_ORDER(start_order) > H2D_GET_V_ORDER(last_order))
          return;

        //get number of sons
        const int num_sons = get_refin_sons(split);

        //initialize orders
        int quad_orders[H2D_MAX_ELEMENT_SONS];
        OrderPermutator quad_perms[H2D_MAX_ELEMENT_SONS];
        for(int i = 0; i < num_sons; i++)
        {
          quad_orders[i] = start_order;
          quad_perms[i] = OrderPermutator(start_order, last_order, iso_p, &quad_orders[i]);
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
      Hermes::vector<Cand> OptimumSelector<Scalar>::create_candidates(Element* e, int quad_order)
      {
        Hermes::vector<Cand> candidates;

        // Get the current order range.
        int current_min_order, current_max_order;
        this->get_current_order_range(e, current_min_order, current_max_order);

        int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);

        if(current_max_order < std::max(order_h, order_v))
          current_max_order = std::max(order_h, order_v);

        bool tri = e->is_triangle();

        //clear list of candidates
        candidates.reserve(H2DRS_ASSUMED_MAX_CANDS);

        //generate all P-candidates (start from intention of generating all possible candidates
        //and restrict it according to the given adapt-type)
        bool iso_p = false;
        int last_order_h = std::min(current_max_order, order_h + H2DRS_MAX_ORDER_INC), last_order_v = std::min(current_max_order, order_v + H2DRS_MAX_ORDER_INC);
        int last_order = H2D_MAKE_QUAD_ORDER(last_order_h, last_order_v);
        switch(cand_list)
        {
        case H2D_H_ISO:
        case H2D_H_ANISO: last_order = quad_order; break; //no P-candidates except the original candidate
        case H2D_P_ISO:
        case H2D_HP_ISO:
        case H2D_HP_ANISO_H: iso_p = true; break; //iso change of orders
        }
        append_candidates_split(candidates, quad_order, last_order, H2D_REFINEMENT_P, tri || iso_p);

        //generate all H-candidates
        iso_p = false;
        int start_order_h = std::max(current_min_order, order_h - 1), start_order_v = std::max(current_min_order, order_v - 1);
        int start_order = H2D_MAKE_QUAD_ORDER(start_order_h, start_order_v);
        last_order_h = std::min(current_max_order, order_h + H2DRS_MAX_ORDER_INC), last_order_v = std::min(current_max_order, order_v + H2DRS_MAX_ORDER_INC);
        last_order = H2D_MAKE_QUAD_ORDER(last_order_h, last_order_v);
        switch(cand_list)
        {
        case H2D_H_ISO:
        case H2D_H_ANISO:
          last_order = start_order = quad_order; break; //no only one candidate will be created
        case H2D_P_ISO:
        case H2D_P_ANISO: last_order = -1; break; //no H-candidate will be generated
        case H2D_HP_ISO:
        case H2D_HP_ANISO_H: iso_p = true; break; //iso change of orders
        }
        append_candidates_split(candidates, start_order, last_order, H2D_REFINEMENT_H, tri || iso_p);

        //generate all ANISO-candidates
        /** \todo Find and why is iro_cache compared with the number 8. What does the number 8 mean? */
        if(!tri && e->iro_cache < 8 && (cand_list == H2D_H_ANISO || cand_list == H2D_HP_ANISO_H || cand_list == H2D_HP_ANISO))
        {
          iso_p = false;
          int start_order_hz = H2D_MAKE_QUAD_ORDER(order_h, std::max(current_min_order, (order_v - 1)));
          int last_order_hz = H2D_MAKE_QUAD_ORDER(std::min(current_max_order, order_h + H2DRS_MAX_ORDER_INC), std::min(order_v, H2D_GET_V_ORDER(start_order) + H2DRS_MAX_ORDER_INC));
          int start_order_vt = H2D_MAKE_QUAD_ORDER(std::max(current_min_order, (order_h - 1)), order_v);
          int last_order_vt = H2D_MAKE_QUAD_ORDER(std::min(order_h, H2D_GET_H_ORDER(start_order) + H2DRS_MAX_ORDER_INC), std::min(current_max_order, order_v + H2DRS_MAX_ORDER_INC));
          switch(cand_list)
          {
          case H2D_H_ANISO:
            last_order_hz = start_order_hz = quad_order;
            last_order_vt = start_order_vt = quad_order;
            break; //only one candidate will be created
          case H2D_HP_ANISO_H: iso_p = true; break; //iso change of orders
          }
          if(iso_p) { //make orders uniform: take mininmum order since nonuniformity is caused by different handling of orders along directions
            int order = std::min(H2D_GET_H_ORDER(start_order_hz), H2D_GET_V_ORDER(start_order_hz));
            start_order_hz = H2D_MAKE_QUAD_ORDER(order, order);
            order = std::min(H2D_GET_H_ORDER(start_order_vt), H2D_GET_V_ORDER(start_order_vt));
            start_order_vt = H2D_MAKE_QUAD_ORDER(order, order);

            order = std::min(H2D_GET_H_ORDER(last_order_hz), H2D_GET_V_ORDER(last_order_hz));
            last_order_hz = H2D_MAKE_QUAD_ORDER(order, order);
            order = std::min(H2D_GET_H_ORDER(last_order_vt), H2D_GET_V_ORDER(last_order_vt));
            last_order_vt = H2D_MAKE_QUAD_ORDER(order, order);
          }
          append_candidates_split(candidates, start_order_hz, last_order_hz, H2D_REFINEMENT_ANISO_H, iso_p);
          append_candidates_split(candidates, start_order_vt, last_order_vt, H2D_REFINEMENT_ANISO_V, iso_p);
        }

        return candidates;
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::update_cands_info(Hermes::vector<Cand>& candidates, CandsInfo& info_h, CandsInfo& info_p, CandsInfo& info_aniso) const
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
      void OptimumSelector<Scalar>::evaluate_cands_dof(Hermes::vector<Cand>& candidates, Element* e, MeshFunction<Scalar>* rsln)
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
                    c.dofs -= num_shapes[HERMES_MODE_TRIANGLE][std::min(H2D_GET_H_ORDER(c.p[j]) , H2D_GET_H_ORDER(c.p[central])) + 1][H2DRS_ORDER_ANY + 1][H2DSI_TRI_EDGE] / 3; //shared edge: since triangle has three edges which are identified by a single order this will find 3 x different edge of a given order
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
      void OptimumSelector<Scalar>::evaluate_candidates(Hermes::vector<Cand>& candidates, Element* e, MeshFunction<Scalar>* rsln)
      {
        evaluate_cands_error(candidates, e, rsln);

        evaluate_cands_dof(candidates, e, rsln);

        evaluate_cands_score(candidates, e);
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::set_dof_score_exponent(double exponent)
      {
        this->dof_score_exponent = exponent;
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::evaluate_cands_score(Hermes::vector<Cand>& candidates, Element* e)
      {
        // Original candidate.
        Cand& unrefined = candidates[0];
        // Original candidate score is zero.
        unrefined.score = 0;

        for (int i = 1; i < candidates.size(); i++)
        {
          Cand& candidate = candidates[i];

          // We are only interested in candidates decreasing the error.
          if(candidate.error < unrefined.error)
            candidate.score = (log(unrefined.error / candidate.error)) / std::pow(candidate.dofs - unrefined.dofs, this->dof_score_exponent);
          else
            candidate.score = 0;
        }
      }

      template<typename Scalar>
      bool OptimumSelector<Scalar>::compare_cand_score(const Cand& a, const Cand& b)
      {
        return a.score > b.score;
      }

      template<typename Scalar>
      void OptimumSelector<Scalar>::select_best_candidate(Hermes::vector<Cand>& candidates, Element* e, Cand*& best_candidate, Cand* best_candidates_specific_type[4])
      {
        //sort according to the score
        const int num_cands = (int)candidates.size();
        
        if(num_cands > 2)
          std::sort(candidates.begin() + 1, candidates.end(), compare_cand_score);

        // Overall best candidate.
        if(candidates[1].score == 0)
          // This means that the best candidate is the 'do-nothing' one.
          best_candidate = &candidates[0];
        else
          // The one with the highest score is the best.
          best_candidate = &candidates[1];

        for(int i = 0; i < num_cands; i++)
        {
          if(candidates[i].split == H2D_REFINEMENT_P)
          {
            best_candidates_specific_type[H2D_REFINEMENT_P] = &candidates[i];
            break;
          }
        }

        for(int i = 0; i < num_cands; i++)
        {
          if(candidates[i].split == H2D_REFINEMENT_H)
          {
            best_candidates_specific_type[H2D_REFINEMENT_H] = &candidates[i];
            break;
          }
        }

        for(int i = 0; i < num_cands; i++)
        {
          if(candidates[i].split == H2D_REFINEMENT_ANISO_H)
          {
            best_candidates_specific_type[H2D_REFINEMENT_ANISO_H] = &candidates[i];
            break;
          }
        }

        for(int i = 0; i < num_cands; i++)
        {
          if(candidates[i].split == H2D_REFINEMENT_ANISO_V)
          {
            best_candidates_specific_type[H2D_REFINEMENT_ANISO_V] = &candidates[i];
            break;
          }
        }
      }

      template<typename Scalar>
      bool OptimumSelector<Scalar>::select_refinement(Element* element, int quad_order, MeshFunction<Scalar>* rsln, ElementToRefine& refinement)
      {
        //make an uniform order in a case of a triangle
        int order_h = H2D_GET_H_ORDER(quad_order), order_v = H2D_GET_V_ORDER(quad_order);
        if(element->is_triangle())
        {
          order_v = order_h;
          quad_order = H2D_MAKE_QUAD_ORDER(order_h, order_v); //in a case of a triangle, order_v is zero. Set it to order_h in order to simplify the routines.
        }

        //build candidates.
        Hermes::vector<Cand> candidates = create_candidates(element, quad_order);
        //there are candidates to choose from
        Cand* best_candidate;
        Cand* best_candidates_specific_type[4];
        memset(best_candidates_specific_type, 0, 4 * sizeof(Cand*));

        if(candidates.size() > 1)
        { 
          // evaluate candidates (sum partial projection errors, calculate dofs)
          evaluate_candidates(candidates, element, rsln);

          //select candidate
          select_best_candidate(candidates, element, best_candidate, best_candidates_specific_type);
        }

        //there is not candidate to choose from, select the original candidate
        else
        { 
          best_candidate = &candidates[0];
        }

        if(best_candidate == &candidates[0])
          return false;

        //copy result to output
        refinement.split = best_candidate->split;
        ElementToRefine::copy_orders(refinement.refinement_polynomial_order, best_candidate->p);
        for(int i = 1; i < 4; i++)
          if(best_candidates_specific_type[i] != NULL)
            ElementToRefine::copy_orders(refinement.best_refinement_polynomial_order_type[i], best_candidates_specific_type[i]->p);

        ElementToRefine::copy_errors(refinement.errors, best_candidate->errors);

        //modify orders in a case of a triangle such that order_v is zero
        if(element->is_triangle())
        {
          for(int i = 0; i < H2D_MAX_ELEMENT_SONS; i++)
          {
#ifdef _DEBUG
            if(!(H2D_GET_V_ORDER(refinement.refinement_polynomial_order[i]) == 0 || H2D_GET_H_ORDER(refinement.refinement_polynomial_order[i]) == H2D_GET_V_ORDER(refinement.refinement_polynomial_order[i])))
              throw Exceptions::Exception("Triangle processed but the resulting order (%d, %d) of son %d is not uniform", H2D_GET_H_ORDER(refinement.refinement_polynomial_order[i]), H2D_GET_V_ORDER(refinement.refinement_polynomial_order[i]), i);
#endif
            refinement.refinement_polynomial_order[i] = H2D_MAKE_QUAD_ORDER(H2D_GET_H_ORDER(refinement.refinement_polynomial_order[i]), 0);

            for(int poly_order_type_i = 1; poly_order_type_i < 4; poly_order_type_i++)
              refinement.best_refinement_polynomial_order_type[poly_order_type_i][i] = H2D_MAKE_QUAD_ORDER(H2D_GET_H_ORDER(refinement.best_refinement_polynomial_order_type[poly_order_type_i][i]), 0);
          }
        }

        return true;
      }

      template class HERMES_API OptimumSelector<double>;
      template class HERMES_API OptimumSelector<std::complex<double> >;
    }
  }
}
