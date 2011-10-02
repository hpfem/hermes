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

#include "hermes2d_common_defs.h"
#include "mesh.h"
#include "transformable.h"
#include "traverse.h"
namespace Hermes
{
  namespace Hermes2D
  {
    const uint64_t ONE = (uint64_t) 1 << 63;

    struct Rect
    {
      uint64_t l, b, r, t;
    };

    struct State
    {
      bool visited;
      Element** e;
      Rect  cr;
      Rect* er;
      bool bnd[3];
      uint64_t lo[3], hi[3];
      int* trans;
    };


    int get_split_and_sons(Element* e, Rect* cr, Rect* er, int4& sons)
    {
      uint64_t hmid = (er->l + er->r) >> 1;
      uint64_t vmid = (er->t + er->b) >> 1;

      if (e->bsplit())
      {
        if (cr->r <= hmid && cr->t <= vmid)
          return (sons[0] = sons[1] = sons[2] = sons[3] = 0), 0;
        else if (cr->l >= hmid && cr->t <= vmid)
          return (sons[0] = sons[1] = sons[2] = sons[3] = 1), 0;
        else if (cr->l >= hmid && cr->b >= vmid)
          return (sons[0] = sons[1] = sons[2] = sons[3] = 2), 0;
        else if (cr->r <= hmid && cr->b >= vmid)
          return (sons[0] = sons[1] = sons[2] = sons[3] = 3), 0;
        else if (cr->r <= hmid)
          return (sons[0] = sons[1] = 0, sons[2] = sons[3] = 3), 1;
        else if (cr->l >= hmid)
          return (sons[0] = sons[1] = 1, sons[2] = sons[3] = 2), 1;
        else if (cr->t <= vmid)
          return (sons[0] = sons[3] = 0, sons[1] = sons[2] = 1), 2;
        else if (cr->b >= vmid)
          return (sons[0] = sons[3] = 3, sons[1] = sons[2] = 2), 2;
        else
          return (sons[0] = 0, sons[1] = 1, sons[2] = 2, sons[3] = 3), 3;
      }
      else if (e->hsplit())
      {
        if (cr->t <= vmid)
          return (sons[0] = sons[1] = sons[2] = sons[3] = 4), 0;
        else if (cr->b >= vmid)
          return (sons[0] = sons[1] = sons[2] = sons[3] = 5), 0;
        else
          return (sons[0] = sons[1] = 4, sons[2] = sons[3] = 5), 1;
      }
      else // e->vsplit()
      {
        if (cr->r <= hmid)
          return (sons[0] = sons[1] = sons[2] = sons[3] = 6), 0;
        else if (cr->l >= hmid)
          return (sons[0] = sons[1] = sons[2] = sons[3] = 7), 0;
        else
          return (sons[0] = sons[3] = 6, sons[1] = sons[2] = 7), 2;
      }
    }


    static void move_to_son(Rect* rnew, Rect* rold, int son)
    {
      uint64_t hmid = (rold->l + rold->r) >> 1;
      uint64_t vmid = (rold->t + rold->b) >> 1;
      if(rnew != rold)
        memcpy(rnew, rold, sizeof(Rect));

      switch (son)
      {
      case 0: rnew->r = hmid; rnew->t = vmid; break;
      case 1: rnew->l = hmid; rnew->t = vmid; break;
      case 2: rnew->l = hmid; rnew->b = vmid; break;
      case 3: rnew->r = hmid; rnew->b = vmid; break;
      case 4: rnew->t = vmid; break;
      case 5: rnew->b = vmid; break;
      case 6: rnew->r = hmid; break;
      case 7: rnew->l = hmid; break;
      }
    }


    void init_transforms(Transformable* fn, Rect* cr, Rect* er)
    {
      Rect r;
      memcpy(&r, er, sizeof(Rect));

      while (cr->l > r.l || cr->r < r.r || cr->b > r.b || cr->t < r.t)
      {
        uint64_t hmid = (r.l + r.r) >> 1;
        uint64_t vmid = (r.t + r.b) >> 1;
        int son;

        if (cr->r <= hmid && cr->t <= vmid) son = 0;
        else if (cr->l >= hmid && cr->t <= vmid) son = 1;
        else if (cr->l >= hmid && cr->b >= vmid) son = 2;
        else if (cr->r <= hmid && cr->b >= vmid) son = 3;
        else if (cr->r <= hmid) son = 6;
        else if (cr->l >= hmid) son = 7;
        else if (cr->t <= vmid) son = 4;
        else if (cr->b >= vmid) son = 5;
        else assert(0);

        fn->push_transform(son);
        move_to_son(&r, &r, son);
      }
    }


    State* Traverse::push_state()
    {
      if (top >= size) error("Stack overflow. Increase stack size.");

      if (stack[top].e == NULL)
      {
        stack[top].e = new Element*[num];
        stack[top].er = new Rect[num];
        stack[top].trans = new int[num];
      }

      stack[top].visited = false;
      memset(stack[top].trans, 0, num * sizeof(int));
      return stack + top++;
    }


    void Traverse::set_boundary_info(State* s, bool* bnd, SurfPos* surf_pos)
    {
      Element* e = NULL;
      for (int i = 0; i < num; i++)
        if ((e = s->e[i]) != NULL) break;

      if (tri)
      {
        for (int i = 0; i < 3; i++)
        {
          if ((bnd[i] = (s->bnd[i] && e->en[i]->bnd)))
          {
            surf_pos[i].lo = (double) s->lo[i] / ONE;
            surf_pos[i].hi = (double) s->hi[i] / ONE;
          }
        }
      }
      else
      {
        bnd[0] = (s->cr.b == 0)   && e->en[0]->bnd;
        bnd[1] = (s->cr.r == ONE) && e->en[1]->bnd;
        bnd[2] = (s->cr.t == ONE) && e->en[2]->bnd;
        bnd[3] = (s->cr.l == 0)   && e->en[3]->bnd;

        if (bnd[0]) { surf_pos[0].lo = (double) s->cr.l / ONE;        surf_pos[0].hi = (double) s->cr.r / ONE; }
        if (bnd[1]) { surf_pos[1].lo = (double) s->cr.b / ONE;        surf_pos[1].hi = (double) s->cr.t / ONE; }
        if (bnd[2]) { surf_pos[2].lo = (double) (ONE-s->cr.r) / ONE;  surf_pos[2].hi = (double) (ONE-s->cr.l) / ONE; }
        if (bnd[3]) { surf_pos[3].lo = (double) (ONE-s->cr.t) / ONE;  surf_pos[3].hi = (double) (ONE-s->cr.b) / ONE; }
      }

      for (unsigned int i = 0; i < base->get_num_surf(); i++)
      {
        int j = base->next_vert(i);
        surf_pos[i].v1 = base->vn[i]->id;
        surf_pos[i].v2 = base->vn[j]->id;
        surf_pos[i].marker = e->en[i]->marker;
        surf_pos[i].surf_num = i;
      }
    }


    Element** Traverse::get_next_state(bool* bnd, SurfPos* surf_pos)
    {
      while (1)
      {
        int i, j, son;
        // When the stack of states is not empty (i.e. not at the beginning) the function starts here.
        // If the top state was visited already, we are returning through it:
        // undo all its transformations, pop it and continue with a non-visited one
        State* s;
        while (top > 0 && (s = stack + top-1)->visited)
        {
          if (fn != NULL)
          {
            // For every mesh.
            for (i = 0; i < num; i++)
              if (s->e[i] != NULL)
              {
                // If the element on i-th mesh is transformed to a subelement.
                if (s->trans[i] > 0)
                {
                  if (fn[i]->get_transform() == subs[i])
                    fn[i]->pop_transform();
                  subs[i] = fn[i]->get_transform();
                }
                // If the element on the i-th mesh is active, we have to reset transformation for all functions on the i-th mesh.
                else if (s->trans[i] < 0)
                  fn[i]->reset_transform();
              }
          }
          // Since we are about to process this state and return its elements, we want to take it out of the stack,
          // i.e. to lower the stack's top.
          top--;
        }

        // The stack is empty, take next base element
        // The process starts here (at the beginning the stack is always empty, i.e. top == 0)
        if (top <= 0)
        {
          // Push the state of a new base element.
          // This function only allocates memory for the new state,
          // with as many Elements* as there are meshes in this stage.
          // (Traverse knows what stage it is, because it is being initialized by calling trav.begin(..)).
          s = push_state();
          static const Rect H2D_UNITY = { 0, 0, ONE, ONE };
          s->cr = H2D_UNITY;
          while (1)
          {
            // No more base elements? we're finished.
            // Id is set to zero at the beginning by the function trav.begin(..).
            if (id >= meshes[0]->get_num_base_elements())
              return NULL;
            int nused = 0;
            // The variable num is the number of meshes in the stage
            for (i = 0; i < num; i++)
            {
              // Retrieve the Element with this id on the i-th mesh.
              s->e[i] = meshes[i]->get_element(id);
              if (!s->e[i]->used)
              {
                s->e[i] = NULL;
                continue;
              }
              if (s->e[i]->active && fn != NULL)
                // Important, sets the active element for all functions that share the i-th mesh
                // (PrecalcShapesets, Solutions from previous time/Newton iterations)
                fn[i]->set_active_element(s->e[i]);
              s->er[i] = H2D_UNITY;
              subs[i] = 0;
              nused++;
              base = s->e[i];
            }
            // If there is any used element in this stage we continue with the calculation
            // (break this cycle looking for such an element id).
            if (nused)
              break;
            id++;
          }

          // Sets necessary things for when the new base element is a triangle.
          tri = base->is_triangle();
          id++;

          if (tri)
          {
            for (i = 0; i < 3; i++)
            {
              s->bnd[i] = true;
              s->lo[i] = 0;
              s->hi[i] = ONE;
            }
          }
        }

        // Entering a new state: perform the transformations for it
        s->visited = true;
        if (fn != NULL)
        {
          // For every mesh of this stage..
          for (i = 0; i < num; i++)
            // ..where the element is used ..
            if (s->e[i] != NULL)
            {
              // ..if the element on i-th mesh is inactive (we have to go down)..
              if (s->trans[i] > 0)
              {
                if (fn[i]->get_transform() == subs[i])
                  fn[i]->push_transform(s->trans[i]-1);
                subs[i] = fn[i]->get_transform();
              }
              // ..and when it is active we have to activate on it all functions on i-th mesh
              // (PrecalcShapesets, Solutions from previous time/Newton iterations)
              else if (s->trans[i] < 0)
              {
                fn[i]->set_active_element(s->e[i]);
                if (!tri)
                  init_transforms(fn[i], &s->cr, s->er + i);
                subs[i] = fn[i]->get_transform();
              }
            }
        }

        // Is this the leaf state?
        bool leaf = true;
        for (i = 0; i < num; i++)
          if (s->e[i] != NULL)
            if (!s->e[i]->active)
            {
              leaf = false;
              break;
            }

            // if yes, set boundary flags and return the state
            if (leaf)
            {
              if (bnd != NULL)
                set_boundary_info(s, bnd, surf_pos);
              return s->e;
            }

            // Triangle: push son states
            if (tri)
            {
              // Triangle always has 4 sons.
              for (son = 0; son <= 3; son++)
              {
                State* ns = push_state();
                // For every mesh..
                for (i = 0; i < num; i++)
                {
                  // ..if the element is not used.
                  if (s->e[i] == NULL)
                  {
                    ns->e[i] = NULL;
                  }
                  // ..if the current element is active.
                  else if (s->e[i]->active)
                  {
                    ns->e[i] = s->e[i];
                    ns->trans[i] = son + 1;
                  }
                  // ..we move to the son.
                  else
                  {
                    ns->e[i] = s->e[i]->sons[son];
                    // If the son's element is active.
                    if (ns->e[i]->active)
                      ns->trans[i] = -1;
                  }
                }

                // Determine boundary flags and positions for the new state.
                if (son < 3)
                {
                  memcpy(ns->bnd, s->bnd, sizeof(ns->bnd));
                  memcpy(ns->lo, s->lo, sizeof(ns->lo));
                  memcpy(ns->hi, s->hi, sizeof(ns->hi));

#define mid(n) ((s->lo[n] + s->hi[n]) / 2)
                  switch (son)
                  {
                  case 0: ns->bnd[1] = false; ns->hi[0] = mid(0); ns->lo[2] = mid(2); break;
                  case 1: ns->bnd[2] = false; ns->lo[0] = mid(0); ns->hi[1] = mid(1); break;
                  case 2: ns->bnd[0] = false; ns->lo[1] = mid(1); ns->hi[2] = mid(2); break;
                  }
                }
                else
                {
                  memset(ns->bnd, 0, sizeof(ns->bnd));
                }
              }
            }
            // Quad: this is a little more complicated, same principle, though.
            else
            {
              // Obtain split types and son numbers for the current rectangle on all elements.
              int split = 0;
              for (i = 0; i < num; i++)
                if (s->e[i] != NULL && !s->e[i]->active)
                  split |= get_split_and_sons(s->e[i], &s->cr, s->er + i, sons[i]);

              // Both splits: recur to four sons, similar to triangles.
              if (split == 3)
              {
                for (son = 0; son <= 3; son++)
                {
                  State* ns = push_state();
                  // Sets the son's "base" rectangle to the correct one.
                  move_to_son(&ns->cr, &s->cr, son);

                  for (i = 0; i < num; i++)
                  {
                    if (s->e[i] == NULL)
                    {
                      ns->e[i] = NULL;
                    }
                    else if (s->e[i]->active)
                    {
                      ns->e[i] = s->e[i];
                      ns->trans[i] = son + 1;
                    }
                    else
                    {
                      ns->e[i] = s->e[i]->sons[sons[i][son] & 3];
                      // Sets the son's "current mesh" rectangle correctly.
                      move_to_son(ns->er + i, s->er + i, sons[i][son]);
                      if (ns->e[i]->active)
                        ns->trans[i] = -1;
                    }
                  }
                }
              }
              // V or h split, recur to two sons.
              else if (split > 0)
              {
                int son0 = 4, son1 = 5;
                if (split == 2) { son0 = 6; son1 = 7; }

                for (son = son0; son <= son1; son++)
                {
                  State* ns = push_state();
                  move_to_son(&ns->cr, &s->cr, son);

                  j = (son == 4 || son == 6) ? 0 : 2;
                  for (i = 0; i < num; i++)
                  {
                    if (s->e[i] == NULL)
                    {
                      ns->e[i] = NULL;
                    }
                    else if (s->e[i]->active)
                    {
                      ns->e[i] = s->e[i];
                      ns->trans[i] = son + 1;
                    }
                    else
                    {
                      ns->e[i] = s->e[i]->sons[sons[i][j] & 3];
                      move_to_son(ns->er + i, s->er + i, sons[i][j]);
                      if (ns->e[i]->active) ns->trans[i] = -1;
                    }
                  }
                }
              }
              // No splits, recur to one son.
              else
              {
                State* ns = push_state();
                memcpy(&ns->cr, &s->cr, sizeof(Rect));

                for (i = 0; i < num; i++)
                {
                  if (s->e[i] == NULL)
                  {
                    ns->e[i] = NULL;
                  }
                  else if (s->e[i]->active)
                  {
                    ns->e[i] = s->e[i];
                  }
                  else
                  {
                    ns->e[i] = s->e[i]->sons[sons[i][0] & 3];
                    move_to_son(ns->er + i, s->er + i, sons[i][0]);
                    if (ns->e[i]->active) ns->trans[i] = -1;
                  }
                }
              }
            }
      }
    }


    void Traverse::begin(int n, Mesh** meshes, Transformable** fn)
    {
      //if (stack != NULL) finish();

      assert(n > 0);
      num = n;

      this->meshes = meshes;
      this->fn = fn;

      top = 0;
      size = 256;
      stack = new State[size];
      memset(stack, 0, size * sizeof(State));

      sons = new int4[num];
      subs = new uint64_t[num];
      id = 0;

#ifndef H2D_DISABLE_MULTIMESH_TESTS
      // Test whether all master meshes have the same number of elements.
      int base_elem_num = meshes[0]->get_num_base_elements();
      for (int i = 1; i < n; i++)
        if (base_elem_num != meshes[i]->get_num_base_elements())
          error("Meshes not compatible in Traverse::begin().");

      // Test whether areas of corresponding elements are the same.
      double *areas = new double [base_elem_num];
      memset(areas, 0, base_elem_num*sizeof(double));
      if (areas == NULL) error("Not enough memory in Traverse::begin().");
      // Read base element areas from the first mesh,
      // Also get minimum element area.
      int counter = 0;
      double min_elem_area = 1e30;
      Element* e;
      for_all_base_elements(e, meshes[0])
      {
        areas[counter] = e->get_area();
        if (areas[counter] < min_elem_area) min_elem_area = areas[counter];
        //printf("base_element[%d].area = %g\n", counter, areas[counter]);
        counter++;
      }
      // take one mesh at a time and compare element areas to the areas[] array
      double tolerance = min_elem_area/100.;
      for (int i = 1; i < n; i++)
      {
        counter = 0;
        for_all_base_elements(e, meshes[i])
        {
          if (fabs(areas[counter] - e->get_area()) > tolerance && areas[counter] != 0)
          {
            printf("counter = %d, area_1 = %g, area_2 = %g.\n", counter, areas[counter], e->get_area());
            error("Meshes not compatible in Traverse::begin().");
          }
          counter++;
        }
      }
      delete [] areas;
#endif
    }


    void free_state(State* state)
    {
      delete [] state->e;
      delete [] state->er;
      delete [] state->trans;
      memset(state, 0, sizeof(State));
    }


    void Traverse::finish()
    {
      if (stack == NULL) return;

      for (int i = 0; i < size; i++)
        if (stack[i].e != NULL)
          free_state(stack + i);

      delete [] stack;
      stack = NULL;

      delete [] subs;
      delete [] sons;
    }



    //// union mesh ////////////////////////////////////////////////////////////////////////////////////

    uint64_t Traverse::init_idx(Rect* cr, Rect* er)
    {
      Rect r;
      memcpy(&r, er, sizeof(Rect));

      uint64_t idx = 0;
      while (cr->l > r.l || cr->r < r.r || cr->b > r.b || cr->t < r.t)
      {
        uint64_t hmid = (r.l + r.r) >> 1;
        uint64_t vmid = (r.t + r.b) >> 1;
        int son = -1;

        if (cr->r <= hmid && cr->t <= vmid) { son = 0; r.r = hmid; r.t = vmid; }
        else if (cr->l >= hmid && cr->t <= vmid) { son = 1; r.l = hmid; r.t = vmid; }
        else if (cr->l >= hmid && cr->b >= vmid) { son = 2; r.l = hmid; r.b = vmid; }
        else if (cr->r <= hmid && cr->b >= vmid) { son = 3; r.r = hmid; r.b = vmid; }
        else if (cr->t <= vmid) { son = 4; r.t = vmid; }
        else if (cr->b >= vmid) { son = 5; r.b = vmid; }
        else if (cr->r <= hmid) { son = 6; r.r = hmid; }
        else if (cr->l >= hmid) { son = 7; r.l = hmid; }
        else assert(0);

        idx = (idx << 3) + son + 1;
      }
      return idx;
    }


    void Traverse::union_recurrent(Rect* cr, Element** e, Rect* er, uint64_t* idx, Element* uni)
    {
      int i, j, son;

      // are we at the bottom?
      bool leaf = true;
      for (i = 0; i < num; i++)
        if (!e[i]->active)
        {
          leaf = false;
          break;
        }

      // if yes, store the element transformation indices
      if (leaf)
      {
        if (udsize <= uni->id)
        {
          if (!udsize) udsize = 1024;
          while (udsize <= uni->id) udsize *= 2;
          for (i = 0; i < num; i++)
            unidata[i] = (UniData*) realloc(unidata[i], udsize * sizeof(UniData));
        }
        for (i = 0; i < num; i++)
        {
          unidata[i][uni->id].e = e[i];
          unidata[i][uni->id].idx = idx[i];
        }
        return;
      }

      // state arrays
      Element** e_new = new Element*[num];
      Rect* er_new = new Rect[num];
      Rect cr_new;

      int4* sons = new int4[num];
      uint64_t* idx_new = new uint64_t[num];
      memcpy(idx_new, idx, num*sizeof(uint64_t));

      if (tri)
      {
        // visit all sons of the triangle
        unimesh->refine_element_id(uni->id);
        for (son = 0; son <= 3; son++)
        {
          for (i = 0; i < num; i++)
          {
            if (e[i]->active)
            {
              e_new[i] = e[i];
              idx_new[i] = (idx[i] << 3) + son + 1;
            } else
              e_new[i] = e[i]->sons[son];
          }
          union_recurrent(NULL, e_new, NULL, idx_new, uni->sons[son]);
        }
      }
      else
      {
        // obtain split types and son numbers for the current rectangle on all elements
        int split = 0;
        for (i = 0; i < num; i++)
          if (!e[i]->active)
            split |= get_split_and_sons(e[i], cr, er + i, sons[i]);

        // both splits: recur to four sons
        if (split == 3)
        {
          unimesh->refine_element_id(uni->id, 0);

          for (son = 0; son <= 3; son++)
          {
            move_to_son(&cr_new, cr, son);
            for (i = 0; i < num; i++)
            {
              if (e[i]->active)
              {
                e_new[i] = e[i];
                idx_new[i] = (idx[i] << 3) + son + 1;
              } else
              {
                e_new[i] = e[i]->sons[sons[i][son] & 3];
                move_to_son(&(er_new[i]), er + i, sons[i][son]);
                if (e_new[i]->active) idx_new[i] = init_idx(&cr_new, &(er_new[i]));
              }
            }
            union_recurrent(&cr_new, e_new, er_new, idx_new, uni->sons[son]);
          }
        }
        // v or h split, recur to two sons
        else if (split > 0)
        {
          unimesh->refine_element_id(uni->id, split);

          int son0 = 4, son1 = 5;
          if (split == 2) { son0 = 6; son1 = 7; }

          for (son = son0; son <= son1; son++)
          {
            move_to_son(&cr_new, cr, son);
            j = (son == 4 || son == 6) ? 0 : 2;
            for (i = 0; i < num; i++)
            {
              if (e[i]->active)
              {
                e_new[i] = e[i];
                idx_new[i] = (idx[i] << 3) + son + 1;
              } else
              {
                e_new[i] = e[i]->sons[sons[i][j] & 3];
                move_to_son(&(er_new[i]), er + i, sons[i][j]);
                if (e_new[i]->active) idx_new[i] = init_idx(&cr_new, &(er_new[i]));
              }
            }
            union_recurrent(&cr_new, e_new, er_new, idx_new, uni->sons[son & 3]);
          }
        }
        // no splits, recur to one son
        else
        {
          memcpy(&cr_new, cr, sizeof(Rect));
          for (i = 0; i < num; i++)
          {
            if (e[i]->active)
              e_new[i] = e[i];
            else
            {
              e_new[i] = e[i]->sons[sons[i][0] & 3];
              move_to_son(&(er_new[i]), er + i, sons[i][0]);
              if (e_new[i]->active) idx_new[i] = init_idx(&cr_new, &(er_new[i]));
            }
          }
          union_recurrent(&cr_new, e_new, er_new, idx_new, uni);
        }
      }

      delete [] e_new;
      delete [] er_new;
      delete [] sons;
      delete [] idx_new;
    }


    UniData** Traverse::construct_union_mesh(Mesh* unimesh)
    {
      int i;
      Element** e = new Element*[num];
      Rect* er = new Rect[num];
      Rect cr;

      this->unimesh = unimesh;
      unimesh->copy_base(meshes[0]);

      udsize = 0;
      unidata = new UniData*[num];
      memset(unidata, 0, sizeof(UniData*) * num);

      uint64_t* idx = new uint64_t[num];
      memset(idx, 0, num*sizeof(uint64_t));

      for (id = 0; id < meshes[0]->get_num_base_elements(); id++)
      {
        for (i = 0; i < num; i++)
        {
          e[i] = meshes[i]->get_element(id);
          static const Rect H2D_UNITY = { 0, 0, ONE, ONE };
          cr = er[i] = H2D_UNITY;
        }
        base = e[0];
        tri = base->is_triangle();
        union_recurrent(&cr, e, er, idx, unimesh->get_element(id));
      }

      delete [] e;
      delete [] er;
      delete [] idx;

      return unidata;
    }
  }
}