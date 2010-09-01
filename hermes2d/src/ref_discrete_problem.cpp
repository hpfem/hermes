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

#include "common.h"
#include "space.h"
#include "weakform.h"
#include "ref_discrete_problem.h"
#include "solution.h"

RefDiscreteProblem::RefDiscreteProblem(DiscreteProblem* base, int order_increase, 
		    int refinement) : DiscreteProblem(base->wf)
{
  this->base = base;
  if (this->base->have_spaces == false) 
    error("RefDiscreteProblem: spaces in base system not up to date.");
  
  this->solver = base->solver;

  if(this->base->get_num_dofs() <= 0) error("RefDiscreteProblem - base system has ndof = 0.");

  this->wf = this->base->wf;
  this->order_increase = order_increase;
  this->refinement = refinement;

  // perform uniform mesh refinement
  global_refinement();

  // internal check
  if (this->have_spaces == false) error("RefDiscreteProblem: missing space(s).");
}

void RefDiscreteProblem::global_refinement() 
{
  // after this, meshes and spaces are NULL
  this->free_spaces();

  // create new meshes and spaces
  this->meshes = new Mesh*[this->wf->neq];
  this->spaces = new Space*[this->wf->neq];
  this->sp_seq = new int[this->wf->neq];

  int i, j;
  // copy meshes from the coarse problem and refine them
  for (i = 0; i < this->wf->neq; i++)
  {
    Mesh* mesh = base->spaces[i]->get_mesh();

    // check if we already have the same mesh
    for (j = 0; j < i; j++)
      if (mesh->get_seq() == base->spaces[j]->get_mesh()->get_seq())
        break;

    if (j < i) // yes
    {
      meshes[i] = meshes[j];
    }
    else // no, copy and refine the coarse one
    {
      Mesh* rmesh = new Mesh;
      rmesh->copy(mesh);
      if (refinement ==  1) rmesh->refine_all_elements();
      if (refinement == -1) rmesh->unrefine_all_elements();
      this->meshes[i] = rmesh;
    }
  }

  // duplicate spaces from the coarse problem, assign reference orders and dofs
  int ndof = 0;
  for (i = 0; i < this->wf->neq; i++)
  {
    this->spaces[i] = this->base->spaces[i]->dup(this->meshes[i]);

    if (refinement == -1)
    {
      Element* re;
      for_all_active_elements(re, meshes[i])
      {
        Mesh* mesh = this->base->spaces[i]->get_mesh();
        Element* e = mesh->get_element(re->id);
        int max_order_h = 0, max_order_v = 0;
        if (e->active) {
          int quad_order = base->spaces[i]->get_element_order(e->id);
          max_order_h = H2D_GET_H_ORDER(quad_order);
          max_order_v = H2D_GET_V_ORDER(quad_order);
        }
        else
        { //find maximum order of sons
          for (int son = 0; son < 4; son++)
          {
            if (e->sons[son] != NULL)
            {
              int quad_order = base->spaces[i]->get_element_order(e->sons[son]->id);
              max_order_h = std::max(max_order_h, H2D_GET_H_ORDER(quad_order));
              max_order_v = std::max(max_order_v, H2D_GET_V_ORDER(quad_order));
            }
          }
        }

        //increase order and set it to element
        max_order_h = std::max(1, max_order_h + order_increase);
        if (re->is_triangle())
          max_order_v = 0;
        else
          max_order_v = std::max(1, max_order_v + order_increase);
        spaces[i]->set_element_order_internal(re->id, H2D_MAKE_QUAD_ORDER(max_order_h, max_order_v));
      }
    }
    else {
      this->spaces[i]->copy_orders(base->spaces[i], order_increase);
    }
  }

  this->assign_dofs();

  for (int i=0; i < this->wf->neq; i++) {
    this->spaces[i]->set_seq(this->base->spaces[i]->get_seq());
    this->sp_seq[i] = this->base->spaces[i]->get_seq();
  }

  this->pss = new PrecalcShapeset*[this->wf->neq];
  for(int i=0; i < this->wf->neq; i++) this->pss[i] = this->base->pss[i];
  this->have_spaces = true;
}

RefDiscreteProblem::~RefDiscreteProblem()
{
}

// just to prevent the user from calling this
void RefDiscreteProblem::set_spaces(Tuple<Space*> spaces)
{
  error("Method set_spaces must not be called in RefDiscreteProblem.");
}

// just to prevent the user from calling this
void RefDiscreteProblem::set_pss(Tuple<PrecalcShapeset*> pss)
{
  error("Method set_pss must not be called in RefDiscreteProblem.");
}

void RefDiscreteProblem::set_order_increase(int order_increase)
{
  this->order_increase = order_increase;
}

bool RefDiscreteProblem::solve_exact(ExactFunction exactfn, Solution* sln)
{
  Space* space = spaces[0];

  // sanity checkz
  if (!space->is_up_to_date()) error("'space' is not up to date.");

  //set mesh and function
  sln->set_exact(meshes[0], exactfn);

  return true;
}

