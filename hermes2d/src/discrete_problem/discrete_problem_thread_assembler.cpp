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

#include "discrete_problem/discrete_problem_thread_assembler.h"
#include "discrete_problem/discrete_problem_selective_assembler.h"
#include "shapeset/precalc.h"
#include "function/solution.h"
#include "weakform/weakform.h"
#include "function/exact_solution.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    DiscreteProblemThreadAssembler<Scalar>::DiscreteProblemThreadAssembler(DiscreteProblemSelectiveAssembler<Scalar>* selectiveAssembler) :
      pss(nullptr), refmaps(nullptr), u_ext(nullptr),
      selectiveAssembler(selectiveAssembler), integrationOrderCalculator(selectiveAssembler),
      ext_funcs(nullptr), ext_funcs_allocated_size(0), ext_funcs_local(nullptr), ext_funcs_local_allocated_size(0)
    {
    }

    template<typename Scalar>
    DiscreteProblemThreadAssembler<Scalar>::~DiscreteProblemThreadAssembler()
    {
      this->free();
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_spaces(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces)
    {
      this->free_spaces();

      this->spaces_size = spaces.size();

      pss = malloc_with_check<PrecalcShapeset*>(spaces_size);
      refmaps = malloc_with_check<RefMap*>(spaces_size);

      for (unsigned int j = 0; j < spaces_size; j++)
      {
        pss[j] = new PrecalcShapeset(spaces[j]->shapeset);
        refmaps[j] = new RefMap();
        refmaps[j]->set_quad_2d(&g_quad_2d_std);
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::set_weak_formulation(WeakForm<Scalar>* wf_)
    {
      this->free_weak_formulation();

      this->wf = wf_->clone();
      this->wf->cloneMembers(wf_);
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_u_ext(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Solution<Scalar>** u_ext_sln)
    {
      assert(this->spaces_size == spaces.size() && this->pss);

      free_u_ext();
      u_ext = malloc_with_check<Solution<Scalar>*>(spaces_size);

      for (unsigned int j = 0; j < spaces_size; j++)
      {
        if (u_ext_sln)
        {
          u_ext[j] = new Solution<Scalar>(spaces[j]->get_mesh());
          u_ext[j]->copy(u_ext_sln[j]);
        }
        else
        {
          if (spaces[j]->get_shapeset()->get_num_components() == 1)
            u_ext[j] = new ZeroSolution<Scalar>(spaces[j]->get_mesh());
          else
            u_ext[j] = new ZeroSolutionVector<Scalar>(spaces[j]->get_mesh());
        }
      }

      this->integrationOrderCalculator.u_ext = this->u_ext;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_assembling(Solution<Scalar>** u_ext_sln, const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, bool nonlinear_, bool add_dirichlet_lift_)
    {
      // Init the memory pool - if PJLIB is linked, it will do the magic, if not, it will initialize the pointer to null.
      this->init_funcs_memory_pool();

      // Basic settings.
      this->nonlinear = nonlinear_;
      this->add_dirichlet_lift = add_dirichlet_lift_;

      // Transformables setup.
      fns.clear();
      // - precalc shapesets.
      for (unsigned j = 0; j < this->spaces_size; j++)
      {
        fns.push_back(pss[j]);
        pss[j]->set_quad_2d(&g_quad_2d_std);
      }
      // - wf->ext.
      for (unsigned j = 0; j < this->wf->ext.size(); j++)
      {
        fns.push_back(this->wf->ext[j].get());
        this->wf->ext[j]->set_quad_2d(&g_quad_2d_std);
      }
      // - forms->ext.
      for (unsigned int form_i = 0; form_i < this->wf->get_forms().size(); form_i++)
      {
        Form<Scalar>* form = this->wf->get_forms()[form_i];
        for (unsigned int ext_i = 0; ext_i < form->ext.size(); ext_i++)
        {
          if (form->ext[ext_i])
          {
            fns.push_back(form->ext[ext_i].get());
            form->ext[ext_i]->set_quad_2d(&g_quad_2d_std);
          }
        }
      }
      // - u_ext.
      if (this->nonlinear)
      {
        init_u_ext(spaces, u_ext_sln);
        for (unsigned j = 0; j < this->wf->get_neq(); j++)
        {
          fns.push_back(u_ext[j]);
          u_ext[j]->set_quad_2d(&g_quad_2d_std);
        }
      }

      // Process markers.
      this->wf->processFormMarkers(spaces);

      // Initialize Func storage.
      this->init_funcs();
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_funcs_memory_pool()
    {
#ifdef WITH_PJLIB
      pj_thread_desc rtpdesc;
      pj_thread_t *thread;
      pj_bzero(rtpdesc, sizeof(rtpdesc));
      if (!pj_thread_is_registered())
        pj_thread_register(NULL, rtpdesc, &thread);

      this->FuncMemoryPool = pj_pool_create(&Hermes2DMemoryPoolCache.factory, // the factory
        "pool-DiscreteProblemThreadAssembler", // pool's name
        sizeof(Func<Scalar>) * H2D_MAX_LOCAL_BASIS_SIZE * (H2D_MAX_NUMBER_EDGES + 1), // initial size
        sizeof(Func<Scalar>) * H2D_MAX_LOCAL_BASIS_SIZE * H2D_MAX_NUMBER_EDGES, // increment size
        NULL);
#else
      this->FuncMemoryPool = nullptr;
#endif
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_funcs()
    {
      // Basis & test fns, u_ext funcs.
      for (unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
      {
        for (unsigned int j = 0; j < H2D_MAX_LOCAL_BASIS_SIZE; j++)
          this->funcs[space_i][j] = preallocate_fn<double>(this->FuncMemoryPool);

        for (int edge_i = 0; edge_i < H2D_MAX_NUMBER_EDGES; edge_i++)
        for (unsigned int j = 0; j < H2D_MAX_LOCAL_BASIS_SIZE; j++)
          this->funcsSurface[edge_i][space_i][j] = preallocate_fn<double>(this->FuncMemoryPool);

        if (this->nonlinear)
          this->u_ext_funcs[space_i] = preallocate_fn<Scalar>(this->FuncMemoryPool);
      }

      // Reallocation of wf-(nonlocal-) ext funcs.
      int ext_size = this->wf->ext.size();
      int u_ext_fns_size = this->wf->u_ext_fn.size();

      if (ext_size + u_ext_fns_size > ext_funcs_allocated_size)
      {
        ext_funcs_allocated_size = ext_size + u_ext_fns_size;
        ext_funcs = realloc_with_check<DiscreteProblemThreadAssembler<Scalar>, Func<Scalar>*>(ext_funcs, ext_funcs_allocated_size, this);
      }

      // Initializaton of wf-(nonlocal-)ext funcs
      if (ext_size > 0 || u_ext_fns_size > 0)
      {
        for (int ext_i = 0; ext_i < u_ext_fns_size; ext_i++)
          this->ext_funcs[ext_i] = preallocate_fn<Scalar>(this->FuncMemoryPool);

        for (int ext_i = 0; ext_i < ext_size; ext_i++)
          this->ext_funcs[u_ext_fns_size + ext_i] = preallocate_fn<Scalar>(this->FuncMemoryPool);
      }

      // Calculating local sizes.
      int local_ext_size = 0;
      int local_u_ext_fns_size = 0;
      for (int form_i = 0; form_i < this->wf->forms.size(); form_i++)
      {
        if (this->wf->forms[form_i]->ext.size() > local_ext_size)
          local_ext_size = this->wf->forms[form_i]->ext.size();

        if (this->wf->forms[form_i]->u_ext_fn.size() > local_u_ext_fns_size)
          local_u_ext_fns_size = this->wf->forms[form_i]->u_ext_fn.size();
      }

      if (local_ext_size > 0 || local_u_ext_fns_size > 0)
      {
        // Reallocation of form-(local-) ext funcs.
        if (local_ext_size + local_u_ext_fns_size > ext_funcs_local_allocated_size)
        {
          ext_funcs_local_allocated_size = local_ext_size + local_u_ext_fns_size;
          ext_funcs_local = realloc_with_check<DiscreteProblemThreadAssembler<Scalar>, Func<Scalar>*>(ext_funcs_local, ext_funcs_local_allocated_size, this);
        }

        // Initializaton of form-(local-)ext funcs
        for (int ext_i = 0; ext_i < local_u_ext_fns_size; ext_i++)
          this->ext_funcs_local[ext_i] = preallocate_fn<Scalar>(this->FuncMemoryPool);

        for (int ext_i = 0; ext_i < local_ext_size; ext_i++)
          this->ext_funcs_local[local_u_ext_fns_size + ext_i] = preallocate_fn<Scalar>(this->FuncMemoryPool);
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_funcs()
    {
#ifdef WITH_PJLIB
      pj_pool_release(this->FuncMemoryPool);
#else

      for (unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
      {
        // Test functions
        for (unsigned int j = 0; j < H2D_MAX_LOCAL_BASIS_SIZE; j++)
          delete this->funcs[space_i][j];

        // Test functions - surface
        for (int edge_i = 0; edge_i < H2D_MAX_NUMBER_EDGES; edge_i++)
        {
          for (unsigned int j = 0; j < H2D_MAX_LOCAL_BASIS_SIZE; j++)
            delete this->funcsSurface[edge_i][space_i][j];
        }

        // UExt
        if (this->nonlinear)
          delete this->u_ext_funcs[space_i];
      }

      // Ext
      int ext_size = this->wf->ext.size();
      int u_ext_fns_size = this->wf->u_ext_fn.size();
      if (ext_size > 0 || u_ext_fns_size > 0)
      {
        for (int ext_i = 0; ext_i < u_ext_fns_size; ext_i++)
          delete this->ext_funcs[ext_i];

        for (int ext_i = 0; ext_i < ext_size; ext_i++)
          delete this->ext_funcs[u_ext_fns_size + ext_i];
      }

      // Ext - local
      int local_ext_size = 0;
      int local_u_ext_fns_size = 0;
      for (int form_i = 0; form_i < this->wf->forms.size(); form_i++)
      {
        if (this->wf->forms[form_i]->ext.size() > local_ext_size)
          local_ext_size = this->wf->forms[form_i]->ext.size();

        if (this->wf->forms[form_i]->u_ext_fn.size() > local_u_ext_fns_size)
          local_u_ext_fns_size = this->wf->forms[form_i]->u_ext_fn.size();
      }

      if (local_ext_size > 0 || local_u_ext_fns_size > 0)
      {
        for (int ext_i = 0; ext_i < local_u_ext_fns_size; ext_i++)
          delete this->ext_funcs_local[ext_i];

        for (int ext_i = 0; ext_i < local_ext_size; ext_i++)
          delete this->ext_funcs_local[local_u_ext_fns_size + ext_i];
      }
#endif
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_assembling_one_state(const Hermes::vector<SpaceSharedPtr<Scalar> >& spaces, Traverse::State* current_state_)
    {
      current_state = current_state_;
      this->integrationOrderCalculator.current_state = this->current_state;

      // Active elements.
      for (int j = 0; j < fns.size(); j++)
      {
        if (current_state->e[j])
        {
          fns[j]->set_active_element(current_state->e[j]);
          fns[j]->set_transform(current_state->sub_idx[j]);
        }
      }

      // Assembly lists & refmaps.
      for (int j = 0; j < this->spaces_size; j++)
      {
        if (current_state->e[j])
        {
          spaces[j]->get_element_assembly_list(current_state->e[j], &als[j]);
          refmaps[j]->set_active_element(current_state->e[j]);
          refmaps[j]->force_transform(pss[j]->get_transform(), pss[j]->get_ctm());
          rep_refmap = refmaps[j];
        }
      }

      // Boundary assembly lists
      if (current_state->isBnd && !(this->wf->mfsurf.empty() && this->wf->vfsurf.empty()))
      {
        for (int j = 0; j < this->spaces_size; j++)
        {
          if (current_state->e[j])
          {
            for (int k = 0; k < current_state->rep->nvert; k++)
            {
              if (current_state->bnd[k])
                spaces[j]->get_boundary_assembly_list(current_state->e[j], k, &alsSurface[k][j]);
            }
          }
        }
      }

      // Volumetric integration order.
      this->order = this->integrationOrderCalculator.calculate_order(spaces, this->refmaps, this->wf);

      // Init the variables (funcs, geometry, ...)
      this->init_calculation_variables();
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_calculation_variables()
    {
      for (unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
      {
        if (current_state->e[space_i] == nullptr)
          continue;

        for (unsigned int j = 0; j < this->als[space_i].cnt; j++)
        {
          pss[space_i]->set_active_shape(this->als[space_i].idx[j]);
          init_fn_preallocated(this->funcs[space_i][j], pss[space_i], refmaps[space_i], this->order);
        }
      }

      this->n_quadrature_points = init_geometry_points_allocated_jwt(this->rep_refmap, this->order, this->geometry, this->jacobian_x_weights);

      if (current_state->isBnd && (this->wf->mfsurf.size() > 0 || this->wf->vfsurf.size() > 0))
      {
        int order_local = this->order;

        for (unsigned int edge_i = 0; edge_i < this->current_state->rep->nvert; edge_i++)
        {
          if (!current_state->bnd[edge_i])
            continue;

          this->n_quadrature_pointsSurface[edge_i] = init_surface_geometry_points_allocated_jwt(this->rep_refmap, this->order, edge_i, current_state->rep->marker, this->geometrySurface[edge_i], this->jacobian_x_weightsSurface[edge_i]);
          this->orderSurface[edge_i] = this->order;
          this->order = order_local;

          for (unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
          {
            if (!current_state->e[space_i])
              continue;

            for (unsigned int j = 0; j < this->alsSurface[edge_i][space_i].cnt; j++)
            {
              pss[space_i]->set_active_shape(this->alsSurface[edge_i][space_i].idx[j]);
              init_fn_preallocated(this->funcsSurface[edge_i][space_i][j], pss[space_i], refmaps[space_i], this->orderSurface[edge_i]);
            }
          }
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_calculation_variables()
    {
      this->geometry->free();
      delete this->geometry;

      if (current_state->isBnd && (this->wf->mfsurf.size() > 0 || this->wf->vfsurf.size() > 0))
      {
        for (unsigned int edge_i = 0; edge_i < this->current_state->rep->nvert; edge_i++)
        {
          if (!current_state->bnd[edge_i])
            continue;

          this->geometrySurface[edge_i]->free();
          delete this->geometrySurface[edge_i];
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_u_ext_values(int order)
    {
      if (this->nonlinear)
      {
        for (int i = 0; i < spaces_size; i++)
        {
          if (u_ext[i]->get_active_element())
            init_fn_preallocated(u_ext_funcs[i], u_ext[i], order);
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::init_ext_values(Func<Scalar>** target_array, Hermes::vector<MeshFunctionSharedPtr<Scalar> >& ext, Hermes::vector<UExtFunctionSharedPtr<Scalar> >& u_ext_fns, int order, Func<Scalar>** u_ext_func, Geom<double>* geometry)
    {
      int ext_size = ext.size();
      int u_ext_fns_size = u_ext_fns.size();

      if (ext_size > 0 || u_ext_fns_size > 0)
      {
        for (int ext_i = 0; ext_i < ext_size; ext_i++)
        {
          if (ext[ext_i])
          {
            if (ext[ext_i]->get_active_element())
              init_fn_preallocated(target_array[u_ext_fns_size + ext_i], ext[ext_i].get(), order);
          }
        }
        
        for (int ext_i = 0; ext_i < u_ext_fns_size; ext_i++)
        {
          if (u_ext_fns[ext_i])
            init_fn_preallocated(target_array[ext_i], u_ext_fns[ext_i].get(), target_array, u_ext_func, order, geometry, current_state->rep->get_mode());
        }
      }

      if (this->rungeKutta)
      {
        for (int ext_i = 0; ext_i < ext.size(); ext_i++)
          u_ext_func[ext_i]->add(target_array[ext.size() - this->RK_original_spaces_count + ext_i]);
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::assemble_one_state()
    {
      // init - u_ext_func
      this->init_u_ext_values(this->order);

      // init - ext
      this->init_ext_values(this->ext_funcs, this->wf->ext, this->wf->u_ext_fn, this->order, this->u_ext_funcs, this->geometry);

      if (this->current_mat || this->add_dirichlet_lift)
      {
        for (int current_mfvol_i = 0; current_mfvol_i < this->wf->mfvol.size(); current_mfvol_i++)
        {
          if (!selectiveAssembler->form_to_be_assembled(this->wf->mfvol[current_mfvol_i], current_state))
            continue;

          int form_i = this->wf->mfvol[current_mfvol_i]->i;
          int form_j = this->wf->mfvol[current_mfvol_i]->j;

          this->assemble_matrix_form(this->wf->mfvol[current_mfvol_i], order, funcs[form_j], funcs[form_i], &als[form_i], &als[form_j], n_quadrature_points, geometry, jacobian_x_weights);
        }
      }
      if (this->current_rhs)
      {
        for (int current_vfvol_i = 0; current_vfvol_i < this->wf->vfvol.size(); current_vfvol_i++)
        {
          if (!selectiveAssembler->form_to_be_assembled(this->wf->vfvol[current_vfvol_i], current_state))
            continue;

          int form_i = this->wf->vfvol[current_vfvol_i]->i;

          this->assemble_vector_form(this->wf->vfvol[current_vfvol_i], order, funcs[form_i], &als[form_i], n_quadrature_points, geometry, jacobian_x_weights);
        }
      }

      // Assemble surface integrals now: loop through surfaces of the element.
      if (current_state->isBnd && (this->wf->mfsurf.size() > 0 || this->wf->vfsurf.size() > 0))
      {
        for (int isurf = 0; isurf < current_state->rep->nvert; isurf++)
        {
          if (!current_state->bnd[isurf])
            continue;

          current_state->isurf = isurf;

          // Edge-wise parameters for WeakForm.
          this->wf->set_active_edge_state(current_state->e, isurf);

          // init - u_ext_func
          this->init_u_ext_values(this->orderSurface[isurf]);

          // init - ext
          this->init_ext_values(this->ext_funcs, this->wf->ext, this->wf->u_ext_fn, this->orderSurface[isurf], this->u_ext_funcs, this->geometrySurface[isurf]);

          if (this->current_mat || this->add_dirichlet_lift)
          {
            for (int current_mfsurf_i = 0; current_mfsurf_i < this->wf->mfsurf.size(); current_mfsurf_i++)
            {
              if (!selectiveAssembler->form_to_be_assembled(this->wf->mfsurf[current_mfsurf_i], current_state))
                continue;

              int form_i = this->wf->mfsurf[current_mfsurf_i]->i;
              int form_j = this->wf->mfsurf[current_mfsurf_i]->j;

              this->assemble_matrix_form(this->wf->mfsurf[current_mfsurf_i], orderSurface[isurf], funcsSurface[isurf][form_j], funcsSurface[isurf][form_i],
                &alsSurface[isurf][form_i], &alsSurface[isurf][form_j], n_quadrature_pointsSurface[isurf], geometrySurface[isurf], jacobian_x_weightsSurface[isurf]);
            }
          }

          if (this->current_rhs)
          {
            for (int current_vfsurf_i = 0; current_vfsurf_i < this->wf->vfsurf.size(); current_vfsurf_i++)
            {
              if (!selectiveAssembler->form_to_be_assembled(this->wf->vfsurf[current_vfsurf_i], current_state))
                continue;

              int form_i = this->wf->vfsurf[current_vfsurf_i]->i;

              this->assemble_vector_form(this->wf->vfsurf[current_vfsurf_i], orderSurface[isurf], funcsSurface[isurf][form_i], &alsSurface[isurf][form_i],
                n_quadrature_pointsSurface[isurf], geometrySurface[isurf], jacobian_x_weightsSurface[isurf]);
            }
          }
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns,
      AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights)
    {
      bool surface_form = (dynamic_cast<MatrixFormVol<Scalar>*>(form) == nullptr);

      double block_scaling_coefficient = this->block_scaling_coeff(form);

      bool tra = (form->i != form->j) && (form->sym != 0);
      bool sym = (form->i == form->j) && (form->sym == 1);

      Func<Scalar>** ext_local = this->ext_funcs;
      // If the user supplied custom ext functions for this form.
      if (form->ext.size() > 0 || form->u_ext_fn.size() > 0)
      {
        this->init_ext_values(this->ext_funcs_local, form->ext, (form->u_ext_fn.size() > 0 ? form->u_ext_fn : this->wf->u_ext_fn), order, this->u_ext_funcs, geometry);
        ext_local = this->ext_funcs_local;
      }

      // Account for the previous time level solution previously inserted at the back of ext.
      Func<Scalar>** u_ext_local = this->u_ext_funcs;
      if (this->rungeKutta)
        u_ext_local += form->u_ext_offset;

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_als_i->cnt; i++)
      {
        if (current_als_i->dof[i] < 0 || std::abs(current_als_i->coef[i]) < Hermes::HermesSqrtEpsilon)
          continue;

        if ((!tra || surface_form) && current_als_i->dof[i] < 0)
          continue;

        for (unsigned int j = 0; j < current_als_j->cnt; j++)
        {
          // Skip symmetric values that do not contribute to Dirichlet lift.
          if (sym && j < i && current_als_j->dof[j] >= 0)
            continue;

          // Skip anything that does not contribute to Dirichlet in the case of just rhs assembling.
          if (current_als_j->dof[j] >= 0 && !this->current_mat)
            continue;

          if (std::abs(current_als_j->coef[j]) < Hermes::HermesSqrtEpsilon)
            continue;

          Func<double>* u = base_fns[j];
          Func<double>* v = test_fns[i];

          Scalar val = block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext_local, u, v, geometry, ext_local) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];

          if (current_als_j->dof[j] >= 0)
          {
            int local_matrix_index_array = i * H2D_MAX_LOCAL_BASIS_SIZE + j;

            if (surface_form)
              local_stiffness_matrix[local_matrix_index_array] = 0.5 * val;
            else
              local_stiffness_matrix[local_matrix_index_array] = val;

            if (sym)
            {
              int local_matrix_index_array_transposed = j * H2D_MAX_LOCAL_BASIS_SIZE + i;
              local_stiffness_matrix[local_matrix_index_array_transposed] = local_stiffness_matrix[local_matrix_index_array];
            }
          }
          else if (this->add_dirichlet_lift && this->current_rhs)
          {
            this->current_rhs->add(current_als_i->dof[i], -val);
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
      if (this->current_mat)
        this->current_mat->add(current_als_i->cnt, current_als_j->cnt, local_stiffness_matrix, current_als_i->dof, current_als_j->dof, H2D_MAX_LOCAL_BASIS_SIZE);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if (tra)
      {
        if (form->sym < 0)
          change_sign(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt, H2D_MAX_LOCAL_BASIS_SIZE);
        transpose(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt, H2D_MAX_LOCAL_BASIS_SIZE);

        if (this->current_mat)
          this->current_mat->add(current_als_j->cnt, current_als_i->cnt, local_stiffness_matrix, current_als_j->dof, current_als_i->dof, H2D_MAX_LOCAL_BASIS_SIZE);

        if (this->add_dirichlet_lift && this->current_rhs)
        {
          for (unsigned int j = 0; j < current_als_i->cnt; j++)
          {
            if (current_als_i->dof[j] < 0)
            {
              for (unsigned int i = 0; i < current_als_j->cnt; i++)
              {
                if (current_als_j->dof[i] >= 0)
                {
                  int local_matrix_index_array = i * H2D_MAX_LOCAL_BASIS_SIZE + j;
                  this->current_rhs->add(current_als_j->dof[i], -local_stiffness_matrix[local_matrix_index_array]);
                }
              }
            }
          }
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns,
      AsmList<Scalar>* current_als_i, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights)
    {
      bool surface_form = (dynamic_cast<VectorFormVol<Scalar>*>(form) == nullptr);

      Func<Scalar>** ext_local = this->ext_funcs;
      // If the user supplied custom ext functions for this form.
      if (form->ext.size() > 0 || form->u_ext_fn.size() > 0)
      {
        this->init_ext_values(this->ext_funcs_local, form->ext, (form->u_ext_fn.size() > 0 ? form->u_ext_fn : this->wf->u_ext_fn), order, this->u_ext_funcs, geometry);
        ext_local = this->ext_funcs_local;
      }

      // Account for the previous time level solution previously inserted at the back of ext.
      Func<Scalar>** u_ext_local = this->u_ext_funcs;
      if (this->rungeKutta)
        u_ext_local += form->u_ext_offset;

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_als_i->cnt; i++)
      {
        if (current_als_i->dof[i] < 0)
          continue;

        // Is this necessary, i.e. is there a coefficient smaller than Hermes::HermesSqrtEpsilon?
        if (std::abs(current_als_i->coef[i]) < Hermes::HermesSqrtEpsilon)
          continue;

        Func<double>* v = test_fns[i];

        Scalar val;
        if (surface_form)
          val = 0.5 * form->value(n_quadrature_points, jacobian_x_weights, u_ext_local, v, geometry, ext_local) * form->scaling_factor * current_als_i->coef[i];
        else
          val = form->value(n_quadrature_points, jacobian_x_weights, u_ext_local, v, geometry, ext_local) * form->scaling_factor * current_als_i->coef[i];

        this->current_rhs->add(current_als_i->dof[i], val);
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_assembling_one_state()
    {
      this->deinit_calculation_variables();
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_assembling()
    {
      this->deinit_funcs();
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::free()
    {
      this->free_spaces();
      this->free_weak_formulation();
      this->free_u_ext();

      free_with_check(ext_funcs, true);
      free_with_check(ext_funcs_local, true);
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::free_spaces()
    {
      if (!this->pss)
        return;

      for (unsigned int j = 0; j < spaces_size; j++)
        delete pss[j];
      free_with_check(pss);

      for (unsigned int j = 0; j < spaces_size; j++)
        delete refmaps[j];
      free_with_check(refmaps);
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::free_weak_formulation()
    {
      if (this->wf)
      {
        this->wf->free_ext();
        delete this->wf;
        this->wf = nullptr;
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::free_u_ext()
    {
      if (u_ext)
      {
        for (unsigned int j = 0; j < spaces_size; j++)
          delete u_ext[j];
        free_with_check(u_ext);
      }
    }

    template class HERMES_API DiscreteProblemThreadAssembler<double>;
    template class HERMES_API DiscreteProblemThreadAssembler<std::complex<double> >;
  }
}