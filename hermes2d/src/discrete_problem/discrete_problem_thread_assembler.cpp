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
#ifdef _DEBUG
    static int cache_searches = 0;
    static int cache_record_found = 0;
    static int cache_record_found_reinit = 0;
    static int cache_record_not_found = 0;
#endif

    template<typename Scalar>
    DiscreteProblemThreadAssembler<Scalar>::DiscreteProblemThreadAssembler(DiscreteProblemSelectiveAssembler<Scalar>* selectiveAssembler) :
      pss(nullptr), refmaps(nullptr), u_ext(nullptr), als(nullptr), alsSurface(nullptr),
      selectiveAssembler(selectiveAssembler), integrationOrderCalculator(selectiveAssembler)
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

      pss = new PrecalcShapeset*[spaces_size];
      refmaps = new RefMap*[spaces_size];
      als = new AsmList<Scalar>*[spaces_size];
      alsSurface = new AsmList<Scalar>**[spaces_size];

      for (unsigned int j = 0; j < spaces_size; j++)
      {
        pss[j] = new PrecalcShapeset(spaces[j]->shapeset);
        refmaps[j] = new RefMap();
        refmaps[j]->set_quad_2d(&g_quad_2d_std);
        als[j] = new AsmList<Scalar>();
        alsSurface[j] = new AsmList<Scalar>*[H2D_MAX_NUMBER_EDGES];
        for (unsigned int k = 0; k < H2D_MAX_NUMBER_EDGES; k++)
          alsSurface[j][k] = new AsmList<Scalar>();
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
      u_ext = new Solution<Scalar>*[spaces_size];

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
      // Basic settings.
      this->nonlinear = nonlinear_;
      this->add_dirichlet_lift = add_dirichlet_lift_;

      // Transformables setup.
      fns.clear();
      for (unsigned j = 0; j < this->spaces_size; j++)
      {
        fns.push_back(pss[j]);
        pss[j]->set_quad_2d(&g_quad_2d_std);
      }
      for (unsigned j = 0; j < this->wf->ext.size(); j++)
      {
        fns.push_back(this->wf->ext[j].get());
        this->wf->ext[j]->set_quad_2d(&g_quad_2d_std);
      }

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

      // Pre-prepare calculation storage.
      this->funcs = (Func<double>***)calloc(this->spaces_size, sizeof(Func<double>**));
      this->asmlistCnt = new int[this->spaces_size];
      this->asmlistIdx = new int*[this->spaces_size];

      if (this->wf->mfsurf.size() > 0 || this->wf->vfsurf.size() > 0)
      {
        this->geometrySurface = new Geom<double>*[H2D_MAX_NUMBER_EDGES];
        this->jacobian_x_weightsSurface = new double*[H2D_MAX_NUMBER_EDGES];

        this->n_quadrature_pointsSurface = new int[H2D_MAX_NUMBER_EDGES];
        this->orderSurface = new int[H2D_MAX_NUMBER_EDGES];

        this->funcsSurface = (Func<double>****)calloc(H2D_MAX_NUMBER_EDGES, sizeof(Func<double>***));
        this->asmlistSurfaceCnt = new int*[H2D_MAX_NUMBER_EDGES];
      }
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
          spaces[j]->get_element_assembly_list(current_state->e[j], als[j]);
          refmaps[j]->set_active_element(current_state->e[j]);
          refmaps[j]->force_transform(pss[j]->get_transform(), pss[j]->get_ctm());
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
                spaces[j]->get_boundary_assembly_list(current_state->e[j], k, alsSurface[j][k]);
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
        {
          this->funcs[space_i] = nullptr;
          continue;
        }

        this->asmlistCnt[space_i] = this->als[space_i]->cnt;
        this->asmlistIdx[space_i] = new int[this->als[space_i]->cnt];
        memcpy(this->asmlistIdx[space_i], this->als[space_i]->idx, this->als[space_i]->cnt * sizeof(int));

        this->funcs[space_i] = new Func<double>*[this->asmlistCnt[space_i]];

        for (unsigned int j = 0; j < this->asmlistCnt[space_i]; j++)
        {
          pss[space_i]->set_active_shape(this->als[space_i]->idx[j]);
          this->funcs[space_i][j] = init_fn(pss[space_i], refmaps[space_i], this->order);
        }
      }

      this->n_quadrature_points = init_geometry_points(refmaps, this->spaces_size, this->order, this->geometry, this->jacobian_x_weights);

      if (current_state->isBnd && (this->wf->mfsurf.size() > 0 || this->wf->vfsurf.size() > 0))
      {
        int order_local = this->order;

        for (current_state->isurf = 0; current_state->isurf < this->current_state->rep->nvert; current_state->isurf++)
        {
          if (!current_state->bnd[current_state->isurf])
          {
            this->funcsSurface[current_state->isurf] = nullptr;
            continue;
          }

          this->n_quadrature_pointsSurface[current_state->isurf] = init_surface_geometry_points(refmaps, this->spaces_size, this->order, current_state->isurf, current_state->rep->marker, this->geometrySurface[current_state->isurf], this->jacobian_x_weightsSurface[current_state->isurf]);
          this->orderSurface[current_state->isurf] = this->order;
          this->order = order_local;

          this->funcsSurface[current_state->isurf] = (Func<double>***)calloc(this->spaces_size, sizeof(Func<double>**));

          this->asmlistSurfaceCnt[current_state->isurf] = new int[this->spaces_size];

          for (unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
          {
            if (!current_state->e[space_i])
              continue;

            unsigned int func_count = this->alsSurface[space_i][current_state->isurf]->cnt;
            this->asmlistSurfaceCnt[current_state->isurf][space_i] = func_count;

            this->funcsSurface[current_state->isurf][space_i] = new Func<double>*[func_count];
            for (unsigned int j = 0; j < func_count; j++)
            {
              pss[space_i]->set_active_shape(this->alsSurface[space_i][current_state->isurf]->idx[j]);
              this->funcsSurface[current_state->isurf][space_i][j] = init_fn(pss[space_i], refmaps[space_i], this->orderSurface[current_state->isurf]);
            }
          }
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_calculation_variables()
    {
      for (unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
      {
        if (!this->funcs[space_i])
          continue;

        for (unsigned int i = 0; i < this->asmlistCnt[space_i]; i++)
        {
          this->funcs[space_i][i]->free_fn();
          delete this->funcs[space_i][i];
        }
        delete[] this->funcs[space_i];
        delete[] this->asmlistIdx[space_i];
      }

      delete[] this->jacobian_x_weights;
      this->geometry->free();
      delete this->geometry;

      if (current_state->isBnd && (this->wf->mfsurf.size() > 0 || this->wf->vfsurf.size() > 0))
      {
        for (unsigned int edge_i = 0; edge_i < this->current_state->rep->nvert; edge_i++)
        {
          if (!this->funcsSurface[edge_i])
            continue;

          this->geometrySurface[edge_i]->free();
          delete this->geometrySurface[edge_i];
          delete[] this->jacobian_x_weightsSurface[edge_i];

          for (unsigned int space_i = 0; space_i < this->spaces_size; space_i++)
          {
            if (!this->funcsSurface[edge_i][space_i])
              continue;
            for (unsigned int i = 0; i < this->asmlistSurfaceCnt[edge_i][space_i]; i++)
            {
              this->funcsSurface[edge_i][space_i][i]->free_fn();
              delete this->funcsSurface[edge_i][space_i][i];
            }
            delete[] this->funcsSurface[edge_i][space_i];
          }
          ::free(this->funcsSurface[edge_i]);
          delete[] this->asmlistSurfaceCnt[edge_i];
        }
      }
    }

    template<typename Scalar>
    Func<Scalar>** DiscreteProblemThreadAssembler<Scalar>::init_u_ext_values(int order)
    {
      Func<Scalar>** u_ext_func = nullptr;
      if (this->nonlinear)
      {
        assert(this->u_ext);

        u_ext_func = new Func<Scalar>*[spaces_size];

        for (int i = 0; i < spaces_size; i++)
        {
          assert(u_ext[i]);

          if (u_ext[i]->get_active_element())
            u_ext_func[i] = init_fn(u_ext[i], order);
          else
            u_ext_func[i] = nullptr;
        }
      }

      return u_ext_func;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_u_ext_values(Func<Scalar>** u_ext_func)
    {
      if (u_ext_func)
      {
        for (int i = 0; i < spaces_size; i++)
        {
          if (u_ext_func[i])
          {
            u_ext_func[i]->free_fn();
            delete u_ext_func[i];
          }
        }
        delete[] u_ext_func;
      }
    }

    template<typename Scalar>
    Func<Scalar>** DiscreteProblemThreadAssembler<Scalar>::init_ext_values(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& ext, Hermes::vector<UExtFunctionSharedPtr<Scalar> >& u_ext_fns, int order, Func<Scalar>** u_ext_func, Geom<double>* geometry)
    {
      Func<Scalar>** ext_func = nullptr;
      int ext_size = ext.size();
      int u_ext_fns_size = u_ext_fns.size();

      if (ext_size > 0 || u_ext_fns_size > 0)
      {
        ext_func = new Func<Scalar>*[ext_size + u_ext_fns_size];
        for (int ext_i = 0; ext_i < u_ext_fns_size; ext_i++)
        {
          if (u_ext_fns[ext_i])
          {
            ext_func[ext_i] = init_fn(u_ext_fns[ext_i].get(), u_ext_func, this->spaces_size, order, geometry, current_state->rep->get_mode());
          }
          else
            ext_func[ext_i] = nullptr;
        }

        for (int ext_i = 0; ext_i < ext_size; ext_i++)
        {
          if (ext[ext_i])
          {
            if (ext[ext_i]->get_active_element())
              ext_func[u_ext_fns_size + ext_i] = init_fn(ext[ext_i].get(), order);
            else
              ext_func[u_ext_fns_size + ext_i] = nullptr;
          }
          else
            ext_func[u_ext_fns_size + ext_i] = nullptr;
        }
      }

      if (this->rungeKutta)
      {
        for (int ext_i = 0; ext_i < ext.size(); ext_i++)
          u_ext_func[ext_i]->add(ext_func[ext.size() - this->RK_original_spaces_count + ext_i]);
      }

      return ext_func;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_ext_values(Hermes::vector<MeshFunctionSharedPtr<Scalar> >& ext, Hermes::vector<UExtFunctionSharedPtr<Scalar> >& u_ext_fns, Func<Scalar>** ext_func)
    {
      if (ext_func)
      {
        for (int ext_i = 0; ext_i < u_ext_fns.size(); ext_i++)
        {
          ext_func[ext_i]->free_fn();
          delete ext_func[ext_i];
        }

        for (int ext_i = 0; ext_i < ext.size(); ext_i++)
        {
          if (ext[ext_i] && ext[ext_i]->get_active_element())
          {
            ext_func[u_ext_fns.size() + ext_i]->free_fn();
            delete ext_func[u_ext_fns.size() + ext_i];
          }
        }

        delete[] ext_func;
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::assemble_one_state()
    {
      // init - u_ext_func
      Func<Scalar>** u_ext_func = this->init_u_ext_values(this->order);

      // init - ext
      Func<Scalar>** ext_func = this->init_ext_values(this->wf->ext, this->wf->u_ext_fn, this->order, u_ext_func, this->geometry);

      if (this->current_mat || this->add_dirichlet_lift)
      {
        for (int current_mfvol_i = 0; current_mfvol_i < this->wf->mfvol.size(); current_mfvol_i++)
        {
          MatrixFormVol<Scalar>* mfv = this->wf->mfvol[current_mfvol_i];

          if (!selectiveAssembler->form_to_be_assembled(mfv, current_state))
            continue;

          int form_i = mfv->i;
          int form_j = mfv->j;

          this->assemble_matrix_form(mfv,
            this->order,
            this->funcs[form_j], this->funcs[form_i],
            ext_func, u_ext_func,
            this->als[form_i], this->als[form_j],
            this->n_quadrature_points, this->geometry, this->jacobian_x_weights);
        }
      }
      if (this->current_rhs)
      {
        for (int current_vfvol_i = 0; current_vfvol_i < this->wf->vfvol.size(); current_vfvol_i++)
        {
          VectorFormVol<Scalar>* vfv = this->wf->vfvol[current_vfvol_i];

          if (!selectiveAssembler->form_to_be_assembled(vfv, current_state))
            continue;

          int form_i = vfv->i;

          this->assemble_vector_form(vfv,
            this->order,
            this->funcs[form_i],
            ext_func, u_ext_func,
            this->als[form_i],
            this->n_quadrature_points, this->geometry, this->jacobian_x_weights);
        }
      }

      // deinit - u_ext_func
      this->deinit_u_ext_values(u_ext_func);

      // deinit - ext
      this->deinit_ext_values(this->wf->ext, this->wf->u_ext_fn, ext_func);

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
          Func<Scalar>** u_ext_funcSurf = this->init_u_ext_values(this->orderSurface[isurf]);

          // init - ext
          Func<Scalar>** ext_funcSurf = this->init_ext_values(this->wf->ext, this->wf->u_ext_fn, this->orderSurface[isurf], u_ext_funcSurf, this->geometrySurface[isurf]);

          if (this->current_mat || this->add_dirichlet_lift)
          {
            for (int current_mfsurf_i = 0; current_mfsurf_i < this->wf->mfsurf.size(); current_mfsurf_i++)
            {
              if (!selectiveAssembler->form_to_be_assembled(this->wf->mfsurf[current_mfsurf_i], current_state))
                continue;

              int form_i = this->wf->mfsurf[current_mfsurf_i]->i;
              int form_j = this->wf->mfsurf[current_mfsurf_i]->j;

              this->assemble_matrix_form(this->wf->mfsurf[current_mfsurf_i],
                this->orderSurface[isurf],
                this->funcsSurface[isurf][form_j],
                this->funcsSurface[isurf][form_i],
                ext_funcSurf, u_ext_funcSurf,
                this->alsSurface[form_i][isurf], this->alsSurface[form_j][isurf],
                this->n_quadrature_pointsSurface[isurf], this->geometrySurface[isurf], this->jacobian_x_weightsSurface[isurf]);
            }
          }

          if (this->current_rhs)
          {
            for (int current_vfsurf_i = 0; current_vfsurf_i < this->wf->vfsurf.size(); current_vfsurf_i++)
            {
              if (!selectiveAssembler->form_to_be_assembled(this->wf->vfsurf[current_vfsurf_i], current_state))
                continue;

              int form_i = this->wf->vfsurf[current_vfsurf_i]->i;

              this->assemble_vector_form(this->wf->vfsurf[current_vfsurf_i],
                this->orderSurface[isurf],
                this->funcsSurface[isurf][form_i],
                ext_funcSurf, u_ext_funcSurf,
                this->alsSurface[form_i][isurf],
                this->n_quadrature_pointsSurface[isurf], this->geometrySurface[isurf], this->jacobian_x_weightsSurface[isurf]);
            }
          }

          // deinit - u_ext_func
          this->deinit_u_ext_values(u_ext_funcSurf);

          // deinit - ext
          this->deinit_ext_values(this->wf->ext, this->wf->u_ext_fn, ext_funcSurf);
        }
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext,
      AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights)
    {
      bool surface_form = (dynamic_cast<MatrixFormVol<Scalar>*>(form) == nullptr);

      double block_scaling_coefficient = this->block_scaling_coeff(form);

      bool tra = (form->i != form->j) && (form->sym != 0);
      bool sym = (form->i == form->j) && (form->sym == 1);

      // Assemble the local stiffness matrix for the form form.
      Scalar **local_stiffness_matrix = new_matrix<Scalar>(std::max(current_als_i->cnt, current_als_j->cnt));

      Func<Scalar>** local_ext = ext;
      // If the user supplied custom ext functions for this form.
      if (form->ext.size() > 0)
        local_ext = this->init_ext_values(form->ext, (form->u_ext_fn.size() > 0 ? form->u_ext_fn : this->wf->u_ext_fn), order, u_ext, geometry);

      // Account for the previous time level solution previously inserted at the back of ext.
      if (this->rungeKutta)
        u_ext += form->u_ext_offset;

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

          Scalar val = block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];

          if (current_als_j->dof[j] >= 0)
          {
            if (surface_form)
              local_stiffness_matrix[i][j] = 0.5 * val;
            else
              local_stiffness_matrix[i][j] = val;

            if (sym)
              local_stiffness_matrix[j][i] = local_stiffness_matrix[i][j];
          }
          else if (this->add_dirichlet_lift && this->current_rhs)
          {
            this->current_rhs->add(current_als_i->dof[i], -val);
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
      if (this->current_mat)
        this->current_mat->add(current_als_i->cnt, current_als_j->cnt, local_stiffness_matrix, current_als_i->dof, current_als_j->dof);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if (tra)
      {
        if (form->sym < 0)
          chsgn(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);
        transpose(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);

        if (this->current_mat)
          this->current_mat->add(current_als_j->cnt, current_als_i->cnt, local_stiffness_matrix, current_als_j->dof, current_als_i->dof);

        if (this->add_dirichlet_lift && this->current_rhs)
        {
          for (unsigned int j = 0; j < current_als_i->cnt; j++)
          if (current_als_i->dof[j] < 0)
          for (unsigned int i = 0; i < current_als_j->cnt; i++)
          if (current_als_j->dof[i] >= 0)
            this->current_rhs->add(current_als_j->dof[i], -local_stiffness_matrix[i][j]);
        }
      }

      if (form->ext.size() > 0)
        this->deinit_ext_values(form->ext, (form->u_ext_fn.size() > 0 ? form->u_ext_fn : this->wf->u_ext_fn), local_ext);

      if (this->rungeKutta)
        u_ext -= form->u_ext_offset;

      // Cleanup.
      delete[] local_stiffness_matrix;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext,
      AsmList<Scalar>* current_als_i, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights)
    {
      bool surface_form = (dynamic_cast<VectorFormVol<Scalar>*>(form) == nullptr);

      Func<Scalar>** local_ext = ext;
      if (form->ext.size() > 0)
        local_ext = this->init_ext_values(form->ext, (form->u_ext_fn.size() > 0 ? form->u_ext_fn : this->wf->u_ext_fn), order, u_ext, geometry);

      // Account for the previous time level solution previously inserted at the back of ext.
      if (this->rungeKutta)
        u_ext += form->u_ext_offset;

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
          val = 0.5 * form->value(n_quadrature_points, jacobian_x_weights, u_ext, v, geometry, local_ext) * form->scaling_factor * current_als_i->coef[i];
        else
          val = form->value(n_quadrature_points, jacobian_x_weights, u_ext, v, geometry, local_ext) * form->scaling_factor * current_als_i->coef[i];

        this->current_rhs->add(current_als_i->dof[i], val);
      }

      if (form->ext.size() > 0)
        this->deinit_ext_values(form->ext, (form->u_ext_fn.size() > 0 ? form->u_ext_fn : this->wf->u_ext_fn), local_ext);

      if (this->rungeKutta)
        u_ext -= form->u_ext_offset;
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_assembling_one_state()
    {
      this->deinit_calculation_variables();
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::deinit_assembling()
    {
      ::free(this->funcs);
      delete[] this->asmlistCnt;
      delete[] this->asmlistIdx;

      this->funcsSurface = nullptr;

      if (this->wf->mfsurf.size() > 0 || this->wf->vfsurf.size() > 0)
      {
        ::free(this->funcsSurface);
        delete[] this->geometrySurface;
        delete[] this->jacobian_x_weightsSurface;

        delete[] this->n_quadrature_pointsSurface;
        delete[] this->orderSurface;
        delete[] this->asmlistSurfaceCnt;
      }
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::free()
    {
      this->free_spaces();
      this->free_weak_formulation();
      this->free_u_ext();
    }

    template<typename Scalar>
    void DiscreteProblemThreadAssembler<Scalar>::free_spaces()
    {
      if (!this->pss)
        return;

      for (unsigned int j = 0; j < spaces_size; j++)
        delete pss[j];
      delete[] pss;
      pss = nullptr;

      for (unsigned int j = 0; j < spaces_size; j++)
        delete refmaps[j];
      delete[] refmaps;

      for (unsigned int j = 0; j < spaces_size; j++)
        delete als[j];
      delete[] als;

      for (unsigned int j = 0; j < spaces_size; j++)
      {
        for (unsigned int k = 0; k < H2D_MAX_NUMBER_EDGES; k++)
          delete alsSurface[j][k];
        delete[] alsSurface[j];
      }
      delete[] alsSurface;
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
        delete[] u_ext;
      }
    }

    template class HERMES_API DiscreteProblemThreadAssembler<double>;
    template class HERMES_API DiscreteProblemThreadAssembler<std::complex<double> >;
  }
}