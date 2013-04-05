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

#include "discrete_problem/discrete_problem_form_assembler.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    bool DiscreteProblemFormAssembler<Scalar>::form_to_be_assembled(MatrixForm<Scalar>* form, Traverse::State* current_state)
    {
      if(current_state->e[form->i] && current_state->e[form->j])
      {
        if(fabs(form->scaling_factor) < 1e-12)
          return false;

        // If a block scaling table is provided, and if the scaling coefficient
        // A_mn for this block is zero, then the form does not need to be assembled.
        if(this->block_weights)
          if(fabs(this->block_weights->get_A(form->i, form->j)) < 1e-12)
            return false;
        return true;
      }
      return false;
    }

    template<typename Scalar>
    bool DiscreteProblemFormAssembler<Scalar>::form_to_be_assembled(MatrixFormVol<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((MatrixForm<Scalar>*)form, current_state))
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    bool DiscreteProblemFormAssembler<Scalar>::form_to_be_assembled(MatrixFormSurf<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((MatrixForm<Scalar>*)form, current_state))
        return false;

      if(current_state->rep->en[current_state->isurf]->marker == 0)
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->en[current_state->isurf]->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    bool DiscreteProblemFormAssembler<Scalar>::form_to_be_assembled(VectorForm<Scalar>* form, Traverse::State* current_state)
    {
      if(!current_state->e[form->i])
        return false;
      if(fabs(form->scaling_factor) < 1e-12)
        return false;

      return true;
    }

    template<typename Scalar>
    bool DiscreteProblemFormAssembler<Scalar>::form_to_be_assembled(VectorFormVol<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((VectorForm<Scalar>*)form, current_state))
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    bool DiscreteProblemFormAssembler<Scalar>::form_to_be_assembled(VectorFormSurf<Scalar>* form, Traverse::State* current_state)
    {
      if(!form_to_be_assembled((VectorForm<Scalar>*)form, current_state))
        return false;

      if(current_state->rep->en[current_state->isurf]->marker == 0)
        return false;

      if(form->assembleEverywhere)
        return true;

      int this_marker = current_state->rep->en[current_state->isurf]->marker;
      for (unsigned int ss = 0; ss < form->areas_internal.size(); ss++)
        if(form->areas_internal[ss] == this_marker)
          return true;

      return false;
    }

    template<typename Scalar>
    void DiscreteProblemFormAssembler<Scalar>::assemble_matrix_form(MatrixForm<Scalar>* form, int order, Func<double>** base_fns, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext,
      AsmList<Scalar>* current_als_i, AsmList<Scalar>* current_als_j, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights)
    {
      bool surface_form = (dynamic_cast<MatrixFormVol<Scalar>*>(form) == NULL);

      double block_scaling_coefficient = this->block_scaling_coeff(form);

      bool tra = (form->i != form->j) && (form->sym != 0);
      bool sym = (form->i == form->j) && (form->sym == 1);

      // Assemble the local stiffness matrix for the form form.
      Scalar **local_stiffness_matrix = new_matrix<Scalar>(std::max(current_als_i->cnt, current_als_j->cnt));

      Func<Scalar>** local_ext = ext;
      // If the user supplied custom ext functions for this form.
      if(form->ext.size() > 0)
      {
        int local_ext_count = form->ext.size();
        local_ext = new Func<Scalar>*[local_ext_count];
        for(int ext_i = 0; ext_i < local_ext_count; ext_i++)
          if(form->ext[ext_i])
            local_ext[ext_i] = current_state->e[ext_i] == NULL ? NULL : init_fn(form->ext[ext_i].get(), order);
          else
            local_ext[ext_i] = NULL;
      }

      // Account for the previous time level solution previously inserted at the back of ext.
      if(rungeKutta)
        u_ext += form->u_ext_offset;

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_als_i->cnt; i++)
      {
        if(current_als_i->dof[i] < 0)
          continue;

        if((!tra || surface_form) && current_als_i->dof[i] < 0)
          continue;
        if(std::abs(current_als_i->coef[i]) < 1e-12)
          continue;
        if(!sym)
        {
          for (unsigned int j = 0; j < current_als_j->cnt; j++)
          {
            if(current_als_j->dof[j] >= 0)
            {
              // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
              if(std::abs(current_als_j->coef[j]) < 1e-12)
                continue;

              Func<double>* u = base_fns[j];
              Func<double>* v = test_fns[i];

              if(surface_form)
                local_stiffness_matrix[i][j] = 0.5 * block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];
              else
                local_stiffness_matrix[i][j] = block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];
            }
          }
        }
        // Symmetric block.
        else
        {
          for (unsigned int j = 0; j < current_als_j->cnt; j++)
          {
            if(j < i && current_als_j->dof[j] >= 0)
              continue;
            if(current_als_j->dof[j] >= 0)
            {
              // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
              if(std::abs(current_als_j->coef[j]) < 1e-12)
                continue;

              Func<double>* u = base_fns[j];
              Func<double>* v = test_fns[i];

              Scalar val = block_scaling_coefficient * form->value(n_quadrature_points, jacobian_x_weights, u_ext, u, v, geometry, local_ext) * form->scaling_factor * current_als_j->coef[j] * current_als_i->coef[i];

              local_stiffness_matrix[i][j] = local_stiffness_matrix[j][i] = val;
            }
          }
        }
      }

      // Insert the local stiffness matrix into the global one.
      current_mat->add(current_als_i->cnt, current_als_j->cnt, local_stiffness_matrix, current_als_i->dof, current_als_j->dof);

      // Insert also the off-diagonal (anti-)symmetric block, if required.
      if(tra)
      {
        if(form->sym < 0)
          chsgn(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);
        transpose(local_stiffness_matrix, current_als_i->cnt, current_als_j->cnt);

        current_mat->add(current_als_j->cnt, current_als_i->cnt, local_stiffness_matrix, current_als_j->dof, current_als_i->dof);
      }

      if(form->ext.size() > 0)
      {
        for(int ext_i = 0; ext_i < form->ext.size(); ext_i++)
          if(form->ext[ext_i])
          {
            local_ext[ext_i]->free_fn();
            delete local_ext[ext_i];
          }
          delete [] local_ext;
      }

      if(rungeKutta)
        u_ext -= form->u_ext_offset;

      // Cleanup.
      delete [] local_stiffness_matrix;
    }

    template<typename Scalar>
    void DiscreteProblemFormAssembler<Scalar>::assemble_vector_form(VectorForm<Scalar>* form, int order, Func<double>** test_fns, Func<Scalar>** ext, Func<Scalar>** u_ext, 
      AsmList<Scalar>* current_als_i, Traverse::State* current_state, int n_quadrature_points, Geom<double>* geometry, double* jacobian_x_weights)
    {
      bool surface_form = (dynamic_cast<VectorFormVol<Scalar>*>(form) == NULL);

      Func<Scalar>** local_ext = ext;

      // If the user supplied custom ext functions for this form.
      if(form->ext.size() > 0)
      {
        int local_ext_count = form->ext.size();
        local_ext = new Func<Scalar>*[local_ext_count];
        for(int ext_i = 0; ext_i < local_ext_count; ext_i++)
          if(form->ext[ext_i])
            local_ext[ext_i] = init_fn(form->ext[ext_i].get(), order);
          else
            local_ext[ext_i] = NULL;
      }

      // Account for the previous time level solution previously inserted at the back of ext.
      if(rungeKutta)
        u_ext += form->u_ext_offset;

      // Actual form-specific calculation.
      for (unsigned int i = 0; i < current_als_i->cnt; i++)
      {
        if(current_als_i->dof[i] < 0)
          continue;

        // Is this necessary, i.e. is there a coefficient smaller than 1e-12?
        if(std::abs(current_als_i->coef[i]) < 1e-12)
          continue;

        Func<double>* v = test_fns[i];

        Scalar val;
        if(surface_form)
          val = 0.5 * form->value(n_quadrature_points, jacobian_x_weights, u_ext, v, geometry, local_ext) * form->scaling_factor * current_als_i->coef[i];
        else
          val = form->value(n_quadrature_points, jacobian_x_weights, u_ext, v, geometry, local_ext) * form->scaling_factor * current_als_i->coef[i];

        current_rhs->add(current_als_i->dof[i], val);
      }

      if(form->ext.size() > 0)
      {
        for(int ext_i = 0; ext_i < form->ext.size(); ext_i++)
          if(form->ext[ext_i])
          {
            local_ext[ext_i]->free_fn();
            delete local_ext[ext_i];
          }
          delete [] local_ext;
      }

      if(rungeKutta)
        u_ext -= form->u_ext_offset;
    }
    
    template class HERMES_API DiscreteProblemFormAssembler<double>;
    template class HERMES_API DiscreteProblemFormAssembler<std::complex<double> >;
  }
}