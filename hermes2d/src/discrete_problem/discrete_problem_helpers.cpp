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

#include "discrete_problem/discrete_problem_helpers.h"
#include "forms.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Mixins
    {
      template<typename Scalar>
      DiscreteProblemRungeKutta<Scalar>::DiscreteProblemRungeKutta() : rungeKutta(false), RK_original_spaces_count(0), force_diagonal_blocks(false), block_weights(nullptr)
      {

      }

      template<typename Scalar>
      void DiscreteProblemRungeKutta<Scalar>::set_RK(int original_spaces_count, bool force_diagonal_blocks_, Table* block_weights_)
      {
        this->rungeKutta = true;
        this->RK_original_spaces_count = original_spaces_count;
        this->force_diagonal_blocks = force_diagonal_blocks_;
        this->block_weights = block_weights_;
      }

      template<typename Scalar>
      double DiscreteProblemRungeKutta<Scalar>::block_scaling_coeff(MatrixForm<Scalar>* form) const
      {
        if (block_weights)
          return block_weights->get_A(form->i / this->RK_original_spaces_count, form->j / this->RK_original_spaces_count);
        return 1.0;
      }

      template<typename Scalar>
      double DiscreteProblemRungeKutta<Scalar>::block_scaling_coeff(MatrixFormDG<Scalar>* form) const
      {
        if (block_weights)
          return block_weights->get_A(form->i / this->RK_original_spaces_count, form->j / this->RK_original_spaces_count);
        return 1.0;
      }

      template<typename Scalar>
      DiscreteProblemWeakForm<Scalar>::DiscreteProblemWeakForm(WeakFormSharedPtr<Scalar> wf_) : wf(wf_)
      {
      }

      template<typename Scalar>
      void DiscreteProblemWeakForm<Scalar>::set_weak_formulation(WeakFormSharedPtr<Scalar> wf_)
      {
        if (!wf_)
          throw Hermes::Exceptions::NullException(0);

        this->wf = wf_;
      }

      template<typename Scalar>
      WeakFormSharedPtr<Scalar> DiscreteProblemWeakForm<Scalar>::get_weak_formulation() const
      {
        return this->wf;
      }

      template<typename Scalar>
      DiscreteProblemMatrixVector<Scalar>::DiscreteProblemMatrixVector() : current_mat(nullptr), current_rhs(nullptr)
      {
      }

      template<typename Scalar>
      bool DiscreteProblemMatrixVector<Scalar>::set_matrix(SparseMatrix<Scalar>* mat)
      {
        this->current_mat = mat;
        return true;
      }

      template<typename Scalar>
      bool DiscreteProblemMatrixVector<Scalar>::set_rhs(Vector<Scalar>* rhs)
      {
        this->current_rhs = rhs;
        return true;
      }

      template class HERMES_API DiscreteProblemRungeKutta<double>;
      template class HERMES_API DiscreteProblemRungeKutta<std::complex<double> >;

      template class HERMES_API DiscreteProblemWeakForm<double>;
      template class HERMES_API DiscreteProblemWeakForm<std::complex<double> >;

      template class HERMES_API DiscreteProblemMatrixVector<double>;
      template class HERMES_API DiscreteProblemMatrixVector<std::complex<double> >;
    }

    int init_geometry_points(RefMap** reference_mapping, int reference_mapping_count, int order, Geom<double>*& geometry, double*& jacobian_x_weights)
    {
      Element* rep_element = nullptr;
      RefMap* rep_reference_mapping = nullptr;
      for (int i = 0; i < reference_mapping_count; i++)
      {
        if (reference_mapping[i])
        {
          if (reference_mapping[i]->get_active_element())
          {
            rep_element = reference_mapping[i]->get_active_element();
            rep_reference_mapping = reference_mapping[i];
            break;
          }
        }
      }

      double3* pt = rep_reference_mapping->get_quad_2d()->get_points(order, rep_element->get_mode());
      int np = rep_reference_mapping->get_quad_2d()->get_num_points(order, rep_element->get_mode());

      // Init geometry and jacobian*weights.
      geometry = init_geom_vol(rep_reference_mapping, order);

      jacobian_x_weights = new double[np];
      if (rep_reference_mapping->is_jacobian_const())
      {
        double jac = rep_reference_mapping->get_const_jacobian();
        for (int i = 0; i < np; i++)
          jacobian_x_weights[i] = pt[i][2] * jac;
      }
      else
      {
        double* jac = rep_reference_mapping->get_jacobian(order);
        for (int i = 0; i < np; i++)
          jacobian_x_weights[i] = pt[i][2] * jac[i];
      }
      return np;
    }

    int init_geometry_points_allocated_jwt(RefMap* rep_reference_mapping, int order, Geom<double>*& geometry, double* jacobian_x_weights)
    {
      Element* e = rep_reference_mapping->get_active_element();
      double3* pt = rep_reference_mapping->get_quad_2d()->get_points(order, e->get_mode());
      int np = rep_reference_mapping->get_quad_2d()->get_num_points(order, e->get_mode());

      // Init geometry and jacobian*weights.
      geometry = init_geom_vol(rep_reference_mapping, order);

      if (rep_reference_mapping->is_jacobian_const())
      {
        double jac = rep_reference_mapping->get_const_jacobian();
        for (int i = 0; i < np; i++)
          jacobian_x_weights[i] = pt[i][2] * jac;
      }
      else
      {
        double* jac = rep_reference_mapping->get_jacobian(order);
        for (int i = 0; i < np; i++)
          jacobian_x_weights[i] = pt[i][2] * jac[i];
      }
      return np;
    }

    int init_surface_geometry_points(RefMap** reference_mapping, int reference_mapping_count, int& order, int isurf, int marker, Geom<double>*& geometry, double*& jacobian_x_weights)
    {
      RefMap* rep_reference_mapping = nullptr;
      for (int i = 0; i < reference_mapping_count; i++)
      {
        if (reference_mapping[i])
        {
          if (reference_mapping[i]->get_active_element())
          {
            rep_reference_mapping = reference_mapping[i];
            break;
          }
        }
      }

      int eo = rep_reference_mapping->get_quad_2d()->get_edge_points(isurf, order, rep_reference_mapping->get_active_element()->get_mode());
      double3* pt = rep_reference_mapping->get_quad_2d()->get_points(eo, rep_reference_mapping->get_active_element()->get_mode());
      int np = rep_reference_mapping->get_quad_2d()->get_num_points(eo, rep_reference_mapping->get_active_element()->get_mode());

      // Init geometry and jacobian*weights.
      double3* tan;
      geometry = init_geom_surf(rep_reference_mapping, isurf, marker, eo, tan);
      jacobian_x_weights = new double[np];
      for (int i = 0; i < np; i++)
        jacobian_x_weights[i] = pt[i][2] * tan[i][2];
      order = eo;
      return np;
    }

    int init_surface_geometry_points_allocated_jwt(RefMap* rep_reference_mapping, int& order, int isurf, int marker, Geom<double>*& geometry, double* jacobian_x_weights)
    {
      Element* e = rep_reference_mapping->get_active_element();
      int eo = rep_reference_mapping->get_quad_2d()->get_edge_points(isurf, order, e->get_mode());
      double3* pt = rep_reference_mapping->get_quad_2d()->get_points(eo, e->get_mode());
      int np = rep_reference_mapping->get_quad_2d()->get_num_points(eo, e->get_mode());

      // Init geometry and jacobian*weights.
      double3* tan;
      geometry = init_geom_surf(rep_reference_mapping, isurf, marker, eo, tan);
      for (int i = 0; i < np; i++)
        jacobian_x_weights[i] = pt[i][2] * tan[i][2];
      order = eo;
      return np;
    }
  }
}
