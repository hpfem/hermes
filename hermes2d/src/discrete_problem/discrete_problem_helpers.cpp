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
      DiscreteProblemRungeKutta<Scalar>::DiscreteProblemRungeKutta() : rungeKutta(false), RK_original_spaces_count(0), force_diagonal_blocks(false), block_weights(NULL)
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
        if(block_weights)
          return block_weights->get_A(form->i, form->j);
        return 1.0;
      }

      template<typename Scalar>
      double DiscreteProblemRungeKutta<Scalar>::block_scaling_coeff(MatrixFormDG<Scalar>* form) const
      {
        if(block_weights)
          return block_weights->get_A(form->i, form->j);
        return 1.0;
      }

      template<typename Scalar>
      DiscreteProblemWeakForm<Scalar>::DiscreteProblemWeakForm(WeakForm<Scalar>* wf_) : wf(wf_)
      {
      }

      template<typename Scalar>
      void DiscreteProblemWeakForm<Scalar>::set_weak_formulation(WeakForm<Scalar>* wf_)
      {
        if(!wf_)
          throw Hermes::Exceptions::NullException(0);

        this->wf = wf_;
      }

      template<typename Scalar>
      WeakForm<Scalar>* DiscreteProblemWeakForm<Scalar>::get_weak_formulation() const
      {
        return this->wf;
      }

      DiscreteProblemCacheSettings::DiscreteProblemCacheSettings() : 
        do_not_use_cache(false),
        report_cache_hits_and_misses(false),
        cache_searches(0),
        cache_record_found(0),
        cache_record_found_reinit(0),
        cache_record_not_found(0)
      {
      }

      void DiscreteProblemCacheSettings::set_do_not_use_cache(bool to_set)
      {
        this->do_not_use_cache = to_set;
      }

      void DiscreteProblemCacheSettings::set_report_cache_hits_and_misses(bool to_set)
      {
        this->report_cache_hits_and_misses = to_set;
      }

      void DiscreteProblemCacheSettings::get_cache_hits_and_misses(int& cache_searches_, int& cache_record_found_, int& cache_record_found_reinit_, int& cache_record_not_found_)
      {
        if(!this->report_cache_hits_and_misses)
          throw Exceptions::Exception("Asked for cache hits and misses, without turning on the calculation.");

        cache_searches_ = this->cache_searches;
        cache_record_found_ = this->cache_record_found;
        cache_record_found_reinit_ = this->cache_record_found_reinit;
        cache_record_not_found_ = this->cache_searches;
      }
      
      void DiscreteProblemCacheSettings::add_cache_hits_and_misses(DiscreteProblemCacheSettings* other)
      {
        this->cache_searches += other->cache_searches;
        this->cache_record_found += other->cache_record_found;
        this->cache_record_found_reinit += other->cache_record_found_reinit;
        this->cache_record_not_found += other->cache_searches;
      }

      void DiscreteProblemCacheSettings::zero_cache_hits_and_misses()
      {
        cache_searches = 0;
        cache_record_found = 0;
        cache_record_found_reinit = 0;
        cache_record_not_found = 0;
      }

      template<typename Scalar>
      DiscreteProblemMatrixVector<Scalar>::DiscreteProblemMatrixVector() : current_mat(NULL), current_rhs(NULL)
      {
      }

      template<typename Scalar>
      void DiscreteProblemMatrixVector<Scalar>::set_matrix(SparseMatrix<Scalar>* mat)
      {
        this->current_mat = mat;
      }

      template<typename Scalar>
      void DiscreteProblemMatrixVector<Scalar>::set_rhs(Vector<Scalar>* rhs)
      {
        this->current_rhs = rhs;
      }

      template class HERMES_API DiscreteProblemRungeKutta<double>;
      template class HERMES_API DiscreteProblemRungeKutta<std::complex<double> >;

      template class HERMES_API DiscreteProblemWeakForm<double>;
      template class HERMES_API DiscreteProblemWeakForm<std::complex<double> >;

      template class HERMES_API DiscreteProblemMatrixVector<double>;
      template class HERMES_API DiscreteProblemMatrixVector<std::complex<double> >;
    }

    int init_geometry_points(RefMap* reference_mapping, int order, Geom<double>*& geometry, double*& jacobian_x_weights)
    {
      double3* pt = reference_mapping->get_quad_2d()->get_points(order, reference_mapping->get_active_element()->get_mode());
      int np = reference_mapping->get_quad_2d()->get_num_points(order, reference_mapping->get_active_element()->get_mode());

      // Init geometry and jacobian*weights.
      geometry = init_geom_vol(reference_mapping, order);
      double* jac = NULL;
      if(!reference_mapping->is_jacobian_const())
        jac = reference_mapping->get_jacobian(order);
      jacobian_x_weights = new double[np];
      for(int i = 0; i < np; i++)
      {
        if(reference_mapping->is_jacobian_const())
          jacobian_x_weights[i] = pt[i][2] * reference_mapping->get_const_jacobian();
        else
          jacobian_x_weights[i] = pt[i][2] * jac[i];
      }
      return np;
    }

    int init_surface_geometry_points(RefMap* reference_mapping, int& order, int isurf, int marker, Geom<double>*& geometry, double*& jacobian_x_weights)
    {
      int eo = reference_mapping->get_quad_2d()->get_edge_points(isurf, order, reference_mapping->get_active_element()->get_mode());
      double3* pt = reference_mapping->get_quad_2d()->get_points(eo, reference_mapping->get_active_element()->get_mode());
      int np = reference_mapping->get_quad_2d()->get_num_points(eo, reference_mapping->get_active_element()->get_mode());

      // Init geometry and jacobian*weights.
      double3* tan;
      geometry = init_geom_surf(reference_mapping, isurf, marker, eo, tan);
      jacobian_x_weights = new double[np];
      for(int i = 0; i < np; i++)
        jacobian_x_weights[i] = pt[i][2] * tan[i][2];
      order = eo;
      return np;
    }
  }
}