/// This file is part of Hermes2D.
///
/// Hermes2D is free software: you can redistribute it and/or modify
/// it under the terms of the GNU General Public License as published by
/// the Free Software Foundation, either version 2 of the License, or
/// (at your option) any later version.
///
/// Hermes2D is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY;without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with Hermes2D. If not, see <http:///www.gnu.org/licenses/>.

#ifndef __H2D_DISCRETE_PROBLEM_CACHE_H
#define __H2D_DISCRETE_PROBLEM_CACHE_H

#include "hermes_common.h"
#include "adapt/adapt.h"
#include "graph.h"
#include "forms.h"
#include "weakform/weakform.h"
#include "function/function.h"
#include "neighbor.h"
#include "refinement_selectors/selector.h"
#include "exceptions.h"
#include "mixins2d.h"
#include "discrete_problem_helpers.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// @ingroup inner
    /// Caching in DiscreteProblem.
    /// Internal.
    template<typename Scalar>
    class DiscreteProblemCache
    {
    public:
      DiscreteProblemCache();

      /// Destructor that uses the clear() method and then deallocates even the internal structures.
      ~DiscreteProblemCache();

      /// Just clears all stored data, leaves the internal structures for further use.
      void free();

      /// In every call to assemble(), the unused hash table entries are stored, so that they can be deallocated (they will not be ever again used).
      void free_unused();

      /// Sets the usage table for each call to assemble().
      void init_assembling();

      /// Storage unit - a record.
      class CacheRecord
      {
      public:
        ~CacheRecord();
        void init(Traverse::State* state, PrecalcShapeset** current_pss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, AsmList<Scalar>*** current_alsSurface, WeakForm<Scalar>* current_wf, int order);
        void free();
        
        int spaceCnt;
        int** asmlistIdx;
        int* asmlistCnt;
        int nvert;
        int order;
        Func<double>*** fns;
        Func<double>**** fnsSurface;
        Geom<double>* geometry;
        Geom<double>** geometrySurface;
        double* jacobian_x_weights;
        double** jacobian_x_weightsSurface;
        int n_quadrature_points;
        int* n_quadrature_pointsSurface;
        int* orderSurface;
        int** asmlistSurfaceCnt;

        friend class DiscreteProblemCache<Scalar>;
      };

      /// Returns the cache record and information whether it is initialized (found in the cache).
      /// \param [out] cache_record The record.
      /// \return Found in cache.
      bool get(Element* rep, int rep_sub_idx, int rep_i, CacheRecord*& cache_record);

    private:
      /// Special handling of adaptivity situtation.
      bool get_adaptivity(Element* rep, int rep_sub_idx, int rep_i, CacheRecord*& cache_record);

      /// Starting size of the recordTable.
      static const int DEFAULT_SIZE = 1e5;
      /// Average number of subelements.
      static const int GUESS_NUMBER_OF_SUBELEMENTS = 16;
      /// Starting size of the hashTable.
      static const int DEFAULT_HASH_TABLE_SIZE = DEFAULT_SIZE * GUESS_NUMBER_OF_SUBELEMENTS;

      int size;
      int hash_table_size;

      CacheRecord **recordTable;
      int recordCount;

      class StateHash
      {
      private:
        /// Hash is created from 4 parameters, cache_record_index is an index to the array recordTable.
        /// \param[in] rep_id Id of the representing element of the Traverse::State at hand.
        /// \param[in] parent_son if dealing with adaptive calculation caching, the rep_id is no longer the id of the representing element, 
        /// but the id of its father and this is the son index that together represent the element at hand.
        /// \param[in] rep_sub_idx The sub-element number of the representing element.
        /// \param[in] rep_i In the case of subdomains calculations, this identifies what space is the representing element in.
        StateHash(int rep_id, int parent_son, int rep_sub_idx, int rep_i, int cache_record_index);

        int rep_id;
        int parent_son;
        int rep_sub_idx;
        int rep_i;
        int cache_record_index;
        friend class DiscreteProblemCache<Scalar>;
      };

      StateHash **hashTable;
      bool *hashTableUsed;

      int get_hash_record(int rep_id, int parent_son, int rep_sub_idx, int rep_i);

      int hashFunction(int rep_id, int parent_son, int rep_sub_idx, int rep_i) const;

      friend class DiscreteProblem<Scalar>;
    };
  }
}
#endif