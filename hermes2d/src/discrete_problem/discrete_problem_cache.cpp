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

#include "global.h"
#include "mesh/traverse.h"
#include "space/space.h"
#include "shapeset/precalc.h"
#include "mesh/refmap.h"
#include "function/solution.h"
#include "discrete_problem/discrete_problem_cache.h"

namespace Hermes
{
  namespace Hermes2D
  {
    template<typename Scalar>
    void DiscreteProblemCache<Scalar>::CacheRecord::init(Traverse::State* current_state, PrecalcShapeset** current_pss, RefMap** current_refmaps, Solution<Scalar>** current_u_ext, AsmList<Scalar>** current_als, AsmList<Scalar>*** current_alsSurface, WeakForm<Scalar>* current_wf, int order)
    {
      this->spaceCnt = current_wf->get_neq();
      this->fns = new Func<double>**[spaceCnt];
      this->fnsSurface = NULL;
      this->asmlistCnt = new int[spaceCnt];
      this->asmlistIdx = new int*[spaceCnt];

      this->nvert = current_state->rep->nvert;
      this->order = order;

      int rep_space_i = -1;

      for(unsigned int space_i = 0; space_i < spaceCnt; space_i++)
      {
        if(current_state->e[space_i] == NULL)
        {
          this->fns[space_i] = NULL;
          continue;
        }
        else
          rep_space_i = space_i;

        this->asmlistCnt[space_i] = current_als[space_i]->cnt;
        this->asmlistIdx[space_i] = new int[current_als[space_i]->cnt];
        memcpy(this->asmlistIdx[space_i], current_als[space_i]->idx, current_als[space_i]->cnt * sizeof(int));

        this->fns[space_i] = new Func<double>*[this->asmlistCnt[space_i]];

        for (unsigned int j = 0; j < this->asmlistCnt[space_i]; j++)
        {
          current_pss[space_i]->set_active_shape(current_als[space_i]->idx[j]);
          this->fns[space_i][j] = init_fn(current_pss[space_i], current_refmaps[space_i], order);
        }
      }

      this->n_quadrature_points = init_geometry_points(current_refmaps[rep_space_i], this->order, this->geometry, this->jacobian_x_weights);

      if(current_state->isBnd && (current_wf->mfsurf.size() > 0 || current_wf->vfsurf.size() > 0))
      {
        this->geometrySurface = new Geom<double>*[this->nvert];
        this->jacobian_x_weightsSurface = new double*[this->nvert];

        this->n_quadrature_pointsSurface = new int[this->nvert];
        this->orderSurface = new int[this->nvert];

        int order = this->order;

        this->fnsSurface = new Func<double>***[this->nvert];
        memset(this->fnsSurface, NULL, sizeof(Func<double>***) * this->nvert);
        this->asmlistSurfaceCnt = new int*[this->nvert];

        for (current_state->isurf = 0; current_state->isurf < this->nvert; current_state->isurf++)
        {
          if(!current_state->bnd[current_state->isurf])
          {
            this->fnsSurface[current_state->isurf] = NULL;
            continue;
          }
          this->n_quadrature_pointsSurface[current_state->isurf] = init_surface_geometry_points(current_refmaps[rep_space_i], order, current_state->isurf, current_state->rep->marker, this->geometrySurface[current_state->isurf], this->jacobian_x_weightsSurface[current_state->isurf]);
          this->orderSurface[current_state->isurf] = order;
          order = this->order;

          this->fnsSurface[current_state->isurf] = new Func<double>**[spaceCnt];
          memset(this->fnsSurface[current_state->isurf], NULL, sizeof(Func<double>**) * spaceCnt);

          this->asmlistSurfaceCnt[current_state->isurf] = new int[spaceCnt];

          for(unsigned int space_i = 0; space_i < spaceCnt; space_i++)
          {
            if(current_state->e[space_i] == NULL)
              continue;

#ifdef _DEBUG
            assert(current_alsSurface[space_i][current_state->isurf]);
#endif

            this->asmlistSurfaceCnt[current_state->isurf][space_i] = current_alsSurface[space_i][current_state->isurf]->cnt;

            this->fnsSurface[current_state->isurf][space_i] = new Func<double>*[current_alsSurface[space_i][current_state->isurf]->cnt];
            for (unsigned int j = 0; j < current_alsSurface[space_i][current_state->isurf]->cnt; j++)
            {
              current_pss[space_i]->set_active_shape(current_alsSurface[space_i][current_state->isurf]->idx[j]);
              this->fnsSurface[current_state->isurf][space_i][j] = init_fn(current_pss[space_i], current_refmaps[space_i], this->orderSurface[current_state->isurf]);
            }
          }
        }
      }
    }

    template<typename Scalar>
    DiscreteProblemCache<Scalar>::CacheRecord::~CacheRecord()
    {
      this->free();
    }

    template<typename Scalar>
    void DiscreteProblemCache<Scalar>::init_assembling()
    {
      memset(this->hashTableUsed, 0, this->hash_table_size);
    }

    template<typename Scalar>
    void DiscreteProblemCache<Scalar>::CacheRecord::free()
    {
      for(unsigned int space_i = 0; space_i < spaceCnt; space_i++)
      {
        if(this->fns[space_i] == NULL)
          continue;
        for(unsigned int i = 0; i < this->asmlistCnt[space_i]; i++)
        {
          this->fns[space_i][i]->free_fn();
          delete this->fns[space_i][i];
        }
        delete [] this->fns[space_i];
        delete [] this->asmlistIdx[space_i];
      }
      delete [] this->fns;
      delete [] this->asmlistCnt;
      delete [] this->asmlistIdx;

      delete [] this->jacobian_x_weights;
      this->geometry->free();
      delete this->geometry;

      if(this->fnsSurface)
      {
        for(unsigned int edge_i = 0; edge_i < nvert; edge_i++)
        {
          if(this->fnsSurface[edge_i] == NULL)
            continue;

          this->geometrySurface[edge_i]->free();
          delete this->geometrySurface[edge_i];
          delete [] this->jacobian_x_weightsSurface[edge_i];

          for(unsigned int space_i = 0; space_i < spaceCnt; space_i++)
          {
            if(this->fnsSurface[edge_i][space_i] == NULL)
              continue;
            for(unsigned int i = 0; i < this->asmlistSurfaceCnt[edge_i][space_i]; i++)
            {
              this->fnsSurface[edge_i][space_i][i]->free_fn();
              delete this->fnsSurface[edge_i][space_i][i];
            }
            delete [] this->fnsSurface[edge_i][space_i];
          }
          delete [] this->fnsSurface[edge_i];
          delete [] this->asmlistSurfaceCnt[edge_i];
        }

        delete [] this->fnsSurface;
        delete [] this->geometrySurface;
        delete [] this->jacobian_x_weightsSurface;

        delete [] this->n_quadrature_pointsSurface;
        delete [] this->orderSurface;
        delete [] this->asmlistSurfaceCnt;

        this->fnsSurface = NULL;
      }
    }

    template<typename Scalar>
    DiscreteProblemCache<Scalar>::DiscreteProblemCache() : recordCount(0), size(DEFAULT_SIZE), hash_table_size(DEFAULT_HASH_TABLE_SIZE)
    {
      recordTable = (CacheRecord**)(malloc(size * sizeof(CacheRecord*)));
      hashTable = (StateHash**)(calloc(hash_table_size, sizeof(StateHash*)));
      hashTableUsed = (bool*)(calloc(hash_table_size, sizeof(bool)));
    }

    template<typename Scalar>
    void DiscreteProblemCache<Scalar>::free()
    {
      for(int i = 0; i < this->hash_table_size; i++)
        if(hashTable[i])
        {
          delete this->recordTable[this->hashTable[i]->cache_record_index];
          delete hashTable[i];
        }
        
      memset(this->recordTable, 0, this->size * sizeof(CacheRecord*));
      memset(this->hashTable, 0, this->hash_table_size * sizeof(StateHash*));
      memset(this->hashTableUsed, 0, this->hash_table_size * sizeof(bool));
    }

    template<typename Scalar>
    void DiscreteProblemCache<Scalar>::free_unused()
    {
      for(int i = 0; i < this->hash_table_size; i++)
      {
        if(!this->hashTableUsed[i] && this->hashTable[i])
        {
          delete this->recordTable[this->hashTable[i]->cache_record_index];
          delete this->hashTable[i];
          this->hashTable[i] = NULL;
        }
      }
    }

    template<typename Scalar>
    DiscreteProblemCache<Scalar>::~DiscreteProblemCache()
    {
      this->free();
      ::free(recordTable);
      ::free(hashTable);
      ::free(hashTableUsed);
    }

    template<typename Scalar>
    int DiscreteProblemCache<Scalar>::get_hash_record(int rep_id, int parent_son, int rep_sub_idx, int rep_i)
    {
      int hash = this->hashFunction(rep_id, parent_son, rep_sub_idx, rep_i);
      while (hashTable[hash] && (hashTable[hash]->rep_id != rep_id || hashTable[hash]->parent_son != parent_son || hashTable[hash]->rep_sub_idx != rep_sub_idx || hashTable[hash]->rep_i != rep_i))
        hash = (hash + 1) % hash_table_size;
      hashTableUsed[hash] = true;
      return hash;
    }

    template<typename Scalar>
    bool DiscreteProblemCache<Scalar>::get_adaptivity(Element* rep, int rep_sub_idx, int rep_i, CacheRecord*& cache_record)
    {
      int parent_son;
      for(int i = 0; i < 4; i++)
      {
        if(rep->parent->sons[i] == rep)
        {
          parent_son = i;
          break;
        }
      }

      int hash = this->get_hash_record(rep->parent->id, parent_son, rep_sub_idx, rep_i);
      if (this->hashTable[hash] == NULL)
      {
        if(recordCount > 0.9 * size)
        {
#pragma omp critical(record_size_increase)
          {
            size *= 1.5;
            this->recordTable = (CacheRecord**)realloc(this->recordTable, size * sizeof(CacheRecord*));
          }
        }
        cache_record = new CacheRecord();
#pragma omp critical(record_count_increase)
        {
          recordTable[recordCount] = cache_record;
          this->hashTable[hash] = new StateHash(rep->parent->id, parent_son, rep_sub_idx, rep_i, recordCount++);
        }
        return false;
      }
      else
      {
        cache_record = this->recordTable[this->hashTable[hash]->cache_record_index];
        return true;
      }
    }

    template<typename Scalar>
    bool DiscreteProblemCache<Scalar>::get(Element* rep, int rep_sub_idx, int rep_i, CacheRecord*& cache_record)
    {
      if(rep->parent)
        return get_adaptivity(rep, rep_sub_idx, rep_i, cache_record);
      else
      {
        int hash = this->get_hash_record(rep->id, -1, rep_sub_idx, rep_i);
        if (this->hashTable[hash] == NULL)
        {
          if(recordCount > 0.9 * size)
          {
#pragma omp critical(record_size_increase)
            {
              size *= 1.5;
              this->recordTable = (CacheRecord**)realloc(this->recordTable, size * sizeof(CacheRecord*));
            }
          }
          cache_record = new CacheRecord();
#pragma omp critical(record_count_increase)
          {
            recordTable[recordCount] = cache_record;
            this->hashTable[hash] = new StateHash(rep->id, -1, rep_sub_idx, rep_i, recordCount++);
          }
          return false;
        }
        else
        {
          cache_record = this->recordTable[this->hashTable[hash]->cache_record_index];
          return true;
        }
      }
    }

    template<typename Scalar>
    DiscreteProblemCache<Scalar>::StateHash::StateHash(int rep_id, int parent_son, int rep_sub_idx, int rep_i, int cache_record_index)
      : rep_id(rep_id), parent_son(parent_son), rep_sub_idx(rep_sub_idx), rep_i(rep_i), cache_record_index(cache_record_index) {};

    template<typename Scalar>
    int DiscreteProblemCache<Scalar>::hashFunction(int rep_id, int parent_son, int rep_sub_idx, int rep_i) const
    {
      if(rep_sub_idx < GUESS_NUMBER_OF_SUBELEMENTS)
        return rep_id * GUESS_NUMBER_OF_SUBELEMENTS + rep_sub_idx + rep_i + (parent_son + 1);
      else
        return (rep_id + 1) * GUESS_NUMBER_OF_SUBELEMENTS + rep_i + (parent_son + 1);
    }

    template class HERMES_API DiscreteProblemCache<double>;
    template class HERMES_API DiscreteProblemCache<std::complex<double> >;
  }
}