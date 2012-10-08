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

#ifndef __H2D_TRAVERSE_H
#define __H2D_TRAVERSE_H

#include "hermes_common.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views{
      class Orderizer;
      class Linearizer;
      class Vectorizer;
    };

    /// @defgroup inner Hermes hp-FEM/hp-DG assembling core

    /// @ingroup inner
    /// \brief Determines the position on an element surface (edge in 2D and Face in 3D).
    /// \details Used for the retrieval of boundary condition values.
    /// \details Same in H2D and H3D.
    ///
    struct SurfPos
    {
      int marker;    ///< surface marker (surface = edge in 2D and face in 3D)
      int surf_num;   ///< local element surface number

      Element *base; ///< for internal use

      int v1, v2;    ///< H2D only: edge endpoint vertex id numbers
      double t;      ///< H2D only: position between v1 and v2 in the range[0..1]
      double lo, hi; ///< H2D only: for internal use
    };

    class Mesh;
    class Transformable;
    struct State;
    struct Rect;

    /// @ingroup inner
    struct UniData
    {
      Element* e;
      uint64_t idx;
    };

    /// @ingroup inner
    static const uint64_t ONE = (uint64_t) 1 << 63;

    /// @ingroup inner
    struct Rect
    {
      uint64_t l, b, r, t;
    };

    /// @ingroup inner
    /// Traverse is a multi-mesh traversal utility class. Given N meshes sharing the
    /// same base mesh it walks through all (pseudo-)elements of the union of all
    /// the N meshes.
    ///
    class HERMES_API Traverse : public Hermes::Mixins::Loggable
    {
    public:
      Traverse(bool master = false);
    private:
      class State
      {
      public:
        Element** e;
        bool bnd[4];
        bool isBnd;
        Element* rep;
        ~State();
      private:
        State();
        void operator=(const State * other);
        void push_transform(int son, int i, bool is_triangle = false);
        uint64_t get_transform(int i);
        bool visited;
        uint64_t* sub_idx;
        Rect  cr;
        Rect* er;
        int num;
        int isurf;
      friend class Traverse;
      friend class Views::Linearizer;
      friend class Views::Vectorizer;
      template<typename Scalar> friend class DiscreteProblem;
      template<typename Scalar> friend class DiscreteProblemLinear;
      };

      void begin(int n, const Mesh** meshes, Transformable** fn = NULL);
      void finish();

      State* get_next_state(int* top_by_ref = NULL, int* id_by_ref = NULL);
      int get_num_states(Hermes::vector<const Mesh*> meshes);
      inline Element*  get_base() const { return base; }

      void init_transforms(State* s, int i);

      UniData** construct_union_mesh(Mesh* unimesh);

      int num;
      const Mesh** meshes;
      Transformable** fn;

      State* stack;
      int top, size;

      int id;
      bool tri;
      Element* base;
      int4* sons;
      uint64_t* subs;

      UniData** unidata;
      int udsize;

      State* push_state(int* top_by_ref = NULL);
      void set_boundary_info(State* s);
      void union_recurrent(Rect* cr, Element** e, Rect* er, uint64_t* idx, Element* uni);
      uint64_t init_idx(Rect* cr, Rect* er);

      void free_state(State* state);

      bool master;

      Mesh* unimesh;
      template<typename T> friend class Adapt;
      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemLinear;
      template<typename T> friend class Filter;
      template<typename T> friend class SimpleFilter;
      template<typename T> friend class Global;
      friend class Views::Orderizer;
      friend class Views::Vectorizer;
      friend class Views::Linearizer;
    };
  }
}
#endif