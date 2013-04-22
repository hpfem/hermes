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
#include "mesh.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      class Orderizer;
      class Linearizer;
      class Vectorizer;
    };

    /** @defgroup inner Hermes hp-FEM/hp-DG assembling core
    * Inner functionality classes that are not for the user to modify.
    */

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
      Traverse(int spaces_size);
      class State
      {
      public:
        Element** e;
        uint64_t* sub_idx;
        bool bnd[H2D_MAX_NUMBER_EDGES];
        bool isBnd;
        Element* rep;
        uint64_t rep_subidx;
        int rep_i;
        ~State();
        int isurf;
        int num;
      private:
        State();
        //void operator=(const State * other);
        static State* clone(const State * other);
        void push_transform(int son, int i, bool is_triangle = false);
        bool is_triangle();
        uint64_t get_transform(int i);
        bool visited;
        Rect  cr;
        Rect* er;
      friend class Traverse;
      friend class Views::Linearizer;
      friend class Views::Vectorizer;
      template<typename T> friend class DiscreteProblemCache;
      template<typename Scalar> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemAssemblyData;
      template<typename T> friend class DiscreteProblemDGAssembler;
      template<typename T> friend class DiscreteProblemThreadAssembler;
      };

      /// Returns all states on the passed meshes.
      /// \param[in] meshes Meshes.
      /// \param[out] num Number of states.
      /// \return The states.
      State** get_states(Hermes::vector<MeshSharedPtr> meshes, int& num);
      
    private:
      /// Used by get_states.
      void begin(int n);
      /// Used by get_states.
      void finish();
      /// Used by get_states.
      void init_transforms(State* s, int i);

#pragma region union-mesh
      UniData** construct_union_mesh(int n, MeshSharedPtr* meshes, MeshSharedPtr unimesh);
      void union_recurrent(Rect* cr, Element** e, Rect* er, uint64_t* idx, Element* uni);
      uint64_t init_idx(Rect* cr, Rect* er);

      UniData** unidata;
      int udsize;
      bool tri;
#pragma endregion

      /// Internal.
      int num;
      /// Internal.
      State* stack;
      /// Internal.
      int top, size;

      /// Internal.
      State* push_state(int* top_by_ref = NULL);
      /// Internal.
      void set_boundary_info(State* s);
      /// Internal.
      void free_state(State* state);
      /// Internal.
      int spaces_size;

      MeshSharedPtr unimesh;
      template<typename T> friend class Adapt;
      template<typename T> friend class KellyTypeAdapt;
      template<typename T> friend class DiscreteProblem;
      template<typename T> friend class DiscreteProblemCache;
      template<typename T> friend class DiscreteProblemDGAssembler;
      template<typename T> friend class DiscreteProblemIntegrationOrderCalculator;
      template<typename T> friend class DiscreteProblemAssemblyData;
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