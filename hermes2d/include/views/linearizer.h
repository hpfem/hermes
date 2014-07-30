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

#ifndef __H2D_LINEARIZER_H
#define __H2D_LINEARIZER_H

#include "thread_linearizer.h"
#include "linearizer_utils.h"
#ifndef NOGLUT
  #include <pthread.h>
#endif
#include "../function/solution.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      /// LinearizerMultidimensional is a utility class which converts a higher-order FEM solution defined on
      /// a curvilinear, irregular mesh to a linear FEM solution defined on a straight-edged,
      /// regular mesh. This is done by adaptive refinement of the higher-order mesh and its
      /// subsequent regularization. The linearized mesh can then be easily displayed or
      /// exported to standard formats. The class correctly handles discontinuities in the
      /// solution (e.g., gradients or in Hcurl) by inserting double vertices where necessary.
      /// LinearizerMultidimensional also serves as a container for the resulting linearized mesh.
      template<typename LinearizerDataDimensions>
      class HERMES_API LinearizerMultidimensional :
        public Hermes::Mixins::TimeMeasurable,
        public Hermes::Mixins::Loggable,
        public Hermes::Hermes2D::Mixins::Parallel
      {
      public:
        LinearizerMultidimensional(LinearizerOutputType linearizerOutputType, bool auto_max = true);
        ~LinearizerMultidimensional();

        /// Main method - processes the solution and stores the data obtained by the process.
        /// \param[in] sln the solution
        /// \param[in] item what item (function value, derivative wrt. x, ..) to use in the solution.
        void process_solution(MeshFunctionSharedPtr<double> sln, int item = H2D_FN_VAL_0);
        void process_solution(MeshFunctionSharedPtr<double>* sln, int* items);

        /// Save a MeshFunction (Solution, Filter) in VTK format.
        void save_solution_vtk(MeshFunctionSharedPtr<double> sln, const char* filename, const char* quantity_name, bool mode_3D = true, int item = H2D_FN_VAL_0);
        void save_solution_vtk(Hermes::vector<MeshFunctionSharedPtr<double> > slns, Hermes::vector<int> items, const char* filename, const char* quantity_name, bool mode_3D = true);
        void save_solution_tecplot(MeshFunctionSharedPtr<double> sln, const char* filename, const char* quantity_name, int item = H2D_FN_VAL_0);
        void save_solution_tecplot(Hermes::vector<MeshFunctionSharedPtr<double> > slns, Hermes::vector<int> items, const char* filename, Hermes::vector<std::string> quantity_names);

        void set_criterion(LinearizerCriterion criterion);

        /// Set the displacement, i.e. set two functions that will deform the domain for visualization, in the x-direction, and the y-direction.
        void set_displacement(MeshFunctionSharedPtr<double> xdisp, MeshFunctionSharedPtr<double> ydisp, double dmult = 1.0);

        /// Iterator class.
        template<typename T>
        class Iterator
        {
        public:
          Iterator(const LinearizerMultidimensional<LinearizerDataDimensions>* linearizer);
          void operator++();
          T& get() const;
          /// For triangle- and edge- markers.
          int& get_marker() const;
          bool end;
        private:
          int current_thread_index;
          int current_thread;
          int current_thread_size;
          Hermes::vector<int> thread_sizes;
          const LinearizerMultidimensional<LinearizerDataDimensions>* linearizer;
          void check_zero_lengths();
          friend class LinearizerMultidimensional;
        };

        /// Begin - iterators.
        Iterator<typename LinearizerDataDimensions::vertex_t> vertices_begin() const;
        Iterator<typename LinearizerDataDimensions::triangle_t> triangles_begin() const;
        Iterator<typename LinearizerDataDimensions::edge_t> edges_begin() const;
        Iterator<triangle_indices_t> triangle_indices_begin() const;

        /// Counts
        int get_vertex_count() const;
        int get_triangle_count() const;
        int get_edge_count() const;
        int get_triangle_index_count() const;

        void set_max_absolute_value(double max_abs);

        double get_min_value() const;
        double get_max_value() const;

        /// Set the threshold for how fine the output for curved elements.
        /// \param[in] curvature_epsilon The 'curvature' epsilon determining the tolerance of catching the shape of curved elements.
        /// The smaller, the finer.
        /// Default value = 1e-3.
        void set_curvature_epsilon(double curvature_epsilon);

        /// Get the 'curvature' epsilon determining the tolerance of catching the shape of curved elements.
        double get_curvature_epsilon() const;

        /// Free the instance.
        void free();

        /// Assigned criterion.
        /// Defaults to a fixed criterion with one level of refinement.
        LinearizerCriterion criterion;

        /// Internal.
        void lock_data() const;
        void unlock_data() const;

        /// Returns axis aligned bounding box (AABB) of vertices. Assumes lock.
        void calc_vertices_aabb(double* min_x, double* max_x, double* min_y, double* max_y) const;

      protected:
        Quad2D *old_quad[LinearizerDataDimensions::dimension], *old_quad_x, *old_quad_y;

        LinearizerOutputType linearizerOutputType;

        /// Before assembling.
        void check_data(MeshFunctionSharedPtr<double>* sln);

        /// Assembly data.
        ThreadLinearizerMultidimensional<LinearizerDataDimensions>** threadLinearizerMultidimensional;

        void init(MeshFunctionSharedPtr<double>* sln, int* item);

        Hermes::vector<MeshSharedPtr > meshes;

        /// Standard and curvature epsilon.
        double curvature_epsilon;

        /// Information if user-supplied displacement functions have been provided.
        bool user_xdisp, user_ydisp;

#ifndef NOGLUT
        mutable pthread_mutex_t data_mutex;
#endif

        /// Displacement functions, default to ZeroFunctions, may be supplied by set_displacement();
        MeshFunctionSharedPtr<double> xdisp, ydisp;

        /// Displacement multiplicator - used e.g. in Elasticity to multiply the displacement to make it more noticeable.
        double dmult;

        /// What kind of information do we want to get out of the solution.
        int item[LinearizerDataDimensions::dimension], component[LinearizerDataDimensions::dimension], value_type[LinearizerDataDimensions::dimension];

        // Finish - contour triangles calculation etc.
        void finish(MeshFunctionSharedPtr<double>* sln);

        Traverse::State** states;

        int num_states;

        double min_val, max_val;

        void find_min_max();

        friend class ThreadLinearizerMultidimensional<LinearizerDataDimensions>;
      };

      typedef LinearizerMultidimensional<ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE> > Linearizer;
      typedef LinearizerMultidimensional<VectorLinearizerDataDimensions<LINEARIZER_DATA_TYPE> > Vectorizer;
    }
  }
}
#endif
