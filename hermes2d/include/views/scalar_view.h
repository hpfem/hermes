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
/*! \file scalar_view.h
\brief File containing ScalarView class.
*/
#ifndef __H2D_SCALAR_VIEW_H
#define __H2D_SCALAR_VIEW_H

#include "view.h"
#include "linearizer.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      // you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

      /// \brief Visualizes a Scalar PDE solution.
      ///
      /// ScalarView is a visualization window for all Scalar-valued PDE solutions.
      ///
      class HERMES_API ScalarView : public View
      {
      public:

        void init();
#ifndef _MSC_VER
        ScalarView(const char* title = "ScalarView", WinGeom* wg = nullptr);
#else
        ScalarView(WinGeom* wg = nullptr);
#endif
        ScalarView(char* title, WinGeom* wg = nullptr);
        ~ScalarView();

        void show(MeshFunctionSharedPtr<double> sln, int item = H2D_FN_VAL_0,
          MeshFunctionSharedPtr<double> xdisp = nullptr, MeshFunctionSharedPtr<double> ydisp = nullptr, double dmult = 1.0);

        void show(MeshFunctionSharedPtr<std::complex<double> > sln, int item = H2D_FN_VAL_0,
          MeshFunctionSharedPtr<double> xdisp = nullptr, MeshFunctionSharedPtr<double> ydisp = nullptr, double dmult = 1.0)
        {
          throw Exceptions::Exception("Visualization of complex 2D solution is not possible, please use a filter that converts the solution into a real function, then display that one.");
        }

        void show_linearizer_data(double eps, int item = H2D_FN_VAL_0);

        inline void show_mesh(bool show = true) { show_edges = show; refresh(); }
        inline void show_bounding_box(bool show = true) { show_aabb = show; refresh(); }
        void show_contours(double step, double orig = 0.0);
        inline void hide_contours() { contours = false; refresh(); }
        void set_3d_mode(bool enable = true);
        /// Sets the scaling on the vertical axis programmatically.
        void set_vertical_scaling(double sc);
        /// Sets the limits on displayed values.
        void set_min_max_range(double min, double max);

        /// Resets 2d and 3d view.
        virtual void reset_view(bool force_reset);

        /// Returns the internal linearizer for the purpose of parameter settings.
        Linearizer* get_linearizer();

        /// Sets the criterion to use for the linearization process.
        /// This criterion is used in ThreadLinearizerMultidimensional class instances (see threadLinearizerMultidimensional array).
        /// \param[in] criterion The instance of the criterion - see the class LinearizerCriterion for details (method split_decision() for the adaptive criterion, process_[triangle|quad] for the fixed one).
        void set_linearizer_criterion(LinearizerCriterion criterion);

      protected:
        /// LinearizerMultidimensional class responsible for obtaining linearized data.
        Linearizer* lin;

      protected:
        struct ElementInfo ///< element info structure
        {
          /// location of center[in physical coordinates]
          float x, y;
          /// width, height of AABB[in physical coordinates]
          float width, height;
          /// element ID
          int id;
          ElementInfo() : x(0), y(0), width(0), height(0), id(-1) {};
          ElementInfo(int id, float x, float y, float width, float height) : x(x), y(y), width(width), height(height), id(id) {};
        };
        /// Element info.
        std::vector<ElementInfo> element_infos;

        /// A GL display-list denoting a element ID widget. The geometry assumes the size of a pixel is 1x1.
        unsigned int element_id_widget;

        /// true, to draw element info (currently ID) in 2D mode
        bool show_element_info;

        /// Creates element info from mesh.
        void init_element_info(MeshSharedPtr mesh);
        /// Creates element ID widgets if not created already.
        void create_element_info_widgets();
        /// Draws elements infos in 2D mode.
        void draw_element_infos_2d();

      protected: //values
#define H2DV_GL_MAX_EDGE_BUFFER 128 ///< A maximum number of pairs per a buffer.
#pragma pack(push)
#pragma pack(1)
        struct GLVertex2 ///< OpenGL vertex. Used to cache vertices prior rendering
        {
          float x, y;
          float coord;
          GLVertex2() {};
          GLVertex2(float x, float y, float coord) : x(x), y(y), coord(coord) {};
          /// Offset of coordinate
          static const size_t H2D_OFFSETOF_COORD = 2 * sizeof(float);
        };
#pragma pack(pop)

        /// true, if lin now contains new_ values
        bool lin_updated;

        /// A maximum allocated number of vertices
        int max_gl_verts;
        /// A maximum allocated number of triangles
        int max_gl_tris;
        /// A number of OpenGL triangles
        int gl_tri_cnt;

        /// true to show values
        bool show_values;

        /// prepares geometry in a form compatible with GL arrays; Data are updated if lin is updated. In a case of a failure (out of memory), gl_verts is nullptr and an old OpenGL rendering method has to be used.
        void prepare_gl_geometry();
        /// draws values
        void draw_values_2d();
        /// draws edges
        void draw_edges_2d();

      protected: //edges
        /// true to show edges of mesh
        bool show_edges;
        /// true to show the bounding box
        bool show_aabb;
        /// color of edges
        float edges_color[3];

        /// A callback function that draws edge using specified vertex indices. Param is user supplied parameter.
        typedef void(*DrawSingleEdgeCallback)(int inx_vert_a, int inx_vert_b, ScalarView* viewer, void* param);

        /// Calculates AABB from edges.
        void calculate_mesh_aabb(double* x_min, double* x_max, double* y_min, double* y_max);

        /// Draws the axes-aligned bounding box of the model. Assumes a model/view matrix to be the current matrix on the OpenGL stack.
        void draw_aabb();

      protected:
        /// true to enable drawing of contours
        bool contours;
        /// contour settings.
        double cont_orig, cont_step;
        /// color of contours (RGB)
        float cont_color[3];
        /// true to automatically translate the view so that the whole model si displayed
        bool do_zoom_to_fit;
        /// true if the function to be displayed is constant
        bool is_constant;

        // Perspective projection parameters.
        /// Field of view in the vertical direction (in degrees).
        static const int fovy;
        /// Distance of the near clipping plane of the viewing frustum from the camera.
        static const double znear;
        /// Distance of the Far clipping plane of the viewing frustum from the camera.
        static const double zfar;

        bool pmode, mode3d, panning;
        double xrot, yrot, xtrans, ytrans, ztrans;
        double xzscale, yscale, xctr, yctr, zctr;

        ///< Information about the range of vertex values.
        double value_irange, value_range_avg;

        /// This function calculates the distance that the model (3D plot of the solution over the whole solution domain) must be
        /// translated along the z-axis of the eye coordinate system, so that it fills the actual viewport without being clipped.
        /// The only case when the model will be clipped is when the user defines his own vertical range limits - unfortunately,
        /// the values beyond these limits are now still included in the displayed model, so the user may e.g. zoom-out to see them
        /// - try the benchmark 'screen' for instance...
        /// The algorithm essentially performs zooming without changing the field of view of the perspective projection and works as follows:
        ///  1. Apply all transformations to the bounding box of the model at the origin.
        ///  2. Compute the distance (along z-axis) from the origin to the center of perspective projection of the point with the
        ///      biggest vertical (y-axis) distance from the origin.
        ///  3. Compute the distance (along z-axis) from the origin to the center of perspective projection of the point with the
        ///      biggest horizontal (x-axis) distance from the origin.
        ///  4. Take the bigger of the two distances and reverse sign (since we will translate the model, not the camera)
        /// Calculates the z-coordinate (in eye coordinates) of the closest viewpoint from which we can still see the whole model. Assumes a model/view matrix to be the current matrix on the OpenGL stack.
        double calculate_ztrans_to_fit_view();
        /// Updates layout, i.e., centers 2d and 3d mesh.
        virtual void update_layout();

        void draw_tri_contours(ScalarLinearizerDataDimensions<LINEARIZER_DATA_TYPE>::triangle_t&);
        void init_lighting();
        /// Updates mesh info. Assumes that data lock is locked.
        void update_mesh_info();

        virtual void on_display();
        void on_display_2d();
        void on_display_3d();

        virtual void on_key_down(unsigned char key, int x, int y);
        virtual void on_mouse_move(int x, int y);
        /// Handles selecting/deselecting of nodes.
        virtual void on_right_mouse_down(int x, int y);
        virtual void on_middle_mouse_down(int x, int y);
        virtual void on_middle_mouse_up(int x, int y);
        virtual const char* get_help_text() const;
        virtual void on_close();
      };
#else
      class HERMES_API ScalarView : public View
      {
      public:
        void init() { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        ScalarView(const char* title = "ScalarView", WinGeom* wg = nullptr) {}
        ScalarView(char* title, WinGeom* wg = nullptr) {}

        void show(MeshFunctionSharedPtr<double> sln, int item = H2D_FN_VAL_0,
          MeshFunctionSharedPtr<double> xdisp = nullptr, MeshFunctionSharedPtr<double> ydisp = nullptr, double dmult = 1.0) {
          throw Hermes::Exceptions::Exception("GLUT disabled.");
        }

        void show(MeshFunctionSharedPtr<std::complex<double> > sln, int item = H2D_FN_VAL_0,
          MeshFunctionSharedPtr<double> xdisp = nullptr, MeshFunctionSharedPtr<double> ydisp = nullptr, double dmult = 1.0)
        {
          throw Exceptions::Exception("Visualization of complex 2D solution is not possible, please use a filter that converts the solution into a real function, then display that one.");
        }

        void show_linearizer_data(double eps, int item = H2D_FN_VAL_0) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        inline void show_mesh(bool show = true) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        inline void show_bounding_box(bool show = true) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void show_contours(double step, double orig = 0.0) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        inline void hide_contours() { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        inline void set_3d_mode(bool enable = true) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void set_vertical_scaling(double sc) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void set_min_max_range(double min, double max) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        Linearizer* get_linearizer()  { throw Hermes::Exceptions::Exception("GLUT disabled."); return nullptr; }
        void set_linearizer_criterion(LinearizerCriterion criterion) { throw Hermes::Exceptions::Exception("GLUT disabled."); return nullptr; }
      };
#endif
    }
  }
}
#endif