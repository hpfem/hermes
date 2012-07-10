// This file is part of Hermes2D.
//
// Copyright 2009-2010 Ivo Hanak <hanak@byte.cz>
// Copyright 2005-2008 Jakub Cerveny <jakub.cerveny@gmail.com>
// Copyright 2005-2008 Lenka Dubcova <dubcova@gmail.com>
// Copyright 2005-2008 Pavel Solin <solin@unr.edu>
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

// $Id: view.h 1086 2008-10-21 09:05:44Z jakub $

#ifndef __H2D_SCALAR_VIEW_H
#define __H2D_SCALAR_VIEW_H

#include "view.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      // you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

      /// GUI inside views (the current GUI uses AntTweakBar and it is still experimental)
      /// to enable view define: ENABLE_VIEWER_GUI
# ifdef ENABLE_VIEWER_GUI
#   include <AntTweakBar.h>
#   define VIEWER_GUI(__def) __def
#   define VIEWER_GUI_CALLBACK(__clbk) if(__clbk) { refresh(); } else
# else
#   define TW_WND_ID_NONE -1
#   define VIEWER_GUI(__def)
#   define VIEWER_GUI_CALLBACK(__cblk)
#   define TwBar void /* avoid necessity to define ENABLE_VIEWER_GUI in underlaying applications */
# endif

      /// \brief Visualizes a Scalar PDE solution.
      ///
      /// ScalarView is a visualization window for all Scalar-valued PDE solutions.
      ///
      class HERMES_API ScalarView : public View
      {
      public:

        void init();
#ifndef _MSC_VER
        ScalarView(const char* title = "ScalarView", WinGeom* wg = NULL);
#endif
        ScalarView(char* title, WinGeom* wg = NULL);
        virtual ~ScalarView();

        void show(MeshFunction<double>* sln, double eps = HERMES_EPS_NORMAL, int item = H2D_FN_VAL_0,
          MeshFunction<double>* xdisp = NULL, MeshFunction<double>* ydisp = NULL, double dmult = 1.0);

        void show_linearizer_data(double eps = HERMES_EPS_NORMAL, int item = H2D_FN_VAL_0);

        inline void show_mesh(bool show = true) { show_edges = show; refresh(); }
        inline void show_bounding_box(bool show = true) { show_aabb = show; refresh(); }
        void show_contours(double step, double orig = 0.0);
        inline void hide_contours() { contours = false; refresh(); }
        inline void set_3d_mode(bool enable = true) { mode3d = enable; refresh(); }
        void set_vertical_scaling(double sc);  ///< Sets the scaling on the vertical axis programmatically.
        void set_min_max_range(double min, double max);  ///< Sets the limits on displayed values.

      public:
        Linearizer* lin;

      protected:
        /// Information about a vertex node.
        struct VertexNodeInfo
        {
          float x, y; ///< location of the node in coordinates of the mesh
          int id; ///< id of the node
          bool selected; ///< true if the node is selected
          TwBar* tw_bar; ///< a pointer to a gui window (CTwBar*) (GUI only).
          VertexNodeInfo() {}; ///< An empty default constructor to limit time
          VertexNodeInfo(int id, float x, float y) : x(x), y(y), id(id), selected(false), tw_bar(NULL) {};
        };
        Hermes::vector<VertexNodeInfo> vertex_nodes; ///< Vertex nodes. Sorted accordin to the X-axis.
        VertexNodeInfo* pointed_vertex_node; ///< A vertex node that is under the mouse cursor. NULL if none.

        bool allow_node_selection; ///> True if node selection is allowed
        unsigned int pointed_node_widget; ///> A GL display-list denoting a pointed vertex node. The geometry assumes the size of a pixel is 1x1.
        unsigned int selected_node_widget; ///> A GL display-list denoting a selected mesh node. The geometry assumes the size of a pixel is 1x1.

        const int node_pixel_radius; ///< A radius of node selection, in pixels.
        const int node_widget_vert_cnt; ///< A number of vertices for a mesh node widget.

        void init_vertex_nodes(const Mesh* mesh); ///< Creates a copy of vertex nodes for purpose of displaying and selection.
        VertexNodeInfo* find_nearest_node_in_range(float x, float y, float radius); ///< Finds nearest node in range.
        static bool compare_vertex_nodes_x(const VertexNodeInfo& a, const VertexNodeInfo& b); ///< Returns true, if a's X-axis coordinate is lower than b's one. Used to sort mesh nodes for searching purposes.
        void draw_vertex_nodes(); ///< Draws vertex nodes.
        void draw_single_vertex_node(const VertexNodeInfo& node); ///< Draws a single vertex node.
        void create_nodes_widgets(); ///< Creates vertex nodes widgets if not created already.

      protected:
        struct ElementInfo ///< element info structure
        {
          float x, y; ///< location of center[in physical coordinates]
          float width, height; ///< width, height of AABB[in physical coordinates]
          int id; ///< element ID
          ElementInfo() : x(0), y(0), width(0), height(0), id(-1) {};
          ElementInfo(int id, float x, float y, float width, float height) : x(x), y(y), width(width), height(height), id(id) {};
        };
        Hermes::vector<ElementInfo> element_infos; ///< Element info.

        unsigned int element_id_widget; ///< A GL display-list denoting a element ID widget. The geometry assumes the size of a pixel is 1x1.

        bool show_element_info; ///< true, to draw element info (currently ID) in 2D mode

        void init_element_info(const Mesh* mesh); ///< Creates element info from mesh.
        void create_element_info_widgets(); ///< Creates element ID widgets if not created already.
        void draw_element_infos_2d(); ///< Draws elements infos in 2D mode.

      protected: //GUI
        int tw_wnd_id; ///< tw window ID
        TwBar* tw_setup_bar; ///< setup bar

        virtual void create_setup_bar(); ///< create setup bar. Assumes that current window is set

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
          static const size_t H2D_OFFSETOF_COORD = 2*sizeof(float); ///< Offset of coordinate
        };
#pragma pack(pop)

        bool lin_updated; ///< true, if lin now contains new values

        unsigned int gl_coord_buffer; ///< Vertex coordinate buffer. (x, y, t)
        unsigned int gl_index_buffer; ///< Index data buffer.
        unsigned int gl_edge_inx_buffer; ///< A buffer for edge indices. The side of the buffer is H2DV_GL_MAX_EDGE_BUFFER pairs of indids.
        int max_gl_verts; ///< A maximum allocated number of vertices
        int max_gl_tris; ///< A maximum allocated number of triangles
        int gl_tri_cnt; ///< A number of OpenGL triangles

        bool show_values; ///< true to show values

        void prepare_gl_geometry(); ///< prepares geometry in a form compatible with GL arrays; Data are updated if lin is updated. In a case of a failure (out of memory), gl_verts is NULL and an old OpenGL rendering method has to be used.
        void draw_values_2d(); ///< draws values
        void draw_edges_2d(); ///< draws edges

        void draw_normals_3d(); ////< Draws normals of the 3d mesh. Used for debugging purposses only.

      protected: //edges
        bool show_edges; ///< true to show edges of mesh
        bool show_aabb;  ///< true to show the bounding box
        float edges_color[3]; ///< color of edges

        typedef void (*DrawSingleEdgeCallback)(int inx_vert_a, int inx_vert_b, ScalarView* viewer, void* param); ///< A callback function that draws edge using specified vertex indices. Param is user supplied parameter.

        void calculate_mesh_aabb(double* x_min, double* x_max, double* y_min, double* y_max); ///< Calculates AABB from edges.

        static void draw_gl_edge(int inx_vert_a, int inx_vert_b, ScalarView* viewer, void* param); ///< Draws edge specified by edge indices using GL. Functions assumes that data are locked.
        void draw_edges(DrawSingleEdgeCallback draw_single_edge, void* param, bool boundary_only); ///< Draws edges of elements and boundary of mesh. Functions assumes that data are locked.
        void draw_aabb(); ///< Draws the axes-aligned bounding box of the model. Assumes a model/view matrix to be the current matrix on the OpenGL stack.

      protected:
        bool contours; ///< true to enable drawing of contours
        double cont_orig, cont_step; ///< contour settings.
        float cont_color[3]; ///< color of contours (RGB)
        bool do_zoom_to_fit; ///< true to automatically translate the view so that the whole model si displayed
        bool is_constant; ///< true if the function to be displayed is constant

        // Perspective projection parameters.
        static const int fovy;        ///< Field of view in the vertical direction (in degrees).
        static const double znear;  ///< Distance of the near clipping plane of the viewing frustum from the camera.
        static const double zfar;     ///< Distance of the Far clipping plane of the viewing frustum from the camera.

        bool pmode, mode3d, panning;
        double xrot, yrot, xtrans, ytrans, ztrans;
        double xzscale, yscale, xctr, yctr, zctr;

        ///< Information about the range of vertex values.
        double value_irange, value_range_avg;

        double3* normals;

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
        double calculate_ztrans_to_fit_view(); ///< Calculates the z-coordinate (in eye coordinates) of the closest viewpoint from which we can still see the whole model. Assumes a model/view matrix to be the current matrix on the OpenGL stack.
        virtual void reset_view(bool force_reset); ///< Resets 2d and 3d view.
        virtual void update_layout(); ///< Updates layout, i.e., centers 2d and 3d mesh.

        void draw_tri_contours(double3* vert, int3* tri);
        void calculate_normals(double3* verts, int num_verts, int3* tris, int num_tris); ///< Initializes normals.
        void init_lighting();
        void update_mesh_info(); ///< Updates mesh info. Assumes that data lock is locked.

        virtual void on_display();
        virtual void on_key_down(unsigned char key, int x, int y);
        virtual void on_special_key(int key, int x, int y);
        virtual void on_mouse_move(int x, int y);
        virtual void on_right_mouse_down(int x, int y); ///< Handles selecting/deselecting of nodes.
        virtual void on_middle_mouse_down(int x, int y);
        virtual void on_middle_mouse_up(int x, int y);
        virtual const char* get_help_text() const;
        virtual void on_close();
        virtual void on_create(int output_id);

        virtual void on_left_mouse_down(int x, int y);
        virtual void on_left_mouse_up(int x, int y);
        virtual void on_right_mouse_up(int x, int y);
        virtual void on_reshape(int width, int height);
      };
#else
class HERMES_API ScalarView : public View
      {
      public:
        void init() { throw Hermes::Exceptions::Exception("GLUT disabled."); }
#ifndef _MSC_VER
        ScalarView(const char* title = "ScalarView", WinGeom* wg = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
#endif
        ScalarView(char* title, WinGeom* wg = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        void show(MeshFunction<double>* sln, double eps = HERMES_EPS_NORMAL, int item = H2D_FN_VAL_0,
          MeshFunction<double>* xdisp = NULL, MeshFunction<double>* ydisp = NULL, double dmult = 1.0) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        void show_linearizer_data(double eps = HERMES_EPS_NORMAL, int item = H2D_FN_VAL_0) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        inline void show_mesh(bool show = true) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        inline void show_bounding_box(bool show = true) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void show_contours(double step, double orig = 0.0) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        inline void hide_contours() { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        inline void set_3d_mode(bool enable = true) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void set_vertical_scaling(double sc) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void set_min_max_range(double min, double max) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
      };
#endif
    }
  }
}
#endif