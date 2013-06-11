// This file is part of Hermes2D.
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

#ifndef __H2D_VIEW_H
#define __H2D_VIEW_H

#include "../global.h"
#include "vectorizer.h"
#include "orderizer.h"

namespace Hermes
{
  namespace Hermes2D
  {
    /// Namespace containing Views classes, Linearizer classes, and support for these.
    namespace Views
    {
      // Constants
#define H2D_DEFAULT_X_POS -1
#define H2D_DEFAULT_Y_POS -1
#define H2D_DEFAULT_WIDTH 600
#define H2D_DEFAULT_HEIGHT 400
#define H2DV_SCALE_LOG_BASE 1.005 ///< Base of the scale coefficient. Scale = base^{mouse move}.

      /// Wait events.
      enum ViewWaitEvent {
        HERMES_WAIT_CLOSE, ///< Wait for all windows to close.
        HERMES_WAIT_KEYPRESS, ///< Wait for any unprocessed keypress to happen.
      };

      struct WinGeom {
        int x;
        int y;
        int width;
        int height;

        WinGeom(int x, int y, int width, int height) {
          this->x = x;
          this->y = y;
          this->width = width;
          this->height = height;
        };
      };

      /// View palette type.
      enum ViewPaletteType {
        H2DV_PT_DEFAULT = -1, ///< Default palette. Depends on viewer.
        H2DV_PT_HUESCALE = 0, ///< A palette based on hue scale.
        H2DV_PT_GRAYSCALE = 1, ///< Greyscale.
        H2DV_PT_INVGRAYSCALE = 2, ///< Inverted grayscale.
        H2DV_PT_MAX_ID = 3 ///< Maximum ID of view palette type.
      };

      // you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

      /// \brief Represents a simple visualization window.
      ///
      /// View is a base class providing a simple OpenGL visualization window.
      /// Its task is to define basic functionality, such as the ability of the
      /// window to be responsive even when the main program thread is busy
      /// with calculations (ie., the windows are run in a background thread),
      /// to provide zooming and panning capabilities for use by the descendant
      /// classes, etc.
      ///
      class HERMES_API View : public Hermes::Mixins::TimeMeasurable, public Hermes::Mixins::Loggable
      {
      public:

        void init();
        View(const char* title, WinGeom* wg = NULL);
        View(char* title, WinGeom* wg = NULL);
        ~View();

        int  create();
        void close();
        void refresh(); ///< Refreshes views

        /// Returns the title.
        const char* get_title() const;
        /// Changes the window name (in its title-bar) to 'title'.
        void set_title(const char* title);

        void set_min_max_range(double min, double max);
        void auto_min_max_range();
        void get_min_max_range(double& min, double& max);

        void show_scale(bool show = true);
        void set_scale_position(int horz, int vert);
        void set_scale_size(int width, int height, int numticks);
        void set_scale_format(const char* fmt);
        void fix_scale_width(int width = 80);

        /// Saves the current content of the window to a .BMP file.
        /// If 'high_quality' is true, an anti-aliased frame is rendered and saved.
        void save_screenshot(const char* bmpname, bool high_quality = false);
        /// Like save_screenshot(), but forms the file name in printf-style using the 'number'
        /// parameter, e.g., format = "screen%03d.bmp" and number = 5 gives the file name "screen005.bmp".
        void save_numbered_screenshot(const char* format, int number, bool high_quality = false);

        void set_palette(ViewPaletteType type);
        void set_num_palette_steps(int num);
        void set_palette_filter(bool linear);

        void wait_for_keypress(const char* text = NULL); ///< Waits for keypress. Deprecated.
        void wait_for_close();
        void wait_for_draw();

        static void wait(const char* text); ///< Closes all views at once.
        static void wait(ViewWaitEvent wait_event = HERMES_WAIT_CLOSE, const char* text = NULL); ///< Waits for an event.
        void draw_help();
        virtual void reset_view(bool force_reset); ///< Resets view based on the axis-aligned bounding box of the mesh. Assumes that the bounding box is set up. Does not reset if view_not_reset is false.

      protected: //FPS measurement
#define FPS_FRAME_SIZE 5
        double rendering_frames[FPS_FRAME_SIZE]; ///< time spend in rendering of frames[in ms]
        int rendering_frames_top; ///< the new location of the next FPS
        void draw_fps(); ///< draws current FPS
        static double get_tick_count(); ///< returns a current time[in ms]

      protected: //view
        bool view_not_reset; ///< True if the view was not reset and therefore it has to be.
        double vertices_min_x, vertices_max_x, vertices_min_y, vertices_max_y; ///< AABB of shown mesh
        double scale, log_scale, trans_x, trans_y;
        double center_x, center_y;
        int margin, lspace, rspace;
        int mouse_x, mouse_y;
        int scx, scy;
        double objx, objy;
        bool dragging, scaling;

        virtual void on_create(int output_id);
        virtual void on_display() {};
        virtual void on_reshape(int width, int height);
        virtual void on_mouse_move(int x, int y);
        virtual void on_left_mouse_down(int x, int y);
        virtual void on_left_mouse_up(int x, int y);
        virtual void on_left_mouse_double_click(int x, int y) {}
        virtual void on_right_mouse_down(int x, int y);
        virtual void on_right_mouse_up(int x, int y);
        virtual void on_right_mouse_double_click(int x, int y) {}
        virtual void on_middle_mouse_down(int x, int y) {}
        virtual void on_middle_mouse_up(int x, int y) {}
        virtual void on_middle_mouse_double_click(int x, int y) {}
        virtual void on_key_down(unsigned char key, int x, int y);
        virtual void on_special_key(int key, int x, int y);
        virtual void on_entry(int state) {}
        virtual void on_close();

        virtual void update_layout(); ///< Updates layout, i.e., centers mesh.

      protected:
        std::string title;
        int output_id;
        int output_x, output_y, output_width, output_height;
        float jitter_x, jitter_y;
        bool hq_frame, frame_ready;

        ViewPaletteType pal_type;
        int pal_steps, pal_filter;
        double tex_scale, tex_shift;
        bool range_auto;
        double range_min, range_max;

        bool b_scale, b_help;
        bool scale_focused, scale_dragging;
        int pos_horz, pos_vert;
        int scale_x, scale_y;
        int scale_width, scale_height, labels_width;
        int scale_numticks, scale_box_height, scale_box_skip;
        char scale_fmt[20];
        int scale_fixed_width;

        bool want_screenshot;
        static int screenshot_no;
        std::string screenshot_filename;

      protected: //palette
        unsigned int gl_pallete_tex_id; ///< OpenGL texture object ID

        void create_gl_palette(); ///< Creates pallete texture in OpenGL. Assumes that view_sync is locked.
        virtual void get_palette_color(double x, float* gl_color); ///< Fills gl_color with palette color. Assumes that gl_color points to a vector of three components (RGB).

      protected: //internal functions
        inline double transform_x(double x) { return (x * scale + trans_x) + center_x; }
        inline double transform_y(double y) { return center_y - (y * scale + trans_y); }
        inline double untransform_x(double x) { return (x - center_x - trans_x) / scale; }
        inline double untransform_y(double y) { return (center_y - y - trans_y) / scale; }

        virtual void clear_background(); ///< Clears background.
        void pre_display();
        void display_antialiased();

        void set_ortho_projection(bool no_jitter = false);
        void set_3d_projection(int fov, double znear, double zfar);

        void draw_text(double x, double y, const char* text, int align = -1);
        int  get_text_width(const char* text);

        char *get_screenshot_file_name();
        void save_screenshot_internal(const char* filename);

        virtual void scale_dispatch();
        virtual int measure_scale_labels();
        void draw_continuous_scale(char* title, bool righttext);
        void draw_discrete_scale(int numboxes, const char* boxnames[], const float boxcolors[][3]);

        void update_tex_adjust();

        virtual const char* get_help_text() const = 0;

        friend void on_display_stub(void);
        friend void on_reshape_stub(int, int);
        friend void on_mouse_move_stub(int, int);
        friend void on_mouse_click_stub(int, int, int, int);
        friend void on_key_down_stub(unsigned char, int, int);
        friend void on_special_key_stub(int, int, int);
        friend void on_entry_stub(int);
        friend void on_idle_stub();
        friend void on_close_stub();
        friend int add_view_in_thread(void*);
        friend int remove_view_in_thread(void*);
        friend void on_create(int);
      };
#else
class HERMES_API View
      {
      public:

        void init() { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        View() { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        View(const char* title, WinGeom* wg = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        View(char* title, WinGeom* wg = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        int  create() { throw Hermes::Exceptions::Exception("GLUT disabled."); return -1; }
        void close() { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void refresh() { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        /// Changes the window name (in its title-bar) to 'title'.
        void set_title(const char* title) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        void set_min_max_range(double min, double max) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void auto_min_max_range() { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void get_min_max_range(double& min, double& max) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        void show_scale(bool show = true) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void set_scale_position(int horz, int vert) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void set_scale_size(int width, int height, int numticks) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void set_scale_format(const char* fmt) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void fix_scale_width(int width = 80) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        /// Saves the current content of the window to a .BMP file.
        /// If 'high_quality' is true, an anti-aliased frame is rendered and saved.
        void save_screenshot(const char* bmpname, bool high_quality = false) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        /// Like save_screenshot(), but forms the file name in printf-style using the 'number'
        /// parameter, e.g., format = "screen%03d.bmp" and number = 5 gives the file name "screen005.bmp".
        void save_numbered_screenshot(const char* format, int number, bool high_quality = false) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        void set_palette(ViewPaletteType type) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void set_num_palette_steps(int num) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void set_palette_filter(bool linear) { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        void wait_for_keypress(const char* text = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void wait_for_close() { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        void wait_for_draw() { throw Hermes::Exceptions::Exception("GLUT disabled."); }

        static void wait(const char* text) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
        static void wait(ViewWaitEvent wait_event = HERMES_WAIT_CLOSE, const char* text = NULL) { throw Hermes::Exceptions::Exception("GLUT disabled."); }
      };
#endif
    }
  }
}
#endif