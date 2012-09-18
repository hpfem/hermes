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

#ifndef NOGLUT

#include <sstream>
#include <GL/glew.h>
#include <GL/freeglut.h>
#ifndef WIN32
# include <sys/time.h>
#endif

#include "global.h"
#include "view_support.h"
#include "view.h"
#include "solution.h"
#include "view_data.cpp"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
#define HERMES_WAIT_CLOSE_MSG "close all views to continue"
#define HERMES_WAIT_KEYPRESS_MSG "press spacebar to continue"

      int View::screenshot_no = 1;

      void View::init()
      {
        jitter_x = jitter_y = 0.0;
        dragging = scaling = false;
        hq_frame = false;
        frame_ready = false;
        range_auto = true;
        range_min = 0;
        range_max = 1;
        pal_type = H2DV_PT_HUESCALE;
        pal_steps = 50;
        pal_filter = GL_NEAREST;
        margin = 15;
        b_scale = true;
        b_help = false;
        scale_focused = scale_dragging = false;
        pos_horz = pos_vert = 0;
        scale_x = scale_y = labels_width = 0;
        scale_width = 16;
        scale_height = 320;
        scale_numticks = 9;
        strcpy(scale_fmt, "%.3g");
        scale_fixed_width = -1;
        want_screenshot = false;

        rendering_frames_top = 0;
        memset(rendering_frames, 0, FPS_FRAME_SIZE * sizeof(double));
      }

      View::View(const char* title, WinGeom* wg) :
      view_not_reset(true),
        vertices_min_x(0),
        vertices_max_x(0),
        vertices_min_y(0),
        vertices_max_y(0),
        title(title),
        output_id(-1),
        gl_pallete_tex_id(0)
      {
        if(wg == NULL)
        {
          output_x = H2D_DEFAULT_X_POS;
          output_y = H2D_DEFAULT_Y_POS;
          output_width = H2D_DEFAULT_WIDTH;
          output_height = H2D_DEFAULT_HEIGHT;
        }
        else
        {
          output_x = wg->x;
          output_y = wg->y;
          output_width = wg->width;
          output_height = wg->height;
        }

        init();
      }

      View::View(char* title, WinGeom* wg) :
      view_not_reset(true),
        vertices_min_x(0),
        vertices_max_x(0),
        vertices_min_y(0),
        vertices_max_y(0),
        title(title),
        output_id(-1),
        gl_pallete_tex_id(0)
      {
        if(wg == NULL)
        {
          output_x = H2D_DEFAULT_X_POS;
          output_y = H2D_DEFAULT_Y_POS;
          output_width = H2D_DEFAULT_WIDTH;
          output_height = H2D_DEFAULT_HEIGHT;
        }
        else
        {
          output_x = wg->x;
          output_y = wg->y;
          output_width = wg->width;
          output_height = wg->height;
        }

        init();
      }

      View::~View()
      {
        if(output_id >= 0)
          close();
      }

      int View::create()
      {
        if(output_id < 0) //does not need thread protection because it is set up by a callback during call of add_view
          return add_view(this, output_x, output_y, output_width, output_height, title.c_str());
        else
          return output_id;
      }

      void View::close()
      {
        if(output_id >= 0) //does not need thread protection because it is set up by a callback during call of add_view
          remove_view(output_id);
      }

      void View::wait(const char* text)
      {
        wait(HERMES_WAIT_CLOSE, text);
      }

      void View::wait(ViewWaitEvent wait_event, const char* text)
      {
        //prepare message
        std::stringstream str;
        str << "  << ";
        if(text != NULL)
          str << text;
        else
        {
          switch(wait_event)
          {
          case HERMES_WAIT_CLOSE: str << HERMES_WAIT_CLOSE_MSG; break;
          case HERMES_WAIT_KEYPRESS: str << HERMES_WAIT_KEYPRESS_MSG; break;
          default: throw Hermes::Exceptions::Exception("Unknown wait event"); break;
          }
        }
        str << " >>" << std::endl;

        //do something
        switch(wait_event)
        {
        case HERMES_WAIT_CLOSE: wait_for_all_views_close(str.str().c_str()); break;
        case HERMES_WAIT_KEYPRESS: wait_for_any_key(str.str().c_str()); break;
        default: throw Hermes::Exceptions::Exception("Unknown wait event"); break;
        }
      }

      void View::refresh()
      {
        bool do_refresh = true;
        view_sync.enter();
        if(output_id < 0)
          do_refresh = false;
        view_sync.leave();
        if(do_refresh)
          refresh_view(output_id);
      }

      void View::reset_view(bool force_reset)
      {
        if(force_reset || view_not_reset)
        {
          double mesh_width  = vertices_max_x - vertices_min_x;
          double mesh_height = vertices_max_y - vertices_min_y;

          double usable_width = output_width - 2*margin - lspace - rspace;
          double usable_height = output_height - 2*margin;

          // align in the proper direction
          if(usable_width / usable_height < mesh_width / mesh_height)
            scale = usable_width / mesh_width;
          else
            scale = usable_height / mesh_height;
          log_scale = log(scale) / log(H2DV_SCALE_LOG_BASE);

          // center of the mesh
          trans_x = -scale * (vertices_min_x + vertices_max_x) / 2;
          trans_y = -scale * (vertices_min_y + vertices_max_y) / 2;

          view_not_reset = false;
        }
      }

      void View::on_create(int output_id)
      {
        this->output_id = output_id; //does not need thread protection because it is during execution of add_view
        create_gl_palette();
        set_palette_filter(pal_filter == GL_LINEAR);
      }

      void View::on_close()
      {
        view_sync.enter();
        output_id = -1;
        view_sync.leave();
      }

      void View::clear_background()
      {
        glClearColor(1.0, 1.0, 1.0, 0.0);
        glClear(GL_COLOR_BUFFER_BIT);
      }

      void View::pre_display()
      {
        //this->info("display: lock");
        view_sync.enter();

        //begin time measuring
        double time_start = get_tick_count();

        //antialising is supported through accumulation buffer (FIXME: use ARB_MULTISAMPLE if available)
        if(!hq_frame)
        {
          clear_background();
          on_display();
        }
        else
        {
          display_antialiased();
          hq_frame = false;
        }

        if(b_help) draw_help();
        else if(b_scale) scale_dispatch();

        //wait to finish
        glFinish();

        //calculate statistics
        rendering_frames[rendering_frames_top] = get_tick_count() - time_start;
        rendering_frames_top = (rendering_frames_top + 1) % FPS_FRAME_SIZE;

        if(want_screenshot)
        {
          glReadBuffer(GL_BACK_LEFT);
          save_screenshot_internal(screenshot_filename.c_str());
          want_screenshot = false;
        }

        glutSwapBuffers();

        //frame synchronization
        frame_ready = true;
        view_sync.signal_drawing_finished();
        view_sync.leave();
      }

      static float jitter16[16][2] =
      {
        { 0.4375, 0.4375 }, { 0.1875, 0.5625 },
        { 0.9375, 1.1875 }, { 0.4375, 0.9375-1 },
        { 0.6875, 0.5625 }, { 0.1875, 0.0625 },
        { 0.6875, 0.3125 }, { 0.1875, 0.3125 },
        { 0.4375, 0.1875 }, { 0.9375-1, 0.4375 },
        { 0.6875, 0.8125 }, { 0.4375, 0.6875 },
        { 0.6875, 0.0625 }, { 0.9375, 0.9375 },
        { 1.1875, 0.8125 }, { 0.9375, 0.6875 }
      };

      void View::display_antialiased()
      {
        glClear(GL_ACCUM_BUFFER_BIT);
        for (int i = 0; i < 16; i++)
        {
          jitter_x = jitter16[i][0];
          jitter_y = jitter16[i][1];
          set_ortho_projection();
          clear_background();
          on_display();
          glAccum(GL_ACCUM, 1.0 / 16);
        }
        glAccum(GL_RETURN, 1.0);
        jitter_x = jitter_y = 0.0;
      }

      void View::set_ortho_projection(bool no_jitter)
      {
        double jx = no_jitter ? 0.0 : jitter_x;
        double jy = no_jitter ? 0.0 : jitter_y;

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(jx, output_width + jx, output_height-1 + jy, -1 + jy, -10, 10);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
      }

      void View::set_3d_projection(int fov, double znear, double zfar)
      {
        double top = znear * tan((double) fov / 2.0 / 180.0 * M_PI);
        double right = (double) output_width / output_height * top;
        double left = -right;
        double bottom = -top;
        double offsx = (right - left) / output_width * jitter_x;
        double offsy = (top - bottom) / output_height * jitter_y;

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glFrustum(left - offsx, right - offsx, bottom - offsy, top - offsy, znear, zfar);
      }

      void View::draw_fps()
      {
        //calculate FPS
        double frame_time_sum = 0;
        for(int i = 0; i < FPS_FRAME_SIZE; i++)
          frame_time_sum += rendering_frames[i];

        //prepare text
        unsigned char buffer[128];
        sprintf((char*)buffer, "avg. frame: %.1f ms", (float)(frame_time_sum / FPS_FRAME_SIZE));

        //prepare environment
        void* font = GLUT_BITMAP_HELVETICA_10;
        int width_px = glutBitmapLength(font, buffer);
        int height_px = glutBitmapHeight(font);
        int edge_thickness = 2;
        set_ortho_projection(false);

        //render background
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
        glBegin(GL_QUADS);
        glColor4f(1.0f, 1.0f, 1.0f, 0.5f);
        glVertex2i(output_width - (width_px + 2*edge_thickness), 0);
        glVertex2i(output_width, 0);
        glVertex2i(output_width, height_px + 2*edge_thickness);
        glVertex2i(output_width - (width_px + 2*edge_thickness), height_px + 2*edge_thickness);
        glEnd();

        // render text
        glDisable(GL_BLEND);
        glColor3f(1.0f, 0.0f, 0.0f);
        glRasterPos2i(output_width - (width_px + edge_thickness), edge_thickness + height_px);
        // iming information is printed into the image.
        glutBitmapString(font, buffer);
      }

      void View::on_reshape(int width, int height)
      {
        glViewport(0, 0, width, height);

        output_width = width;
        output_height = height;
        update_layout();
      }

      void View::on_mouse_move(int x, int y)
      {
        if(dragging)
        {
          trans_x += (x - mouse_x);
          trans_y += (mouse_y - y);
          refresh();
        }
        else if(scaling)
        {
          log_scale += (mouse_y - y);
          scale = pow(H2DV_SCALE_LOG_BASE, log_scale);
          trans_x = scx - objx * scale - center_x;
          trans_y = center_y - objy * scale - scy;
          refresh();
        }
        else if(scale_dragging)
        {
          int oldv = pos_vert, oldh = pos_horz;
          pos_horz = (x > output_width/2);
          pos_vert = (y < output_height/2);
          if(pos_horz != oldh || pos_vert != oldv)
          {
            update_layout();
            refresh();
          }
        }
        else
        {
          bool oldf = scale_focused;
          scale_focused = (x >= scale_x && x <= scale_x + scale_width &&
            y >= scale_y && y <= scale_y + scale_height);
          if(oldf != scale_focused)
            refresh();
        }
        mouse_x = x;
        mouse_y = y;
      }

      void View::on_left_mouse_down(int x, int y)
      {
        if(scale_focused)
          scale_dragging = true;
        else
          dragging = true;
        scaling = false;
        mouse_x = x;
        mouse_y = y;
      }

      void View::on_left_mouse_up(int x, int y)
      {
        scaling = dragging = scale_dragging = false;
        on_mouse_move(x, y);
      }

      void View::on_right_mouse_down(int x, int y)
      {
        scaling = true;
        dragging = false;
        scx = x;
        scy = y;
        objx = (x - center_x - trans_x) / scale;
        objy = (center_y - y - trans_y) / scale;
        mouse_x = x;
        mouse_y = y;
      }

      void View::on_right_mouse_up(int x, int y)
      {
        scaling = dragging = false;
      }

      void View::on_key_down(unsigned char key, int x, int y)
      {
        const char *file_name = NULL;

        switch (key)
        {
        case 'h':
          {
            hq_frame = true;
            refresh();
            break;
          }

        case 27:
        case 'q':
          {
            close();
            break;
          }

        case 's':
          {
            const char *file_name = get_screenshot_file_name();
            glReadBuffer(GL_FRONT_LEFT);
            save_screenshot_internal(file_name);
            break;
          }

        case 'p':
          {
            // There used to be a type called default, but it caused some weird behavior.
            /*
            switch(pal_type)
            {
            case H2DV_PT_DEFAULT: pal_type = H2DV_PT_HUESCALE; break;
            case H2DV_PT_HUESCALE: pal_type = H2DV_PT_GRAYSCALE; break;
            case H2DV_PT_GRAYSCALE: pal_type = H2DV_PT_INVGRAYSCALE; break;
            case H2DV_PT_INVGRAYSCALE: pal_type = H2DV_PT_DEFAULT; break;
            default: throw Hermes::Exceptions::Exception("Invalid palette type");
            }
            */
            switch(pal_type)
            {
            case H2DV_PT_HUESCALE: pal_type = H2DV_PT_GRAYSCALE; break;
            case H2DV_PT_GRAYSCALE: pal_type = H2DV_PT_INVGRAYSCALE; break;
            case H2DV_PT_INVGRAYSCALE: pal_type = H2DV_PT_HUESCALE; break;
            default: throw Hermes::Exceptions::Exception("Invalid palette type");
            }
            create_gl_palette();
            refresh();
            break;
          }

        default:
          view_sync.enter();
          view_sync.signal_keypress();
          view_sync.leave();
          break;
        }
      }

      void View::on_special_key(int key, int x, int y)
      {
        switch (key)
        {
        case GLUT_KEY_F1:
          b_help = !b_help;
          refresh();
          break;
        }
      }

      void View::wait_for_keypress(const char* text)
      {
        this->warn("Function View::wait_for_keypress deprecated: use View::wait instead");
        View::wait(HERMES_WAIT_KEYPRESS, text);
      }

      void View::wait_for_close()
      {
        view_sync.enter();
        if(output_id >= 0)
          view_sync.wait_close();
        view_sync.leave();
      }

      void View::wait_for_draw()
      {
        // For some reason, this function removes the signal handlers. So we just
        // remember them and restore them. Unfortunately, this doesn't work for some
        // reason:
        //sighandler_t old_segv, old_abrt;
        //old_segv = signal(SIGSEGV, SIG_DFL);
        //old_abrt = signal(SIGABRT, SIG_DFL);

        view_sync.enter();
        if(output_id >= 0 && !frame_ready)
          view_sync.wait_drawing_fisnihed();
        view_sync.leave();

        // Restore the old signal handlers -- doesn't work for some reason:
        //signal(SIGSEGV, old_segv);
        //signal(SIGABRT, old_abrt);
        // So we just restore it by calling the original handler:
      }

      double View::get_tick_count()
      {
#ifdef WIN32
        LARGE_INTEGER freq, ticks;
        QueryPerformanceFrequency(&freq);
        QueryPerformanceCounter(&ticks);
        return (1000.0 * (double)ticks.QuadPart) / (double)freq.QuadPart;
#else
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return (double) tv.tv_sec * 1000 + (double) tv.tv_usec / 1000;
#endif
      }

      void View::set_title(const char* title)
      {
        bool do_set_title = true;

        // Always set the title property.
        this->title = title;

        view_sync.enter();
        if(output_id < 0)
          // If the window does not exist, do nothing else and wait until it is created.
          do_set_title = false;

        view_sync.leave();

        // If the window already exists, show the new title in its header.
        if(do_set_title)
          set_view_title(output_id, title);
      }

      void View::get_palette_color(double x, float* gl_color)
      {
        if(pal_type == H2DV_PT_HUESCALE || pal_type == H2DV_PT_DEFAULT) { //default color
          if(x < 0.0) x = 0.0;
          else if(x > 1.0) x = 1.0;
          x *= num_pal_entries;
          int n = (int)x;
          gl_color[0] = palette_data[n][0];
          gl_color[1] = palette_data[n][1];
          gl_color[2] = palette_data[n][2];
        }
        else if(pal_type == H2DV_PT_GRAYSCALE)
          gl_color[0] = gl_color[1] = gl_color[2] = (float)x;
        else if(pal_type == H2DV_PT_INVGRAYSCALE)
          gl_color[0] = gl_color[1] = gl_color[2] = (float)(1.0 - x);
        else
          gl_color[0] = gl_color[1] = gl_color[2] = 1.0f;
      }

      void View::set_num_palette_steps(int num)
      {
        if(num < 2) num = 2;
        if(num > 256) num = 256;
        pal_steps = num;
        update_tex_adjust();

        view_sync.enter();
        if(output_id >= 0)
          create_gl_palette();
        view_sync.leave();

        refresh();
      }

      void View::create_gl_palette()
      {
        int i;
        unsigned char palette[256][3];
        for (i = 0; i < pal_steps; i++)
        {
          float gl_color[3];
          get_palette_color((double) i / pal_steps, gl_color);
          palette[i][0] = (unsigned char) (gl_color[0] * 255);
          palette[i][1] = (unsigned char) (gl_color[1] * 255);
          palette[i][2] = (unsigned char) (gl_color[2] * 255);
        }
        for (i = pal_steps; i < 256; i++)
          memcpy(palette[i], palette[pal_steps-1], 3);

        if(gl_pallete_tex_id == 0)
          glGenTextures(1, &gl_pallete_tex_id);
        glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
        glTexImage1D(GL_TEXTURE_1D, 0, 3, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, palette);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
      }

      void View::set_palette_filter(bool linear)
      {
        view_sync.enter(); //lock to prevent simultaneuous rendering

        pal_filter = linear ? GL_LINEAR : GL_NEAREST;

        if(gl_pallete_tex_id == 0)
          glGenTextures(1, &gl_pallete_tex_id);
        glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, pal_filter);
        glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, pal_filter);
        update_tex_adjust();

        view_sync.leave(); //unlock

        refresh();
      }

      void View::set_palette(ViewPaletteType type)
      {
        view_sync.enter();
        pal_type = type;
        if(output_id >= 0)
          create_gl_palette();
        view_sync.leave();

        //redisplay
        refresh();
      }

      void View::update_tex_adjust()
      {
        if(pal_filter == GL_LINEAR)
        {
          tex_scale = (double) (pal_steps-1) / 256.0;
          tex_shift = 0.5 / 256.0;
        }
        else
        {
          tex_scale = (double) pal_steps / 256.0;
          tex_shift = 0.0;
        }
      }

      void View::set_min_max_range(double min, double max)
      {
        if(max < min)
        {
          std::swap(min, max);
          this->warn("Upper bound set below the lower bound: reversing to (%f, %f).", min, max);
        }
        view_sync.enter();
        range_min = min;
        range_max = max;
        range_auto = false;
        if(output_id >= 0)
          update_layout();
        view_sync.leave();
        refresh();
      }

      void View::auto_min_max_range()
      {
        view_sync.enter();
        range_auto = true;
        if(output_id >= 0)
          update_layout();
        view_sync.leave();

        refresh();
      }

      void View::get_min_max_range(double& min, double& max)
      {
        view_sync.enter();
        min = range_min;
        max = range_max;
        view_sync.leave();
      }

      void View::draw_text(double x, double y, const char* text, int align)
      {
        void* font = GLUT_BITMAP_9_BY_15;
        if(align > -1)
        {
          int width = glutBitmapLength(font, (const unsigned char*) text);
          if(align == 1) x -= width; // align right
          else x -= (double) width / 2; // center
        }
        y += 5; //(double) glutBitmapHeight(font) / 2 - 1;

        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);

        glRasterPos2d((int) (x + 0.5), (int) (y + 0.5));
        glutBitmapString(font, (const unsigned char*) text);
      }

      int View::get_text_width(const char* text)
      {
        void* font = GLUT_BITMAP_9_BY_15;
        return glutBitmapLength(font, (const unsigned char*) text);
      }

      void View::draw_help()
      {
        view_sync.enter();
        set_ortho_projection(true);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glDisable(GL_TEXTURE_1D);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        const char* text = get_help_text();

        int n = 1;
        for (const char* p = text; *p; p++)
          if(*p == '\n') n++;

        int width = get_text_width(text);
        int height = n * glutBitmapHeight(GLUT_BITMAP_9_BY_15);
        int x = 10, y = 10, b = 6;

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(1.0f, 1.0f, 1.0f, 0.65f);
        glBegin(GL_QUADS);
        glVertex2d(x, y + height + 2*b);
        glVertex2d(x + width + 2*b, y + height + 2*b);
        glVertex2d(x + width + 2*b, y);
        glVertex2d(x, y);
        glEnd();

        glDisable(GL_BLEND);
        glColor3f(0, 0, 0);
        draw_text(x + b, y + b+7, text);
        view_sync.leave();
        refresh();
      }

      char* View::get_screenshot_file_name()
      {
        static char file_name[1024] = { 0 };
        bool got_file_name = false;
        do
        {
          sprintf(file_name, "screen%03d.bmp", screenshot_no);
          FILE *f = fopen(file_name, "r");
          if(f == NULL)
            got_file_name = true;
          else
            fclose(f);
          screenshot_no++;
        }
        while (!got_file_name);
        return file_name;
      }

      typedef unsigned int dword;
      typedef unsigned short word;

      const word BITMAP_ID = 0x4D42;

#pragma pack(1)

      struct BitmapFileHeader
      {
        word  type;
        dword size;
        word  reserved1;
        word  reserved2;
        dword off_bits;
      };

      struct BitmapInfoHeader
      {
        dword size;
        dword width;
        dword height;
        word  planes;
        word  bit_count;
        dword compression;
        dword size_image;
        dword xdpi;
        dword ydpi;
        dword clr_used;
        dword clr_important;
      };

      void View::save_screenshot_internal(const char *file_name)
      {
        BitmapFileHeader file_header;
        BitmapInfoHeader info_header;

        // alloc memory for pixel data (4 bytes per pixel)
        char* pixels = NULL;
        if((pixels = (char*) malloc(4 * output_width * output_height)) == NULL)
          throw Hermes::Exceptions::Exception("Could not allocate memory for pixel data");

        // get pixels from framebuffer
#ifdef GL_BGRA_EXT
        glReadPixels(0, 0, output_width, output_height, GL_BGRA_EXT, GL_UNSIGNED_BYTE, pixels);
#else
        glReadPixels(0, 0, output_width, output_height, GL_RGBA, GL_UNSIGNED_BYTE, pixels); // FIXME!!!
        this->warn("BGRA format not supported. Saved image will have inverted colors");
#endif
        // opening file for binary writing
        FILE* file = fopen(file_name, "wb");
        if(file == NULL)
          throw Hermes::Exceptions::Exception("Could not open '%s' for writing", file_name);

        // fill in bitmap header
        file_header.type = BITMAP_ID;
        file_header.size = sizeof(BitmapFileHeader) +  sizeof(BitmapInfoHeader) +
          4 * output_width * output_height;
        file_header.reserved1 = file_header.reserved2 = 0;
        file_header.off_bits = 14 + 40; // length of both headers

        if(fwrite(&file_header, sizeof(file_header), 1, file) != 1)
        {
          fclose(file);
          throw Hermes::Exceptions::Exception("Error writing bitmap header");
        }

        // fill in bitmap info header
        info_header.size = sizeof(BitmapInfoHeader);
        info_header.width = output_width;
        info_header.height = output_height;
        info_header.planes = 1;
        info_header.bit_count = 32; // 4 bytes per pixel = 32 bits
        info_header.compression = 0;
        info_header.size_image = output_width * output_height * 4;
        info_header.xdpi = 2835; // 72 dpi
        info_header.ydpi = 2835; // 72 dpi
        info_header.clr_used = 0;
        info_header.clr_important = 0;

        if(fwrite(&info_header, sizeof(info_header), 1, file) != 1)
        {
          fclose(file);
          throw Hermes::Exceptions::Exception("Error writing bitmap header");
        }

        // write image pixels
        if(fwrite((GLubyte*) pixels, 1, info_header.size_image, file) != info_header.size_image)
        {
          fclose(file);
          throw Hermes::Exceptions::Exception("Error writing pixel data");
        }

        fclose(file);
        free((void*) pixels);
        printf("Image \"%s\" saved.\n", file_name);
      }

      void View::save_screenshot(const char* bmpname, bool high_quality)
      {
        view_sync.enter();
        if(output_id >= 0) { //set variable neccessary to create a screenshot
          hq_frame = high_quality;
          want_screenshot = true;
          screenshot_filename = bmpname;
        }
        view_sync.leave();

        //request redraw
        refresh();
      }

      void View::save_numbered_screenshot(const char* format, int number, bool high_quality)
      {
        char buffer[1000];
        sprintf(buffer, format, number);
        save_screenshot(buffer, high_quality);
      }

      int View::measure_scale_labels()
      {
        int result = 0;
        for (int i = 0; i <= scale_numticks + 1; i++)
        {
          double value = range_min + (double) i * (range_max - range_min) / (scale_numticks + 1);
          if(fabs(value) < 1e-8) value = 0.0;
          char text[50];
          sprintf(text, scale_fmt, value);
          int w = get_text_width(text);
          if(w > result) result = w;
        }
        return result;
      }

      void View::draw_continuous_scale(char* title, bool righttext)
      {
        int i;
        double y0 = scale_y + scale_height;

        set_ortho_projection(true);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glDisable(GL_TEXTURE_1D);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // background
        const int b = 5;
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(1.0f, 1.0f, 1.0f, 0.65f);
        int rt = righttext ? 0 : labels_width + 8;
        glBegin(GL_QUADS);
        glVertex2d(scale_x - b - rt, y0 + 5 + b);
        glVertex2d(scale_x + scale_width + 8 + labels_width + b - rt, y0 + 5 + b);
        glVertex2d(scale_x + scale_width + 8 + labels_width + b - rt, scale_y - 5 - b);
        glVertex2d(scale_x - b - rt, scale_y - 5 - b);
        glEnd();

        // palette
        glDisable(GL_BLEND);
        glColor3f(0.0f, 0.0f, 0.0f);
        glBegin(GL_QUADS);
        glVertex2d(scale_x, scale_y);
        glVertex2d(scale_x, scale_y + scale_height + 1);
        glVertex2d(scale_x + scale_width + 1, scale_y + scale_height + 1);
        glVertex2d(scale_x + scale_width + 1, scale_y);
        glEnd();

        glEnable(GL_TEXTURE_1D);
        glBindTexture(GL_TEXTURE_1D, gl_pallete_tex_id);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
        glBegin(GL_QUADS);
        glTexCoord1d(tex_scale + tex_shift);
        glVertex2d(scale_x + 1, scale_y + 1);
        glVertex2d(scale_x + scale_width, scale_y + 1);
        glTexCoord1d(tex_shift);
        glVertex2d(scale_x + scale_width, scale_y + scale_height);
        glVertex2d(scale_x + 1, scale_y + scale_height);
        glEnd();

        // focus
        glDisable(GL_TEXTURE_1D);
        if(scale_focused)
        {
          glEnable(GL_BLEND);
          glColor4f(1.0f, 1.0f, 1.0f, 0.3f);
          glBegin(GL_QUADS);
          glVertex2d(scale_x + 1, scale_y + 1);
          glVertex2d(scale_x + scale_width, scale_y + 1);
          glVertex2d(scale_x + scale_width, scale_y + scale_height);
          glVertex2d(scale_x + 1, scale_y + scale_height);
          glEnd();
        }

        // ticks
        glColor3f(0, 0, 0);
        glDisable(GL_BLEND);
        glDisable(GL_LINE_STIPPLE);
        glLineWidth(1.0);
        glBegin(GL_LINES);
        for (i = 0; i < scale_numticks; i++)
        {
          y0 = scale_y + scale_height - (double) (i + 1) * scale_height / (scale_numticks + 1);
          glVertex2d(scale_x, y0);
          glVertex2d(scale_x + 0.2 * scale_width + 1, y0);
          glVertex2d(scale_x + 0.8 * scale_width, y0);
          glVertex2d(scale_x + scale_width, y0);
        }
        glEnd();

        // labels
        for (i = 0; i <= scale_numticks + 1; i++)
        {
          double value = range_min + (double) i * (range_max - range_min) / (scale_numticks + 1);
          if(fabs(value) < 1e-8) value = 0.0;
          char text[50];
          sprintf(text, scale_fmt, value);
          y0 = scale_y + scale_height - (double) i * scale_height / (scale_numticks + 1);
          if(righttext)
            draw_text(scale_x + scale_width + 8, y0, text);
          else
            draw_text(scale_x - 8, y0, text, 1);
        }
      }

      void View::draw_discrete_scale(int numboxes, const char* boxnames[], const float boxcolors[][3])
      {
        set_ortho_projection(true);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glDisable(GL_TEXTURE_1D);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // background
        const int b = 5;
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(1.0f, 1.0f, 1.0f, 0.65f);
        glBegin(GL_QUADS);
        glVertex2d(scale_x - b, scale_y - b);
        glVertex2d(scale_x - b, scale_y + scale_height + b + 1);
        glVertex2d(scale_x + scale_width + b + 1, scale_y + scale_height + b + 1);
        glVertex2d(scale_x + scale_width + b + 1, scale_y - b);
        glEnd();

        // boxes
        glDisable(GL_BLEND);
        int y = scale_y;
        for (int i = 0; i < numboxes; i++)
        {
          glColor3f(0.0, 0.0, 0.0);
          glBegin(GL_QUADS);
          glVertex2d(scale_x, y);
          glVertex2d(scale_x, y + scale_box_height + 1);
          glVertex2d(scale_x + scale_width + 1, y + scale_box_height + 1);
          glVertex2d(scale_x + scale_width + 1, y);
          glEnd();

          const float* color = boxcolors[numboxes-1-i];
          float bcolor[3] = { color[0], color[1], color[2] };
          if(scale_focused)
          {
            bcolor[0] = color[0]*0.7f + 1.0f*0.3f;
            bcolor[1] = color[1]*0.7f + 1.0f*0.3f;
            bcolor[2] = color[2]*0.7f + 1.0f*0.3f;
          }

          glColor3f(bcolor[0], bcolor[1], bcolor[2]);
          glBegin(GL_QUADS);
          glVertex2d(scale_x + 1, y + 1);
          glVertex2d(scale_x + 1, y + scale_box_height);
          glVertex2d(scale_x + scale_width, y + scale_box_height);
          glVertex2d(scale_x + scale_width, y + 1);
          glEnd();

          if((color[0] + color[1] + color[2]) / 3 > 0.5)
            glColor3f(0, 0, 0);
          else
            glColor3f(1, 1, 1);

          int a = scale_x + scale_width/2;
          int b = y + scale_box_height/2;
          draw_text(a, b, boxnames[numboxes-1-i], 0);
          draw_text(a + 1, b, boxnames[numboxes-1-i], 0);

          y += scale_box_height + scale_box_skip;
        }
      }

      void View::scale_dispatch()
      {
        draw_continuous_scale(NULL, !pos_horz);
      }

      void View::update_layout()
      {
        lspace = rspace = labels_width = 0;
        if(b_scale)
        {
          labels_width = scale_fixed_width;
          if(labels_width < 0) labels_width = measure_scale_labels();
          int space = scale_width + 8 + labels_width + margin;
          if(pos_horz == 0)
          { lspace = space;  scale_x = margin; }
          else
          { rspace = space;  scale_x = output_width - margin - scale_width; }

          if(pos_vert == 0)
            scale_y = output_height - margin - scale_height;
          else
            scale_y = margin;
        }

        center_x = ((double) output_width - 2*margin - lspace - rspace) / 2 + margin + lspace;
        center_y = (double) output_height / 2;
      }

      void View::show_scale(bool show)
      {
        view_sync.enter();
        b_scale = show;
        if(output_id >= 0)
          update_layout();
        view_sync.leave();
        refresh();
      }

      void View::set_scale_position(int horz, int vert)
      {
        view_sync.enter();
        pos_horz = horz;
        pos_vert = vert;
        if(output_id >= 0)
          update_layout();
        view_sync.leave();
        refresh();
      }

      void View::set_scale_size(int width, int height, int numticks)
      {
        view_sync.enter();
        scale_width = width;
        scale_height = height;
        scale_numticks = numticks;
        update_layout();
        view_sync.leave();
        refresh();
      }

      void View::set_scale_format(const char* fmt)
      {
        view_sync.enter();
        strncpy(scale_fmt, fmt, 19);
        update_layout();
        view_sync.leave();
        refresh();
      }

      void View::fix_scale_width(int width)
      {
        view_sync.enter();
        scale_fixed_width = width;
        update_layout();
        view_sync.leave();
        refresh();
      }
    }
  }
}
#endif