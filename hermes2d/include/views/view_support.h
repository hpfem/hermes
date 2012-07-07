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

#ifndef __H2D_VIEW_SUPPORT_H
#define __H2D_VIEW_SUPPORT_H

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      // you can define NOGLUT to turn off all OpenGL stuff in Hermes2D
#ifndef NOGLUT

      /* types */
      class HERMES_API ViewMonitor ///< A monitor used to synchronize thread in views.
      {
      protected:
        pthread_mutexattr_t mutex_attr; ///< Mutext attributes.
        pthread_mutex_t mutex; ///< Mutex that protects monitor.
        pthread_cond_t cond_cross_thread_call; ///< Condition used to signal a cross-thread call
        pthread_cond_t cond_keypress; ///< Condition used to signal a keypress.
        pthread_cond_t cond_close; ///< Condition used to signal close of a window.
        pthread_cond_t cond_drawing_finished; ///< Condition used to signal that drawing has finished
      public:
        ViewMonitor();
        ~ViewMonitor();
        inline void enter() { pthread_mutex_lock(&mutex); }; ///< enters protected section
        inline void leave() { pthread_mutex_unlock(&mutex); }; ///< leaves protected section
        inline void signal_keypress() { pthread_cond_broadcast(&cond_keypress); }; ///< signals keypress inside a protected section
        inline void wait_keypress() { pthread_cond_wait(&cond_keypress, &mutex); }; ///< waits for keypress inside a protected section
        inline void signal_close() { pthread_cond_broadcast(&cond_close); }; ///< signals close inside a protected section
        inline void wait_close() { pthread_cond_wait(&cond_close, &mutex); }; ///< waits for close inside a protected section
        inline void signal_drawing_finished() { pthread_cond_broadcast(&cond_drawing_finished); }; ///< signals drawing finished inside a protected section
        inline void wait_drawing_fisnihed() { pthread_cond_wait(&cond_drawing_finished, &mutex); }; ///< waits for drawing finished inside a protected section
        inline void signal_cross_thread_call() { pthread_cond_broadcast(&cond_cross_thread_call); }; ///< signals that cross-thread-call finished
        inline void wait_cross_thread_call() { pthread_cond_wait(&cond_cross_thread_call, &mutex); }; ///< waits for finishing of a cross-thread call
      };

      /* exported types */
      extern ViewMonitor view_sync; ///< synchronization between all views. Used to access OpenGL and signal a window close event and a keypress event.

      /* exported functions */
      class View;

      HERMES_API bool init_glut(); ///< Initialize GLUT.
      HERMES_API bool shutdown_glut(); ///< Shutdown GLUT.
      extern "C"
      {
        int add_view(View* view, int x, int y, int width, int height, const char* title); ///< Adds a view.
      }
      extern void set_view_title(int view_id, const char* title); ///< Sets title of a view.
      extern void refresh_view(int view_id); ///< Forces redisplay of a view.
      extern void remove_view(int view_id); ///< Removes a view.
      //extern void force_view_thread_shutdown(); ///< Forces view thread to shutdown.
      extern void wait_for_all_views_close(const char* text); ///< Waits for all views to close.
      extern void wait_for_any_key(const char* text); ///< Waits for a keypress which is not processed.

      extern void on_display_stub(void);
      extern void on_reshape_stub(int width, int height);
      extern void on_mouse_move_stub(int x, int y);
      extern void on_key_down_stub(unsigned char key, int x, int y);
      extern void on_special_key_stub(int key, int x, int y);
      extern void on_entry_stub(int state);
      extern void on_mouse_click_stub(int button, int state, int x, int y);
      extern void on_close_stub();

#endif
    }
  }
}
#endif