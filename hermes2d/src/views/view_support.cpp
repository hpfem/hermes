// This file is part of Hermes2D.
//
// Copyright 2009 Ivo Hanak <hanak@byte.cz>
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

#ifndef NOGLUT

#include <GL/glew.h>
#include <GL/freeglut.h>
#ifndef WIN32
# include <sys/time.h>
#endif

#include <map>
#include <vector>
#include "global.h"
#include "view.h"
#include "view_support.h"

namespace Hermes
{
  namespace Hermes2D
  {
    namespace Views
    {
      ViewMonitor view_sync; ///< Synchronization.

      /* monitor */
      ViewMonitor::ViewMonitor()
      {
        pthread_mutexattr_init(&mutex_attr);
        pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_RECURSIVE);
        pthread_mutex_init(&mutex, &mutex_attr);

        pthread_cond_init(&cond_keypress, NULL);
        pthread_cond_init(&cond_close, NULL);
        pthread_cond_init(&cond_drawing_finished, NULL);
        pthread_cond_init(&cond_cross_thread_call, NULL);
      }

      ViewMonitor::~ViewMonitor()
      {
        pthread_mutex_destroy(&mutex);
        pthread_mutexattr_destroy(&mutex_attr);
        pthread_cond_destroy(&cond_keypress);
        pthread_cond_destroy(&cond_close);
        pthread_cond_destroy(&cond_drawing_finished);
        pthread_cond_destroy(&cond_cross_thread_call);
      }

      /* forward definitions */
      static bool glut_initialized = false; ///< True if GLUT is initialized
      static bool glew_initialized = false; ///< True if GLEW is initialized

      /* threading */
      struct ThreadInfo { ///< Thread information
        pthread_t thread; ///< View thread.
        bool should_quit; ///< True if view thread should end
        ThreadInfo() : should_quit(false) {};
      };

      static ThreadInfo* view_thread; ///< Current view thread.

      typedef int (*CTC_FUNC)(void*);
      static CTC_FUNC ctc_function = NULL;
      static void* ctc_param;
      static int ctc_result;
      static std::map<int, View*> view_instances; ///< Instances of views.

      struct ViewParams { ///< Parameters used to initialize view
        View* view;
        int x, y;
        int width, height;
        const char* title;
        ViewParams(View* view, int x, int y, int width, int height, const char* title) : view(view), x(x), y(y), width(width), height(height), title(title) {};
      };

      struct RemoveParams { ///< Parameters used to remove view
        int view_id;
        bool destroy_glut_window; ///< True to destroy glut window.
        RemoveParams(int view_id, bool destroy_glut_window) : view_id(view_id), destroy_glut_window(destroy_glut_window) {};
      };

      struct TitleParams { ///< Parameters of a title
        int view_id;
        const char* title;
        TitleParams(int view_id, const char* title) : view_id(view_id), title(title) {};
      };

      /* view thread function */
      static void* view_thread_func(void* param)
      {
        ThreadInfo* thread_info = (ThreadInfo*)param;

        //initize GLUT
        init_glut();

        //run message loop
        while(1)
        {
          //handle glut messages
          glutMainLoopEvent();

          //check whether to quit
          if(thread_info->should_quit)
            break;

          //handle CTC
          view_sync.enter();
          if(ctc_function != NULL)
          {
            ctc_result = ctc_function(ctc_param);
            ctc_function = NULL;
            view_sync.signal_cross_thread_call();
          }
          view_sync.leave();

          //wait for 10 ms
#ifdef WIN32
          Sleep(10);
#else
          usleep(10*1000);
#endif
        }

        //cleanup
        shutdown_glut();

        //cleanup
        delete thread_info;
        return NULL;
      }

      /// Returns true, if outside the thread.
      static bool need_safe_call() { return pthread_equal(view_thread->thread, pthread_self()) == 0; }

      /// Does a cross thread call.
      static int call_in_thread(CTC_FUNC func, void* param)
      {
        //check whether the thread is running if not start it
        view_sync.enter();
        if(view_thread == NULL)
        {
          ThreadInfo* new_thread_info = NULL;
          try { new_thread_info = new ThreadInfo(); }
          catch(std::bad_alloc&) { throw Hermes::Exceptions::Exception("Failed to allocate structure for view thread"); }
          int err = pthread_create(&new_thread_info->thread, NULL, view_thread_func, new_thread_info);
          if(err)
          {
            delete new_thread_info;
            throw Hermes::Exceptions::Exception("Failed to create main thread, error: %d", err);
          }
          view_thread = new_thread_info;
        }
        view_sync.leave();

        //make call
        if(need_safe_call())
        {
          view_sync.enter();
          ctc_function = func;
          ctc_param = param;
          view_sync.wait_cross_thread_call();
          int result = ctc_result;
          view_sync.leave();
          return result;
        }
        else
          return func(param);
      }

      /// Sets a title of a view. Function has to be called just from the inside of view thread with a locked sync_view.
      static int set_view_title_in_thread(void* title_pars_ptr)
      {
        TitleParams& title_params = *((TitleParams*)title_pars_ptr);
        std::map<int, View*>::iterator found_view = view_instances.find(title_params.view_id);
        if(found_view == view_instances.end())
        {
          throw Exceptions::Exception("Settings title of a view that is not registered.");
          return -1;
        }

        //create GLUT window
        glutSetWindow(title_params.view_id);
        glutSetWindowTitle(title_params.title);

        return 0;
      }

      /// Adds a new view. Function has to be called just from the inside of view thread with a locked sync_view.
      int add_view_in_thread(void* view_pars_ptr)
      {
        ViewParams& view_params = *((ViewParams*)view_pars_ptr);

        //create GLUT window
        glutInitWindowPosition(view_params.x, view_params.y);
        glutInitWindowSize(view_params.width, view_params.height);
        int view_id = glutCreateWindow(view_params.title);
        glutSetWindowData(view_params.view);

        //initialize GLEW
        GLenum err = glewInit();
        if(err != GLEW_OK)
          throw Exceptions::Exception("GLEW error: %s", glewGetErrorString(err));
        glew_initialized = true;

        //register callbacks
        glutDisplayFunc(on_display_stub);
        glutReshapeFunc(on_reshape_stub);
        glutMotionFunc(on_mouse_move_stub);
        glutPassiveMotionFunc(on_mouse_move_stub);
        glutMouseFunc(on_mouse_click_stub);
        glutKeyboardFunc(on_key_down_stub);
        glutSpecialFunc(on_special_key_stub);
        glutEntryFunc(on_entry_stub);
        glutCloseFunc(on_close_stub);

        //add to structures
        view_instances[view_id] = view_params.view;

        //call handlers
        view_params.view->on_create(view_id);

        return view_id;
      }

      /// Removes a new view. Function has to be called just from the inside of view thread with a locked sync_view.
      int remove_view_in_thread(void* remove_params_ptr)
      {
        RemoveParams& params = *(RemoveParams*)remove_params_ptr;
        std::map<int, View*>::iterator found_view = view_instances.find(params.view_id);
        if(found_view == view_instances.end())
        {
          throw Exceptions::Exception("Removing of a view that is not registered");
          return -1;
        }

        //destroy window if requested (it will not be requested when remove is called as a reaction to on_close)
        if(params.destroy_glut_window)
        {
          //remove window from GLUT
          glutSetWindow(params.view_id);
          glutSetWindowData(NULL); //prevent stubs from being executed if there is still some message waiting for the window

          //call on-close event
          found_view->second->on_close();

          //finish removal of window from GLUT
          glutDestroyWindow(params.view_id);
        }

        //remove from structures
        view_instances.erase(found_view);

        //thread cleanup
        if(view_instances.empty())
        {
          view_thread->should_quit = true;
          view_thread = NULL;

          //signal all events
          view_sync.signal_close();
          view_sync.signal_keypress();
          view_sync.signal_drawing_finished();
        }

        return 0;
      }

      /// Forces a redisplay of a view. Function has to be called just from the inside of view thread.
      static int refresh_view_in_thread(void* view_id_ptr)
      {
        int view_id = *((int*)view_id_ptr);
        std::map<int, View*>::iterator found_view = view_instances.find(view_id);
        if(found_view == view_instances.end())
          throw Exceptions::Exception("Refreshing a view that is not registered");

        //redisplay
        if(found_view != view_instances.end())
        {
          glutSetWindow(view_id);
          glutPostRedisplay();
        }

        return 0;
      }

      /* GLUT */
      static long double_click_delay_ms = 300; ///< Length of the double-click time. (FIXME: get double-click time for Linux)

#define STUB_GET_VIEW() View* wnd = (View*)glutGetWindowData() /* retrieves view for the current GLUT callback */
#define STUB_CALL(__call) STUB_GET_VIEW(); if(wnd != NULL) __call; /* calls a method of a view for the current GLUT callbakc */

      void on_display_stub(void) { STUB_CALL( wnd->pre_display() ); }
      void on_reshape_stub(int width, int height) { STUB_CALL( wnd->on_reshape(width, height) ); }
      void on_mouse_move_stub(int x, int y) { STUB_CALL( wnd->on_mouse_move(x, y) ); }
      void on_key_down_stub(unsigned char key, int x, int y) { STUB_CALL( wnd->on_key_down(key, x, y) ); }
      void on_special_key_stub(int key, int x, int y) { STUB_CALL( wnd->on_special_key(key, x, y) ); }
      void on_entry_stub(int state) { STUB_CALL( wnd->on_entry(state) ); }
      void on_mouse_click_stub(int button, int state, int x, int y)
      {
        STUB_GET_VIEW();
        if(wnd == NULL)
          return;

        // emulate double-click messages
        if(state == GLUT_DOWN)
        {
          static double last_tick = 0;
          double tick = View::get_tick_count();
          //if(tick < last_tick) //todo
          if(tick - last_tick < double_click_delay_ms)
          {
            if(button == GLUT_LEFT_BUTTON)
              wnd->on_left_mouse_double_click(x, y);
            else if(button == GLUT_RIGHT_BUTTON)
              wnd->on_right_mouse_double_click(x, y);
            else
              wnd->on_middle_mouse_double_click(x, y);

            last_tick = 0;
            return;
          }
          last_tick = tick;
        }

        // call proper click handler
        if(button == GLUT_LEFT_BUTTON)
        {
          if(state == GLUT_DOWN)
            wnd->on_left_mouse_down(x, y);
          else
            wnd->on_left_mouse_up(x, y);
        }
        else if(button == GLUT_RIGHT_BUTTON)
        {
          if(state == GLUT_DOWN)
            wnd->on_right_mouse_down(x, y);
          else
            wnd->on_right_mouse_up(x, y);
        }
        else
        {
          if(state == GLUT_DOWN)
            wnd->on_middle_mouse_down(x, y);
          else
            wnd->on_middle_mouse_up(x, y);
        }
      }
      void on_close_stub()
      {
        STUB_GET_VIEW();
        if(wnd == NULL)
          return;

        //call callback
        wnd->on_close();

        //remove view from system
        view_sync.enter();
        RemoveParams params(glutGetWindow(), false);
        remove_view_in_thread(&params);
        view_sync.signal_close();
        view_sync.leave();
      }

      ///initialize GLUT
      bool init_glut()
      {
        static int argc = 1;
        static const char* argv[1] = { "x" };

        //prepare GLUT environment
        if(!glut_initialized)
        {
          glutInit(&argc, (char**) argv);
          glut_initialized = true;
        }
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ACCUM);
        glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

        //obtain other parameters
#ifdef WIN32
        double_click_delay_ms = GetDoubleClickTime();
#endif

        return true;
      }

      /// shutdowns GLUT
      bool shutdown_glut()
      {
        return true;
      }

      /* public functions */
      void set_view_title(int view_id, const char* title)
      {
        TitleParams params(view_id, title);
        call_in_thread(set_view_title_in_thread, &params);
      }

      int add_view(View* view, int x, int y, int width, int height, const char* title)
      {
        ViewParams params(view, x, y, width, height, title);
        return call_in_thread(add_view_in_thread, &params);
      }

      void refresh_view(int view_id) { call_in_thread(refresh_view_in_thread, &view_id); }

      void remove_view(int view_id)
      {
        RemoveParams params(view_id, true);
        call_in_thread(remove_view_in_thread, &params);
      }

      void force_view_thread_shutdown()
      {
        pthread_t current_thread;
        bool should_wait = false;

        view_sync.enter();

        //destroy all views
        std::map<int, View*>::const_iterator iter = view_instances.begin();
        while (iter != view_instances.end())
        {
          glutDestroyWindow(iter->first);
          ++iter;
        }
        view_instances.clear();

        //tell thread to finish
        if(view_thread != NULL)
        {
          current_thread = view_thread->thread;
          view_thread->should_quit = true;
          view_thread = NULL;
        }
        view_sync.leave();

        //wait for thread to finish
        if(should_wait)
        {
          glew_initialized = false;
          pthread_join(current_thread, NULL);
        }
      }

      void wait_for_all_views_close(const char* text)
      {
        pthread_t current_thread;
        bool should_wait = false;

        //tell thread to finish
        view_sync.enter();
        if(view_thread != NULL)
        {
          current_thread = view_thread->thread;
          should_wait = true;
        }
        view_sync.leave();

        //wait for thread to finish
        if(should_wait)
        {
          fprintf(stdout, "%s", text); fflush(stdout);
          pthread_join(current_thread, NULL);
        }
      }

      void wait_for_any_key(const char* text)
      {
        //wait for key
        view_sync.enter();
        if(view_thread != NULL)
        {
          fprintf(stdout, "%s", text); fflush(stdout);
          view_sync.wait_keypress();
        }
        view_sync.leave();
      }
    }
  }
}
#endif