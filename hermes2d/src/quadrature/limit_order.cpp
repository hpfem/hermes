#include "global.h"
#include "quad_all.h"
#include "limit_order.h"

namespace Hermes
{
  namespace Hermes2D
  {
    static int default_order_table_tri[] =
    {
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
      17, 18, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
      20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
      20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
      20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20
    };

#ifdef EXTREME_QUAD
    static int default_order_table_quad[] =
    {
      1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 13, 13, 15, 15,
      17, 17, 19, 19, 21, 21, 23, 23, 25, 25, 27, 27, 29, 29, 31, 31,
      33, 33, 35, 35, 37, 37, 39, 39, 41, 41, 43, 43, 45, 45, 47, 47,
      49, 49, 51, 51, 53, 53, 55, 55, 57, 57, 59, 59, 61, 61, 63, 63,
      65, 65, 67, 67, 69, 69, 71, 71, 73, 73, 75, 75, 77, 77, 79, 79,
      81, 81, 83, 83, 85, 85, 87, 87, 89, 89, 91, 91, 93, 93, 95, 95,
      97, 97, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99, 99
    };
#else
    static int default_order_table_quad[] =
    {
      1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 13, 13, 15, 15, 17,
      17, 19, 19, 21, 21, 23, 23, 24, 24, 24, 24, 24, 24, 24,
      24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
      24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
      24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24
    };
#endif

    static int* g_order_table_quad = default_order_table_quad;
    static int* g_order_table_tri  = default_order_table_tri;
    static bool warned_order = false;

    HERMES_API int  g_max_order;
    HERMES_API int  g_safe_max_order;

    HERMES_API void set_order_limit_table(int* tri_table, int* quad_table, int n)
    {
      if(n < 24) throw Hermes::Exceptions::Exception("Hermes::Order limit tables must have at least 24 entries.");
      g_order_table_tri  = tri_table;
      g_order_table_quad = quad_table;
    }

    HERMES_API void update_limit_table(ElementMode2D mode)
    {
      g_max_order = g_quad_2d_std.get_max_order(mode);
      g_safe_max_order = g_quad_2d_std.get_safe_max_order(mode);
    }

    HERMES_API void reset_warn_order()
    {
      warned_order = false;
    }

    HERMES_API void warn_order()
    {
      if(HermesCommonApi.get_integral_param_value(Hermes::showInternalWarnings))
        if(!warned_order)
        {
  #pragma omp critical (warn_oder)
          if(!warned_order)
          {
            /// \todo Fix this, so that it complies with the rest of the code.
            Hermes::Mixins::Loggable::Static::warn("Warning: Not enough integration rules for exact integration.");
            warned_order = true;
          }
        }
    }

    HERMES_API void limit_order(int& o, ElementMode2D mode)
    {
      if(o > g_quad_2d_std.get_safe_max_order(mode))
      {
        o = g_quad_2d_std.get_safe_max_order(mode);
        warn_order();
      }
      if(mode == HERMES_MODE_TRIANGLE)
        o = g_order_table_tri[o];
      else
        o = g_order_table_quad[o];
    }

    HERMES_API void limit_order_nowarn(int& o, ElementMode2D mode)
    {
      if(o > g_quad_2d_std.get_safe_max_order(mode))
        o = g_quad_2d_std.get_safe_max_order(mode);
      if(mode == HERMES_MODE_TRIANGLE)
        o = g_order_table_tri[o];
      else
        o = g_order_table_quad[o];
    }
  }
}