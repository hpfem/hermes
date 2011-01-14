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

#ifndef __H2D_LIMIT_ORDER_H
#define __H2D_LIMIT_ORDER_H

// can be called to set a custom order limiting table
extern HERMES_API void set_order_limit_table(int* tri_table, int* quad_table, int n);

// limit_order is used in integrals
extern HERMES_API int  g_safe_max_order;
extern HERMES_API int  g_max_order;
extern HERMES_API int* g_order_table;

#ifndef DEBUG_ORDER
  #define limit_order(o) \
    if (o > g_safe_max_order) { o = g_safe_max_order; warn_order(); } \
    o = g_order_table[o];
  #define limit_order_nowarn(o) \
    if (o > g_safe_max_order) o = g_safe_max_order; \
    o = g_order_table[o];
#else
  #define limit_order(o) \
    if (o > g_max_order) warn_order(); \
    o = g_safe_max_order;
  #define limit_order_nowarn(o) \
    o = g_safe_max_order;
#endif

extern HERMES_API void reset_warn_order(); ///< Resets warn order flag.
extern HERMES_API void warn_order(); ///< Warns about integration order iff ward order flags it not set. Sets warn order flag.
extern HERMES_API void update_limit_table(int mode);

#endif

