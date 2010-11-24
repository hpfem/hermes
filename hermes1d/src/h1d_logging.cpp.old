#include "h1d_logging.h"

void intro() {
  printf("\n-------------------------------------------\n");
  printf(" This is Hermes1D - a free ODE solver\n");
  printf(" based on the hp-FEM and Newton's method,\n");
  printf(" developed by the hp-FEM group at UNR\n");
  printf(" and distributed under the BSD license.\n");
  printf(" For more details visit http://hpfem.org/.\n");
  printf("-------------------------------------------\n");
}

//// runtime report control varibles //////////////////////////////////////////////////////////////////////////////////
#if defined(H1D_REPORT_WARNING)
# define __H1D_REP_WARN true
#else
# define __H1D_REP_WARN false
#endif
#if defined(H1D_REPORT_INTR_WARNING)
# define __H1D_REP_WARN_INTR true
#else
# define __H1D_REP_WARN_INTR false
#endif
#if defined(H1D_REPORT_INFO)
# define __H1D_REP_INFO true
#else
# define __H1D_REP_INFO false
#endif
#if defined(H1D_REPORT_VERBOSE)
# define __H1D_REP_VERB true
#else
# define __H1D_REP_VERB false
#endif
#if defined(H1D_REPORT_TRACE)
# define __H1D_REP_TRAC true
#else
# define __H1D_REP_TRAC false
#endif
#if defined(H1D_REPORT_TIME)
# define __H1D_REP_TIME true
#else
# define __H1D_REP_TIME false
#endif
#if defined(_DEBUG) || !defined(NDEBUG)
# define __H1D_REP_DEBG true
#else
# define __H1D_REP_DEBG false
#endif

HERMES_API bool __h1d_report_warn = __H1D_REP_WARN;
HERMES_API bool __h1d_report_warn_intr = __H1D_REP_WARN_INTR;
HERMES_API bool __h1d_report_info = __H1D_REP_INFO;
HERMES_API bool __h1d_report_verbose = __H1D_REP_VERB;
HERMES_API bool __h1d_report_trace = __H1D_REP_TRAC;
HERMES_API bool __h1d_report_time = __H1D_REP_TIME;
HERMES_API bool __h1d_report_debug = __H1D_REP_DEBG;

