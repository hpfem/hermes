#include "h2d_logging.h"

//// logo //////////////////////////////////////////////////////////////////////////////////
#ifndef H2D_NO_LOGO
/// Generates a logo when the library is loaded and the logo is enabled (a preprocessor directive ::H2D_NO_LOGO). \internal \ingroup g_logging
class Hermes2DLogoMessage {
public:
  Hermes2DLogoMessage() {
    printf("\n-------------------------------------------------\n");
    printf("         This application uses Hermes2D.\n");
    printf("       Hermes2D is a C++ library for rapid \n");
    printf("  development of adaptive FEM and hp-FEM solvers\n");
    printf("      developed by the hp-FEM group at UNR\n");
    printf("     and distributed under the GPL license.\n");
    printf("    For more details visit http://hpfem.org/.\n");
    printf("-------------------------------------------------\n");
    fflush(stdout);
  }
};
Hermes2DLogoMessage hermes2d_logo_message;
#endif

//// runtime report control varibles //////////////////////////////////////////////////////////////////////////////////
#if defined(H2D_REPORT_WARNING)
# define __H2D_REP_WARN true
#else
# define __H2D_REP_WARN false
#endif
#if defined(H2D_REPORT_INTR_WARNING)
# define __H2D_REP_WARN_INTR true
#else
# define __H2D_REP_WARN_INTR false
#endif
#if defined(H2D_REPORT_INFO)
# define __H2D_REP_INFO true
#else
# define __H2D_REP_INFO false
#endif
#if defined(H2D_REPORT_VERBOSE)
# define __H2D_REP_VERB true
#else
# define __H2D_REP_VERB false
#endif
#if defined(H2D_REPORT_TRACE)
# define __H2D_REP_TRAC true
#else
# define __H2D_REP_TRAC false
#endif
#if defined(H2D_REPORT_TIME)
# define __H2D_REP_TIME true
#else
# define __H2D_REP_TIME false
#endif
#if defined(_DEBUG) || !defined(NDEBUG)
# define __H2D_REP_DEBG true
#else
# define __H2D_REP_DEBG false
#endif

HERMES_API bool __h2d_report_warn = __H2D_REP_WARN;
HERMES_API bool __h2d_report_warn_intr = __H2D_REP_WARN_INTR;
HERMES_API bool __h2d_report_info = __H2D_REP_INFO;
HERMES_API bool __h2d_report_verbose = __H2D_REP_VERB;
HERMES_API bool __h2d_report_trace = __H2D_REP_TRAC;
HERMES_API bool __h2d_report_time = __H2D_REP_TIME;
HERMES_API bool __h2d_report_debug = __H2D_REP_DEBG;

