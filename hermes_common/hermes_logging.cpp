#include "hermes_logging.h"

//// logo //////////////////////////////////////////////////////////////////////////////////
#ifndef HERMES_NO_LOGO
/// Generates a logo when the library is loaded and the logo is enabled (a preprocessor directive ::HERMES_NO_LOGO). \internal \ingroup g_logging
class HermesLogoMessage {
public:
  HermesLogoMessage() {
    printf("\n-------------------------------------------------\n");
    printf("         This application uses Hermes.\n");
    printf("       Hermes is a C++ library for rapid \n");
    printf("  development of adaptive FEM and hp-FEM solvers\n");
    printf("      developed by the hp-FEM group at UNR\n");
    printf("     and distributed under the GPL license.\n");
    printf("    For more details visit http://hpfem.org/.\n");
    printf("-------------------------------------------------\n");
    fflush(stdout);
  }
};
// FIXME: Currently, this can't be disabled from cmake (due to some bugs), so
// cmake has to be fixed first, then this can be enabled again.
//HermesLogoMessage hermes_logo_message;
#endif

//// runtime report control varibles //////////////////////////////////////////////////////////////////////////////////
#if defined(HERMES_REPORT_WARNING)
# define __HERMES_REP_WARN true
#else
# define __HERMES_REP_WARN false
#endif
#if defined(HERMES_REPORT_INTR_WARNING)
# define __HERMES_REP_WARN_INTR true
#else
# define __HERMES_REP_WARN_INTR false
#endif
#if defined(HERMES_REPORT_INFO)
# define __HERMES_REP_INFO true
#else
# define __HERMES_REP_INFO false
#endif
#if defined(HERMES_REPORT_VERBOSE)
# define __HERMES_REP_VERB true
#else
# define __HERMES_REP_VERB false
#endif
#if defined(HERMES_REPORT_TRACE)
# define __HERMES_REP_TRAC true
#else
# define __HERMES_REP_TRAC false
#endif
#if defined(HERMES_REPORT_TIME)
# define __HERMES_REP_TIME true
#else
# define __HERMES_REP_TIME false
#endif
#if defined(_DEBUG) || !defined(NDEBUG)
# define __HERMES_REP_DEBG true
#else
# define __HERMES_REP_DEBG false
#endif

HERMES_API bool __hermes_report_warn = __HERMES_REP_WARN;
HERMES_API bool __hermes_report_warn_intr = __HERMES_REP_WARN_INTR;
HERMES_API bool __hermes_report_info = __HERMES_REP_INFO;
HERMES_API bool __hermes_report_verbose = __HERMES_REP_VERB;
HERMES_API bool __hermes_report_trace = __HERMES_REP_TRAC;
HERMES_API bool __hermes_report_time = __HERMES_REP_TIME;
HERMES_API bool __hermes_report_debug = __HERMES_REP_DEBG;

