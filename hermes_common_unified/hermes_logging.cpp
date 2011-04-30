#include "hermes_logging.h"
#include "third_party_codes/trilinos-teuchos/Teuchos_stacktrace.hpp"
#ifdef WIN32 //Windows platform
	#include <Windows.h>
#endif

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



void hermes_exit_if(bool cond, int code) {
  if (cond)
    exit(code);
}

/// Logging output monitor. \internal \ingroup g_logging
/** This class protects a logging function __hermes_log_message_if() in multithreded environment. */
class LoggerMonitor 
{
  pthread_mutexattr_t mutex_attr; ///< Mutext attributes.
  pthread_mutex_t mutex; ///< Mutex that protects monitor.

public:
  /// Constructor. Creates a mutex.
  LoggerMonitor() {
    pthread_mutexattr_init(&mutex_attr);
    pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_RECURSIVE);
    pthread_mutex_init(&mutex, &mutex_attr);
  };
  /// Destructor. Deletes a mutex.
  ~LoggerMonitor() {
    pthread_mutex_destroy(&mutex);
    pthread_mutexattr_destroy(&mutex_attr);
  };

  /// Enters protected section.
  void enter() { pthread_mutex_lock(&mutex); };

  /// Leaves protected section.
  void leave() { pthread_mutex_unlock(&mutex); };
};

static LoggerMonitor logger_monitor; ///< A monitor that protects logging function. \internal \ingroup g_logging

static std::map<std::string, bool> logger_written; ///< A list of all log files that were used to write a log. Used to write a log header to a log file. \internal \ingroup g_logging

/// Writes a fancy formatted text to a console. \internal \ingroup g_logging
/** \param[in] code An event code, e.g., ::HERMES_EC_ERROR.
 *  \param[in] emphasize True if the message should be emphasized.
 *  \param[in] text A message. A C-style string.
 *  \return True if the message was written. False if it failed due to some reasone. */
static bool write_console(const char code, const bool emphasize, const char* text) {
#ifdef WIN32 //Windows platform
  HANDLE h_console = GetStdHandle(STD_OUTPUT_HANDLE);
  if (h_console == INVALID_HANDLE_VALUE)
    return false;

  //read current console settings
  CONSOLE_SCREEN_BUFFER_INFO console_info;
  if (!GetConsoleScreenBufferInfo(h_console, &console_info))
    return false;

  //generate console settings
  WORD console_attr_red = FOREGROUND_RED, console_attr_green = FOREGROUND_GREEN, console_attr_blue = FOREGROUND_BLUE;
  if (emphasize) { //invert foreground and background
    console_attr_red = BACKGROUND_RED;
    console_attr_green = BACKGROUND_GREEN;
    console_attr_blue = BACKGROUND_BLUE;
  }
  WORD console_attrs = 0;
  bool console_bold = false;
  switch(code) {
    case HERMES_EC_ERROR:
    case HERMES_EC_ASSERT: console_attrs |= console_attr_red; break;
    case HERMES_EC_WARNING: console_attrs |= console_attr_red | console_attr_green; break;
    case HERMES_EC_INFO: console_bold = true;
    case HERMES_EC_VERBOSE: console_attrs |= console_attr_red | console_attr_green | console_attr_blue; break;
    case HERMES_EC_TRACE: console_attrs |= console_attr_blue; break;
    case HERMES_EC_TIME: console_attrs |= console_attr_green | console_attr_blue; break;
    case HERMES_EC_DEBUG: console_attrs |= console_attr_red | console_attr_blue; break;
    default: printf("Unknown error code: '%c'", code); exit(-1);
  }
  if (console_bold && !emphasize)
    console_attrs |= FOREGROUND_INTENSITY;

  //set new console settings
  SetConsoleTextAttribute(h_console, console_attrs);

  //write text
  DWORD num_written;
  BOOL write_success = WriteConsoleA(h_console, text, strlen(text), &num_written, NULL);

  //return previous settings
  SetConsoleTextAttribute(h_console, console_info.wAttributes);

  if (write_success)
    return true;
  else
    return false;
#else //Linux platform
# define FOREGROUND_RED 1
# define FOREGROUND_GREEN 2
# define FOREGROUND_BLUE 4
  //console color code
  int console_attrs = 0;
  bool console_bold = false;
  switch(code) {
    case HERMES_EC_ERROR:
    case HERMES_EC_ASSERT: console_attrs |= FOREGROUND_RED; break;
    case HERMES_EC_WARNING: console_attrs |= FOREGROUND_RED | FOREGROUND_GREEN; break;
    case HERMES_EC_INFO: console_bold = true;
    case HERMES_EC_VERBOSE: console_attrs |= FOREGROUND_BLUE; break;
    case HERMES_EC_TRACE: console_attrs |= FOREGROUND_BLUE; break;
    case HERMES_EC_TIME: console_attrs |= FOREGROUND_GREEN | FOREGROUND_BLUE; break;
    case HERMES_EC_DEBUG: console_attrs |= FOREGROUND_RED | FOREGROUND_BLUE; break;
    default: printf("Unknown error code: '%c'", code); exit(-1);
  }
  printf("\033[%dm", console_attrs + 30);

  //emphasize and console bold
  if (emphasize)
    printf("\033[7m");
  else if (console_bold)
    printf("\033[1m");

  //print text and reset settings
  printf("%s\033[0m", text);

  return true;
#endif
}

bool hermes_log_message_if(bool cond, const HermesLogEventInfo& info, const char* msg, ...) {
  if (cond) {
    logger_monitor.enter();

    //print message to a buffer (since vfprintf modifies arglist such that it becomes unusable)
    //not safe, but C does not offer any other multiplatform solution. Since vsnprintf modifies the var-args, it cannot be called repeatedly.
    #define BUF_SZ 2048
    bool emphasize = false;
    bool new_block = true;
    char text[BUF_SZ];
    char* text_contents = text + 1;
    if (msg[0] == '!') {
      emphasize = true;
      msg++;
    }
    if (msg[0] == ' ') {
      text[0] = ' ';
      new_block = false;
    }
    else {
      if (info.code == 'E')
          Teuchos::show_stacktrace();
      text[0] = info.code;
      text[1] = ' ';
      text_contents++;
      new_block = true;
    }

    //print the message
    va_list arglist;
    va_start(arglist, msg);
    vsprintf(text_contents, msg, arglist);
    va_end(arglist);

    //print the message
    if (emphasize && new_block)
      printf("\n");
    if (!write_console(info.code, emphasize, text))
      printf("%s", text); //safe fallback
    printf("\n"); //write a new line

    //print to file
    if (info.log_file != NULL) {
      FILE* file = fopen(info.log_file, "at");
      if (file != NULL)
      {
        //check whether log file was already written
        std::map<std::string, bool>::const_iterator found = logger_written.find(info.log_file);
        if (found == logger_written.end()) { //first write, write delimited to a file
          logger_written[info.log_file] = true;
          fprintf(file, "\n");
          for(int i = 0; i < HERMES_LOG_FILE_DELIM_SIZE; i++)
            fprintf(file, "-");
          fprintf(file, "\n\n");
        }

        //build a long version of location
        std::ostringstream location;
        location << '(';
        if (info.src_function != NULL) {
          location << info.src_function;
          if (info.src_file != NULL)
            location << '@';
        }
        if (info.src_file != NULL)
          location << info.src_file << ':' << info.src_line;
        location << ')';

        //get time
        time_t now;
        time(&now);
        struct tm* now_tm = gmtime(&now);
        char time_buf[BUF_SZ];
        strftime(time_buf, BUF_SZ, "%y%m%d-%H:%M", now_tm);

        //write
        if (emphasize && new_block)
          fprintf(file, "\n\n");
        fprintf(file, "%s\t%s %s\n", time_buf, text, location.str().c_str());
        fclose(file);
      }
    }

    logger_monitor.leave();
  }
  return cond;
}

void __hermes_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream, const HermesLogEventInfo& err_info)
{
  if (fwrite(ptr, size, nitems, stream) != nitems || ferror(stream))
    hermes_exit_if(hermes_log_message_if(true, err_info, "Error writing to file: %s", strerror(ferror(stream))));
}

void __hermes_fread(void* ptr, size_t size, size_t nitems, FILE* stream, const HermesLogEventInfo& err_info)
{
  size_t ret = fread(ptr, size, nitems, stream);
  if (ret < nitems)
    hermes_exit_if(hermes_log_message_if(true, err_info, "Premature end of file."));
  else if (ferror(stream))
    hermes_exit_if(hermes_log_message_if(true, err_info, "Error reading file: %s", strerror(ferror(stream))));
}

HermesLogEventInfo::HermesLogEventInfo(const char code, const char* log_file, const char* src_function, const char* src_file, const int src_line)
    : code(code), log_file(log_file), src_function(src_function), src_file(src_file), src_line(src_line) 
{}
