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

#include "common.h"

using namespace std;

/// A size of a delimiter in a log file. \internal \ingroup g_logging
#define H2D_LOG_FILE_DELIM_SIZE 80

H2D_API const std::string get_quad_order_str(const int quad_order) {
  std::stringstream str;
  str << "(H:" << H2D_GET_H_ORDER(quad_order) << ";V:" << H2D_GET_V_ORDER(quad_order) << ")";
  return str.str();
}

H2D_API int make_edge_order(int mode, int edge, int encoded_order)
{
  assert(edge < 4);
  
  if (mode == H2D_MODE_TRIANGLE || edge == 0 || edge == 2)
    return H2D_GET_H_ORDER(encoded_order);
  else
    return H2D_GET_V_ORDER(encoded_order);
}

H2D_API void hermes2d_exit_if(bool cond, int code) {
  if (cond)
    exit(code);
}

/// Logging output monitor. \internal \ingroup g_logging
/** This class protects a logging function __h2d_log_message_if() in multithreded environment. */
class LoggerMonitor {
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

static map<string, bool> logger_written; ///< A list of all log files that were used to write a log. Used to write a log header to a log file. \internal \ingroup g_logging

/// Writes a fancy formatted text to a console. \internal \ingroup g_logging
/** \param[in] code An event code, e.g., ::H2D_EC_ERROR.
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
    case H2D_EC_ERROR:
    case H2D_EC_ASSERT: console_attrs |= console_attr_red; break;
    case H2D_EC_WARNING: console_attrs |= console_attr_red | console_attr_green; break;
    case H2D_EC_INFO: console_bold = true;
    case H2D_EC_VERBOSE: console_attrs |= console_attr_red | console_attr_green | console_attr_blue; break;
    case H2D_EC_TRACE: console_attrs |= console_attr_blue; break;
    case H2D_EC_TIME: console_attrs |= console_attr_green | console_attr_blue; break;
    case H2D_EC_DEBUG: console_attrs |= console_attr_red | console_attr_blue; break;
    default: error("Unknown error code: '%c'", code); break;
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
    case H2D_EC_ERROR:
    case H2D_EC_ASSERT: console_attrs |= FOREGROUND_RED; break;
    case H2D_EC_WARNING: console_attrs |= FOREGROUND_RED | FOREGROUND_GREEN; break;
    case H2D_EC_INFO: console_bold = true;
    case H2D_EC_VERBOSE: console_attrs |= FOREGROUND_BLUE; break;
    case H2D_EC_TRACE: console_attrs |= FOREGROUND_BLUE; break;
    case H2D_EC_TIME: console_attrs |= FOREGROUND_GREEN | FOREGROUND_BLUE; break;
    case H2D_EC_DEBUG: console_attrs |= FOREGROUND_RED | FOREGROUND_BLUE; break;
    default: error("Unknown error code: '%c'", code); break;
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

H2D_API bool hermes2d_log_message_if(bool cond, const Hermes2DLogEventInfo& info, const char* msg, ...) {
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
        map<string, bool>::const_iterator found = logger_written.find(info.log_file);
        if (found == logger_written.end()) { //first write, write delimited to a file
          logger_written[info.log_file] = true;
          fprintf(file, "\n");
          for(int i = 0; i < H2D_LOG_FILE_DELIM_SIZE; i++)
            fprintf(file, "-");
          fprintf(file, "\n\n");
        }

        //build a long version of location
        ostringstream location;
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

void __hermes2d_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream, const Hermes2DLogEventInfo& err_info)
{
  if (fwrite(ptr, size, nitems, stream) != nitems || ferror(stream))
    hermes2d_exit_if(hermes2d_log_message_if(true, err_info, "Error writing to file: %s", strerror(ferror(stream))));
}

void __hermes2d_fread(void* ptr, size_t size, size_t nitems, FILE* stream, const Hermes2DLogEventInfo& err_info)
{
  size_t ret = fread(ptr, size, nitems, stream);
  if (ret < nitems)
    hermes2d_exit_if(hermes2d_log_message_if(true, err_info, "Premature end of file."));
  else if (ferror(stream))
    hermes2d_exit_if(hermes2d_log_message_if(true, err_info, "Error reading file: %s", strerror(ferror(stream))));
}

//// logo //////////////////////////////////////////////////////////////////////////////////
#ifndef H2D_NO_LOGO
/// Generates a logo when the library is loaded and the logo is enabled (a preprocessor directive ::H2D_NO_LOGO). \internal \ingroup g_logging
class Hermes2DLogoMessage {
public:
  Hermes2DLogoMessage() {
    printf("\n-------------------------------------------------\n");
    printf("          This application uses Hermes2D\n");
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

H2D_API bool __h2d_report_warn = __H2D_REP_WARN;
H2D_API bool __h2d_report_warn_intr = __H2D_REP_WARN_INTR;
H2D_API bool __h2d_report_info = __H2D_REP_INFO;
H2D_API bool __h2d_report_verbose = __H2D_REP_VERB;
H2D_API bool __h2d_report_trace = __H2D_REP_TRAC;
H2D_API bool __h2d_report_time = __H2D_REP_TIME;
H2D_API bool __h2d_report_debug = __H2D_REP_DEBG;

//// python support //////////////////////////////////////////////////////////////////////////////////
H2D_API void throw_exception(char *text)
{
  throw std::runtime_error(text);
}

