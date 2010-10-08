// This file is part of Hermes3D
//
// Copyright (c) 2009 hp-FEM group at the University of Nevada, Reno (UNR).
// Email: hpfem-group@unr.edu, home page: http://hpfem.org/.
//
// Hermes3D is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation; either version 2 of the License,
// or (at your option) any later version.
//
// Hermes3D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes3D; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "common.h"

using namespace std;

/// A size of a delimiter in a log file. \internal \ingroup g_logging
#define H3D_LOG_FILE_DELIM_SIZE 80

void hermes3d_exit_if(bool cond, int code) {
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
/** \param[in] code An event code, e.g., ::h3d_EC_ERROR.
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
    case H3D_EC_ERROR:
    case H3D_EC_ASSERT: console_attrs |= console_attr_red; break;
    case H3D_EC_WARNING: console_attrs |= console_attr_red | console_attr_green; break;
    case H3D_EC_INFO: console_bold = true;
    case H3D_EC_VERBOSE: console_attrs |= console_attr_red | console_attr_green | console_attr_blue; break;
    case H3D_EC_TRACE: console_attrs |= console_attr_blue; break;
    case H3D_EC_TIME: console_attrs |= console_attr_green | console_attr_blue; break;
    case H3D_EC_DEBUG: console_attrs |= console_attr_red | console_attr_blue; break;
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
    case H3D_EC_ERROR:
    case H3D_EC_ASSERT: console_attrs |= FOREGROUND_RED; break;
    case H3D_EC_WARNING: console_attrs |= FOREGROUND_RED | FOREGROUND_GREEN; break;
    case H3D_EC_INFO: console_bold = true;
    case H3D_EC_VERBOSE: console_attrs |= FOREGROUND_BLUE; break;
    case H3D_EC_TRACE: console_attrs |= FOREGROUND_BLUE; break;
    case H3D_EC_TIME: console_attrs |= FOREGROUND_GREEN | FOREGROUND_BLUE; break;
    case H3D_EC_DEBUG: console_attrs |= FOREGROUND_RED | FOREGROUND_BLUE; break;
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

bool hermes3d_log_message_if(bool cond, const Hermes3DLogEventInfo& info, const char* msg, ...) {
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
          for(int i = 0; i < H3D_LOG_FILE_DELIM_SIZE; i++)
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