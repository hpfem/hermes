// This file is part of HermesCommon
//
// Hermes is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// Hermes is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes; if not, see <http://www.gnu.prg/licenses/>.
#include "mixins.h"
#include <map>
#include <string>
#include "common.h"

namespace Hermes
{
  namespace Mixins
  {
    Loggable::LoggerMonitor Loggable::logger_monitor;

    std::map<std::string, bool> Loggable::logger_written;

    Loggable::Loggable(bool verbose_output, callbackFn verbose_callback) : verbose_output(verbose_output), verbose_callback(verbose_callback)
    {
    }

    bool Loggable::get_verbose_output() const
    {
      return this->verbose_output;
    }

    Loggable::callbackFn Loggable::get_verbose_callback() const
    {
      return this->verbose_callback;
    }

    void Loggable::Static::info(const char* msg, ...)
    {
      char text[BUF_SZ];
      char* text_contents = text + 1;

      text[0] = HERMES_EC_INFO;
      text[1] = ' ';
      text_contents++;

      //print the message
      va_list arglist;
      va_start(arglist, msg);
      vsprintf(text_contents, msg, arglist);
      va_end(arglist);

      //Windows platform
    #ifdef WIN32
      HANDLE h_console = GetStdHandle(STD_OUTPUT_HANDLE);

      //read current console settings
      CONSOLE_SCREEN_BUFFER_INFO console_info;

      //generate console settings
      WORD console_attr_red = FOREGROUND_RED, console_attr_green = FOREGROUND_GREEN, console_attr_blue = FOREGROUND_BLUE;

      WORD console_attrs = 0;
      console_attrs |= console_attr_green;

      //set new console settings
      SetConsoleTextAttribute(h_console, console_attrs);

      //write text
      DWORD num_written;
      BOOL write_success = WriteConsoleA(h_console, text_contents, strlen(text_contents), &num_written, NULL);
      std::cout << std::endl;

    //Linux platform
    #else
    # define FOREGROUND_RED 1
    # define FOREGROUND_GREEN 2
    # define FOREGROUND_BLUE 4
      //console color code
      int console_attrs = 0;
      bool console_bold = true;

      printf("\033[%dm", console_attrs + 30);

      //emphasize and console bold
      if(console_bold)
        printf("\033[1m");

      //print text and reset settings
      printf("%s\033[0m\n", text_contents);

    #endif
    }

    void Loggable::Static::warn(const char* msg, ...)
    {
      char text[BUF_SZ];
      char* text_contents = text + 1;

      text[0] = HERMES_EC_WARNING;
      text[1] = ' ';
      text_contents++;

      //print the message
      va_list arglist;
      va_start(arglist, msg);
      vsprintf(text_contents, msg, arglist);
      va_end(arglist);

      //Windows platform
    #ifdef WIN32
      HANDLE h_console = GetStdHandle(STD_OUTPUT_HANDLE);

      //read current console settings
      CONSOLE_SCREEN_BUFFER_INFO console_info;

      //generate console settings
      WORD console_attr_red = FOREGROUND_RED, console_attr_green = FOREGROUND_GREEN, console_attr_blue = FOREGROUND_BLUE;

      WORD console_attrs = 0;
      console_attrs |= console_attr_red | console_attr_green;

      //set new console settings
      SetConsoleTextAttribute(h_console, console_attrs);

      //write text
      DWORD num_written;
      BOOL write_success = WriteConsoleA(h_console, text_contents, strlen(text_contents), &num_written, NULL);
      std::cout << std::endl;
    //Linux platform
    #else
    # define FOREGROUND_RED 1
    # define FOREGROUND_GREEN 2
    # define FOREGROUND_BLUE 4
      //console color code
      int console_attrs = 1;
      console_attrs |= FOREGROUND_RED | FOREGROUND_GREEN;
      bool console_bold = false;

      printf("\033[%dm", console_attrs + 30);

      //emphasize and console bold
      if(console_bold)
        printf("\033[1m");

      //print text and reset settings
      printf("%s\033[0m\n", text_contents);

    #endif
    }

    void Loggable::warn(const char* msg, ...) const
    {
      if(!this->verbose_output)
        return;

      char text[BUF_SZ];
      char* text_contents = text + 1;

      text[0] = HERMES_EC_WARNING;
      text[1] = ' ';
      text_contents++;

      //print the message
      va_list arglist;
      va_start(arglist, msg);
      vsprintf(text_contents, msg, arglist);
      va_end(arglist);
      hermes_log_message(HERMES_EC_WARNING, text_contents);
    }

    void Loggable::warn_if(bool cond, const char* msg, ...) const
    {
      if(!this->verbose_output)
        return;

      if(cond)
      {
        char text[BUF_SZ];
        char* text_contents = text + 1;

        text[0] = HERMES_EC_WARNING;
        text[1] = ' ';
        text_contents++;

        //print the message
        va_list arglist;
        va_start(arglist, msg);
        vsprintf(text_contents, msg, arglist);
        va_end(arglist);
        hermes_log_message(HERMES_EC_WARNING, text_contents);
      }
    }
    void Loggable::info(const char* msg, ...) const
    {
      if(!this->verbose_output)
        return;

      char text[BUF_SZ];
      char* text_contents = text + 1;

      text[0] = HERMES_EC_INFO;
      text[1] = ' ';
      text_contents++;

      //print the message
      va_list arglist;
      va_start(arglist, msg);
      vsprintf(text_contents, msg, arglist);
      va_end(arglist);
      hermes_log_message(HERMES_EC_INFO, text_contents);
    }
    void Loggable::info_if(bool cond, const char* msg, ...) const
    {
      if(!this->verbose_output)
        return;

      if(cond)
      {
        char text[BUF_SZ];
        char* text_contents = text + 1;

        text[0] = HERMES_EC_INFO;
        text[1] = ' ';
        text_contents++;

        //print the message
        va_list arglist;
        va_start(arglist, msg);
        vsprintf(text_contents, msg, arglist);
        va_end(arglist);
        hermes_log_message(HERMES_EC_INFO, text_contents);
      }
    }

    bool Loggable::write_console(const char code, const char* text) const
    {
      //Windows platform
    #ifdef WIN32
      HANDLE h_console = GetStdHandle(STD_OUTPUT_HANDLE);
      if(h_console == INVALID_HANDLE_VALUE)
        return false;

      //read current console settings
      CONSOLE_SCREEN_BUFFER_INFO console_info;
      if(!GetConsoleScreenBufferInfo(h_console, &console_info))
        return false;

      //generate console settings
      WORD console_attr_red = FOREGROUND_RED, console_attr_green = FOREGROUND_GREEN, console_attr_blue = FOREGROUND_BLUE;

      WORD console_attrs = 0;
      switch(code)
      {
        case HERMES_EC_WARNING: console_attrs |= console_attr_red | console_attr_green; break;
        case HERMES_EC_INFO:console_attrs |= console_attr_green;  break;
        default: throw Hermes::Exceptions::Exception("Unknown error code: '%c'", code);
      }

      //set new console settings
      SetConsoleTextAttribute(h_console, console_attrs);

      //write text
      DWORD num_written;
      BOOL write_success = WriteConsoleA(h_console, text, strlen(text), &num_written, NULL);

      //return previous settings
      SetConsoleTextAttribute(h_console, console_info.wAttributes);

      if(write_success)
        return true;
      else
        return false;
    //Linux platform
    #else
    # define FOREGROUND_RED 1
    # define FOREGROUND_GREEN 2
    # define FOREGROUND_BLUE 4
      //console color code
      int console_attrs = 0;
      bool console_bold = false;
      switch(code)
      {
      case HERMES_EC_WARNING: console_attrs |= FOREGROUND_RED | FOREGROUND_GREEN; break;
      case HERMES_EC_INFO: console_bold = true; break;
      default: throw Hermes::Exceptions::Exception("Unknown error code: '%c'", code);
      }

      printf("\033[%dm", console_attrs + 30);

      //emphasize and console bold
      if(console_bold)
        printf("\033[1m");

      //print text and reset settings
      printf("%s\033[0m", text);

      return true;
    #endif
    }

    Loggable::HermesLogEventInfo* Loggable::hermes_build_log_info(char event) const
    {
        return new Loggable::HermesLogEventInfo(event, HERMES_LOG_FILE, __CURRENT_FUNCTION, __FILE__, __LINE__);
    }

    void Loggable::hermes_fwrite(const void* ptr, size_t size, size_t nitems, FILE* stream) const
    {
      if(fwrite(ptr, size, nitems, stream) != nitems || ferror(stream))
        throw Hermes::Exceptions::Exception("Error writing to file: %s", strerror(ferror(stream)));
    }

    void Loggable::hermes_fread(void* ptr, size_t size, size_t nitems, FILE* stream) const
    {
      size_t ret = fread(ptr, size, nitems, stream);
      if(ret < nitems)
        throw Hermes::Exceptions::Exception("Premature end of file.");
      else if(ferror(stream))
        throw Hermes::Exceptions::Exception("Error reading file: %s", strerror(ferror(stream)));
    }

    Loggable::HermesLogEventInfo::HermesLogEventInfo(const char code, const char* log_file, const char* src_function, const char* src_file, const int src_line)
      : code(code), log_file(log_file), src_function(src_function), src_file(src_file), src_line(src_line)
    {}

    void Loggable::hermes_log_message(const char code, const char* msg) const
    {
      logger_monitor.enter();

      //print the message
      if(!write_console(code, msg))
        printf("%s", msg);  //safe fallback
      printf("\n");  //write a new line

      HermesLogEventInfo* info = this->hermes_build_log_info(code);

      //print to file
      if(info->log_file != NULL)
      {
        FILE* file = fopen(info->log_file, "at");
        if(file != NULL)
        {
          //check whether log file was already written
          std::map<std::string, bool>::const_iterator found = logger_written.find(info->log_file);
          if(found == logger_written.end()) {  //first write, write delimited to a file
            logger_written[info->log_file] = true;
            fprintf(file, "\n");
            for(int i = 0; i < HERMES_LOG_FILE_DELIM_SIZE; i++)
              fprintf(file, "-");
            fprintf(file, "\n\n");
          }

          //build a long version of location
          std::ostringstream location;
          location << '(';
          if(info->src_function != NULL)
          {
            location << info->src_function;
            if(info->src_file != NULL)
              location << '@';
          }
          if(info->src_file != NULL)
            location << info->src_file << ':' << info->src_line;
          location << ')';

          //get time
          time_t now;
          time(&now);
          struct tm* now_tm = gmtime(&now);
          char time_buf[BUF_SZ];
          strftime(time_buf, BUF_SZ, "%y%m%d-%H:%M", now_tm);

          //write
          fprintf(file, "%s\t%s %s\n", time_buf, msg, location.str().c_str());
          fclose(file);

          if(this->verbose_callback != NULL)
            this->verbose_callback(msg);
        }
      }

      delete info;

      logger_monitor.leave();
    }

    void Loggable::set_verbose_output(bool to_set)
    {
      this->verbose_output = to_set;
    }

    void Loggable::set_verbose_callback(callbackFn callback)
    {
      this->verbose_callback = callback;
    }

    TimeMeasurable::TimeMeasurable(const char *name) : period_name(name == NULL ? "unnamed" : name)
    {
      //initialization
  #ifdef WIN32  //Windows
      LARGE_INTEGER freq;
      if(QueryPerformanceFrequency(&freq))
        frequency = (double)freq.QuadPart;
      else
        frequency = -1;
  #endif  //Linux
      tick_reset();
    }

    TimeMeasurable::SysTime TimeMeasurable::get_time() const
    {
  #ifdef WIN32  //Windows
      if(frequency > 0)
      {
        LARGE_INTEGER ticks;
        QueryPerformanceCounter(&ticks);
        return ticks.QuadPart;
      }
      else
      {
        return clock();
      }
  #elif defined(__APPLE__)  //Mac
      // FIXME: implement time measurement on Mac
      timespec tm;
      return tm;
  #else  //Linux
      timespec tm;
      clock_gettime(CLOCK_REALTIME, &tm);
      return tm;
  #endif
    }

    double TimeMeasurable::period_in_seconds(const SysTime& begin, const SysTime& end) const
    {
  #ifdef WIN32  //Windows
      uint64_t period = end - begin;
      if(frequency > 0)
        return period / frequency;
      else
        return period / (double)CLOCKS_PER_SEC;
  #else  //Linux
      int sec_corr = 0;
      long period_nsec = end.tv_nsec - begin.tv_nsec;
      if(period_nsec < 0)
      {
        sec_corr += -1;
        period_nsec += 1000000000UL;
      }
      long period_sec = (long)(end.tv_sec - begin.tv_sec) + sec_corr;
      return period_sec + (period_nsec / 1E9);
  #endif
    }

    const TimeMeasurable& TimeMeasurable::tick(TimerPeriodTickType type)
    {
      SysTime cur_time = get_time();
      if(type == HERMES_ACCUMULATE)
      {
        double secs = period_in_seconds(last_time, cur_time);
        accum += secs;
        last_period = secs;
      }
      else
        last_period = 0.0;

      last_time = cur_time;
      return *this;
    }

    const std::string& TimeMeasurable::name() const
    {
      return period_name;
    }

    double TimeMeasurable::accumulated() const
    {
      return accum;
    }

    std::string TimeMeasurable::accumulated_str() const {
      return to_string(accum);
    }

    double TimeMeasurable::last() const
    {
      return last_period;
    }

    std::string TimeMeasurable::last_str() const
    {
      return to_string(last_period);
    }

    const TimeMeasurable& TimeMeasurable::tick_reset()
    {
      tick(HERMES_SKIP);
      reset();
      return *this;
    }

    const TimeMeasurable& TimeMeasurable::reset()
    {
      accum = 0;
      last_time = get_time();
      last_period = 0.0;
      return *this;
    }

    std::string TimeMeasurable::to_string(double secs) const
    {
      if(secs < 0)
        return "NO TIME";
      else
      {
        int hours = (int) secs / (3600);
        int mins = (int) fmod(secs, 3600) / 60;
        secs = fmod(secs, 60);

        std::stringstream str;
        if(hours > 0)
          str << hours << "h ";
        if(hours > 0 || mins > 0)
          str << mins << "m ";
        str << secs << "s";

        return str.str();
      }
    }
  }
}