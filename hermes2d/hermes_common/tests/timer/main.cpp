#include <string>
#include <sstream>
#include <fstream>

#include "common_time_period.h"

// This test makes sure that timer works properly.

#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

void print_result(bool value) {
  if (value) {
    printf("OK\n");
  }
  else {
    printf("failed\n");
  }
}

bool test_timer_rough() {
  #define SECS            1

  // test #1 running timer fo 3 secs
  printf("* Running timer for %d secs...", SECS);
  fflush(stdout);

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  sleep(SECS);
  cpu_time.tick();

  bool result = (int) cpu_time.accumulated() == SECS;
  print_result(result);
  if (!result) return false;

  return true;
}

bool test_timer_hiprecision() {
  // test #1 running timer fo 123 ms
  printf("* Running timer for 123 ms...");
  fflush(stdout);

  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  usleep(123000);
  cpu_time.tick();

  int v = (int) (cpu_time.accumulated() * 100);
  bool result = v == 12;
  print_result(result);
  if (!result) return false;

  return true;
}

int main(int argc, char* argv[])
{
  if (!test_timer_rough())
    return ERROR_FAILURE;

  if (!test_timer_hiprecision()) 
    return ERROR_FAILURE;

  return ERROR_SUCCESS;
}
