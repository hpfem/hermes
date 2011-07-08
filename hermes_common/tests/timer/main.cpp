#include "hermes_common.h"
// For the function usleep().
#include <time.h>
#include <unistd.h>

using namespace Hermes;
using namespace std;

// This test makes sure that timer works properly.

void print_result(bool value) 
{
  if (value)
    printf("OK\n");
  else
    printf("failed\n");
}

bool test_timer_hiprecision() 
{
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
  if (!test_timer_hiprecision())
    printf("OK\n");
  else
    printf("failed\n");
}
