#include <iostream>
#include <vector>
#include <unordered_map>
#include <sys/types.h>
#include <set>

#include "utils.hh"

namespace utils
{
std::mutex logging_mutex;

timestamp_t get_timestamp() {
  struct timeval now;
  gettimeofday (&now, NULL);
  return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}
} // namespace utils
