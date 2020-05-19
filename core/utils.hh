#ifndef UTILS_HH
#define UTILS_HH

#include <iostream>
#include <vector>
#include <unordered_map>
#include <sys/types.h>
#include <set>
#include <algorithm>

#include <sys/time.h>
#include <mutex>

namespace utils
{

extern std::mutex logging_mutex;
typedef unsigned long long timestamp_t;
timestamp_t get_timestamp();

struct Log
{
  Log() {}
  ~Log() { std::cout << std::flush; }

#ifdef TESTING
  template <typename T>
  const Log &operator<<(const T &) const
  {
    return *this;
  }
#else
  template <typename T>
  const Log &operator<<(const T &t) const
  {
    std::cout << t;
    return *this;
  }
#endif
};

template <typename T>
void print_alist(const std::unordered_map<T, std::vector<T>> &a_list)
{
    for (auto v : a_list)
    {
        std::cout << v.first << " ";
        for (auto vertex_id : v.second)
        {
            std::cout << vertex_id << " ";
        }
        std::cout << std::endl;
    }
}
template <typename T>
void print_pairs(const std::vector<std::pair<T, T>> &v)
{
    for (auto x : v)
    {
        std::cout << x.first << " " << x.second << std::endl;
    }
}

template <typename T>
void print_vector(const std::vector<T> &v)
{
    for (auto vertex_id : v)
    {
        std::cout << vertex_id << " ";
    }
    std::cout << std::endl;
}

template <typename T>
void print_set(const std::set<T> &v)
{
    for (auto vertex_id : v)
    {
        std::cout << vertex_id << " ";
    }
    std::cout << std::endl;
}

template <typename T>
bool search(const std::vector<T> &vlist, T key){
  return std::find(vlist.begin(), vlist.end(), key) != vlist.end();
}

// much slower than linear search unless vlist is large
// binary search has a misprediction every branch
template <typename T>
bool bsearch(const std::vector<T> &vlist, T key) {
  return std::binary_search(vlist.begin(), vlist.end(), key);
}

} // namespace utils

#endif
