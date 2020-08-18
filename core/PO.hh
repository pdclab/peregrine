#ifndef PO_HH
#define PO_HH

#include <vector>
#include <stdint.h>

#include "./bliss-0.73/graph.hh"

namespace PO {
  std::vector<std::pair<uint32_t, uint32_t>> findPOs(const bliss::Graph &H);
  std::vector<std::vector<uint32_t>> automorphicSets(const bliss::Graph &H);
}

#endif
