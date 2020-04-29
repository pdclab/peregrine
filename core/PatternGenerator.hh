#ifndef PATTERN_GENERATOR_HH
#define PATTERN_GENERATOR_HH

#include <vector>
#include "Graph.hh"

namespace Peregrine
{
  using Edges = std::vector<std::pair<uint32_t,uint32_t>>;
  struct PatternGenerator
  {
    PatternGenerator();
    static SmallGraph clique(uint32_t sz);
    static SmallGraph star(uint32_t sz);
    static std::vector<SmallGraph> all(uint32_t sz, bool vertex_based, bool anti_edges);
    static std::vector<SmallGraph> extend(const std::vector<SmallGraph> &from, bool vertex_based, bool overwrite_anti_edges = OVERWRITE_ANTI_EDGES);
    bool is_connected_pattern(Edges edge_list);

    static const bool VERTEX_BASED = true;
    static const bool EDGE_BASED = false;
    static const bool INCLUDE_ANTI_EDGES = true;
    static const bool EXCLUDE_ANTI_EDGES = false;
    static const bool OVERWRITE_ANTI_EDGES = true;
    static const bool MAINTAIN_ANTI_EDGES = false;
  };
}

#endif
