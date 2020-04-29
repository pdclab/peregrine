#ifndef BUFFER_HH
#define BUFFER_HH

#include <vector>
#include <string>
#include <unordered_map>
#include <map>
#include <list>
#include <unordered_set>
#include <stdint.h> // uint32_t

#include "Graph.hh"

namespace Peregrine
{

  struct adjlist {
    adjlist() : length(0), ptr(nullptr) {}
    adjlist(uint32_t l, uint32_t *p) : length(l), ptr(p) {}
    uint32_t length;
    uint32_t *ptr;
  };


  class DataGraph {
   public:
    DataGraph(std::string data_graph_path);
    DataGraph(const SmallGraph &p);
    DataGraph(const DataGraph *other);
    DataGraph(DataGraph &&other);
    DataGraph(DataGraph &) = delete;
    ~DataGraph();

    struct page_graph_stats
    {
      uint32_t page_size;
      uint32_t page_count;
      uint32_t vertex_count;
      uint32_t edge_count;
    } page_graph_info;

    void set_rbi(const AnalyzedPattern &rbi);
    void set_known_labels(const std::vector<SmallGraph> &patterns);
    void set_known_labels(const std::vector<uint32_t> &labels);
    bool known_label(const uint32_t label) const;

    const std::vector<uint32_t> &get_upper_bounds(uint32_t v) const;
    const std::vector<uint32_t> &get_lower_bounds(uint32_t v) const;
    const adjlist &get_adj(uint32_t v) const;
    uint32_t get_vgs_count() const;
    const SmallGraph &get_vgs(unsigned fi) const;
    const SmallGraph &get_pattern() const;
    const std::vector<std::vector<uint32_t>> &get_qs(unsigned fi) const;
    uint32_t vmap_at(unsigned fi, uint32_t v, unsigned qsi) const;
    uint32_t label(uint32_t dv) const;
    const std::vector<uint32_t> &get_qo(uint32_t fi) const;

    AnalyzedPattern rbi;
    uint32_t new_label;
   private:
    unsigned forest_count;
    bool labelled_graph = false;
    uint32_t *labels;
    adjlist *data_graph;
    uint32_t *graph_in_memory;
    std::unordered_set<uint32_t> known_labels;
  };
}

#endif
