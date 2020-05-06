#include <iostream>

#include <string>
#include <cstring>
#include <algorithm>
#include <unordered_set>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "utils.hh"
#include "Graph.hh"
#include "DataGraph.hh"

namespace Peregrine
{
  DataGraph::DataGraph(const SmallGraph &pp)
  {
    SmallGraph p(pp);
    graph_in_memory = new uint32_t[2 * p.num_true_edges()];
    data_graph = new adjlist[p.num_vertices()+1];

    uint32_t cursor = 0;
    for (uint32_t v = 1; v <= p.num_vertices(); ++v) {
      std::sort(p.true_adj_list.at(v).begin(), p.true_adj_list.at(v).end());

      std::memcpy(&graph_in_memory[cursor], &p.true_adj_list.at(v)[0], p.true_adj_list.at(v).size() * sizeof(uint32_t));
      data_graph[v-1].ptr = &graph_in_memory[cursor];
      data_graph[v-1].length = p.true_adj_list.at(v).size();
      cursor += p.true_adj_list.at(v).size();
    }

    vertex_count = p.num_vertices();
    edge_count = p.num_true_edges();

    labels = new uint32_t[p.num_vertices()+1]();
    if (p.labelling == Graph::LABELLED)
    {
      labelled_graph = true;
      for (uint32_t u = 0; u < p.labels.size(); ++u)
      {
        labels[u+1] = p.labels[u];
      }
    }
  }

  DataGraph::DataGraph(std::string data_graph_path)
  {
    struct stat st;
    std::string data_path(data_graph_path + "/data.bin");
    if (stat(data_path.c_str(), &st) != 0)
    {
      std::cerr << "ERROR: Data graph could not be opened." << std::endl;
      exit(1);
    }


    std::ifstream input_graph(data_path.c_str(), std::ios::binary);
    input_graph.read(reinterpret_cast<char *>(&vertex_count), sizeof(vertex_count));
    input_graph.read(reinterpret_cast<char *>(&edge_count), sizeof(edge_count));

    // don't count the header which was 3 32-bit integers
    assert(st.st_size % 4 == 0);
    uint64_t data_count = st.st_size / 4;
    graph_in_memory = new uint32_t[data_count - 3];

    uint64_t curr_read_offset = 0;
    uint32_t read_batch_size = 2147479552;
    while (input_graph.read(reinterpret_cast<char *>(&graph_in_memory[curr_read_offset/4]), read_batch_size))
    {
      curr_read_offset += read_batch_size;
    }

    data_graph = new adjlist[vertex_count];

    uint32_t cursor = 0;
    for (uint32_t i = 0; i < vertex_count; i++)
    {
      data_graph[i].length = graph_in_memory[cursor];
      data_graph[i].ptr = &graph_in_memory[++cursor];
      cursor += data_graph[i].length;
    }

    labels = new uint32_t[vertex_count+1]();
    std::ifstream labels_file((data_graph_path + "/labels.bin").c_str(), std::ios::binary);
    if (labels_file)
    {
      uint32_t min_label = 0;
      uint32_t max_label = 0;
      labelled_graph = true;

      for (uint32_t i = 0; i < vertex_count; ++i)
      {
        uint32_t buf[2];
        labels_file.read(reinterpret_cast<char *>(buf), 2*sizeof(uint32_t));

        uint32_t v = buf[0];
        uint32_t label = buf[1];
        min_label = std::min(min_label, label);
        max_label = std::max(max_label, label);
        labels[v] = label;
      }

      label_range = {min_label, max_label};
    }

    {
      std::ifstream ids_file((data_graph_path + "/ids.bin").c_str(), std::ios::binary);
      assert(ids_file);
      ids = new uint32_t[vertex_count+1]();
      for (uint32_t i = 1; i <= vertex_count; ++i)
      {
        ids_file.read(reinterpret_cast<char *>(&ids[i]), sizeof(uint32_t));
      }
    }
  }

  DataGraph::DataGraph(DataGraph &&other)
    : rbi(other.rbi),
      vertex_count(other.vertex_count),
      edge_count(other.edge_count),
      forest_count(other.forest_count),
      labelled_graph(other.labelled_graph),
      labels(other.labels),
      label_range(other.label_range),
      ids(other.ids),
      data_graph(other.data_graph),
      graph_in_memory(other.graph_in_memory),
      known_labels(other.known_labels)
  {
    other.labels = nullptr;
    other.data_graph = nullptr;
    other.graph_in_memory = nullptr;
    other.vertex_count = 0;
    other.edge_count = 0;
  }

  DataGraph::~DataGraph() {
    delete[] graph_in_memory;
    delete[] data_graph;
    delete[] labels;
  }

  void DataGraph::set_rbi(const AnalyzedPattern &new_rbi)
  {
    rbi = new_rbi;
    forest_count = rbi.vgs.size();
    if (rbi.labelling_type() == Graph::PARTIALLY_LABELLED)
    {
      new_label = std::distance(rbi.query_graph.labels.cbegin(),
          std::find(rbi.query_graph.labels.cbegin(),
            rbi.query_graph.labels.cend(), static_cast<uint32_t>(-1))
          );
    }
  }

  void DataGraph::set_known_labels(const std::vector<SmallGraph> &patterns)
  {
    for (const auto &p : patterns)
    {
      known_labels.insert(p.labels.cbegin(), p.labels.cend());
    }
  }

  void DataGraph::set_known_labels(const std::vector<uint32_t> &labels)
  {
    known_labels.insert(labels.cbegin(), labels.cend());
  }

  bool DataGraph::known_label(const uint32_t label) const
  {
    return known_labels.contains(label);
  }

  const std::vector<uint32_t> &DataGraph::get_upper_bounds(uint32_t v) const
  {
    return rbi.upper_bounds[v];
  }

  const std::vector<uint32_t> &DataGraph::get_lower_bounds(uint32_t v) const
  {
    return rbi.lower_bounds[v];
  }

  const std::pair<uint32_t, uint32_t> &DataGraph::get_label_range() const
  {
    return label_range;
  }

  const adjlist &DataGraph::get_adj(uint32_t v) const
  {
    return data_graph[v-1];
  }

  const SmallGraph &DataGraph::get_vgs(unsigned vgsi) const
  {
    return rbi.vgs[vgsi];
  }

  uint32_t DataGraph::original_id(uint32_t v) const
  {
    return ids[v];
  }

  const SmallGraph &DataGraph::get_pattern() const {
    return rbi.query_graph;
  }

  const std::vector<std::vector<uint32_t>> &DataGraph::get_qs(unsigned vgsi) const
  {
    return rbi.qs[vgsi];
  }

  uint32_t DataGraph::label(uint32_t dv) const
  {
    return labels[dv];
  }

  uint32_t DataGraph::vmap_at(unsigned vgsi, uint32_t v, unsigned qsi) const
  {
    return rbi.vmap[vgsi][v][qsi];
  }

  const std::vector<uint32_t> &DataGraph::get_qo(uint32_t vgsi) const
  {
    return rbi.qo_book[vgsi];
  }

  uint32_t DataGraph::get_vgs_count() const
  {
    return rbi.vgs.size();
  }

  uint32_t DataGraph::get_vertex_count() const
  {
    return vertex_count;
  }

  uint64_t DataGraph::get_edge_count() const
  {
    return edge_count;
  }
}
