#include <iostream>

#include <string>
#include <cstring>
#include <algorithm>
#include <unordered_set>

#include <unistd.h>

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

    page_graph_info.vertex_count = p.num_vertices();
    page_graph_info.edge_count = p.num_true_edges();

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

      FILE* input_graph = fopen((data_graph_path+"/data.bin").c_str(), "rb+");
      if (!input_graph)
      {
        std::cerr << "Opening page graph file for reading failed" << std::endl;
        exit(1);
      }

      uint64_t offset = 0;
      uint32_t data_count;
      pread(fileno(input_graph),&data_count,4,offset);
      graph_in_memory = new uint32_t[data_count];
      offset+=4;
      pread(fileno(input_graph),&page_graph_info.vertex_count,4,offset);
      offset+=4;
      pread(fileno(input_graph),&page_graph_info.edge_count,4,offset);
      offset+=4;
      uint64_t graph_size = data_count;
      graph_size*=4;
      uint64_t curr_read_offset = 0;
      uint32_t MAX_READ_SIZE = 2147479552;
      while(true){
        if(curr_read_offset + MAX_READ_SIZE > graph_size){
          pread(fileno(input_graph),graph_in_memory+(curr_read_offset/4),graph_size-curr_read_offset,curr_read_offset);
          break;
        }
        else{
          pread(fileno(input_graph),graph_in_memory+(curr_read_offset/4),MAX_READ_SIZE,curr_read_offset);
          curr_read_offset+=(MAX_READ_SIZE);
        }
      }

      data_graph = new adjlist[page_graph_info.vertex_count];

      for(uint32_t i=0;i<page_graph_info.vertex_count;i++){
          uint32_t curr_vertex_offset = graph_in_memory[3+i];
          data_graph[i].length = graph_in_memory[curr_vertex_offset+1];
          data_graph[i].ptr = &graph_in_memory[curr_vertex_offset+2];
      }


      page_graph_info.page_size = data_graph[0].length;
      labels = new uint32_t[page_graph_info.vertex_count+1]();
      FILE *labels_file = fopen((data_graph_path + "/labels.bin").c_str(), "rb");
      if (labels_file) {
        labelled_graph = true;
        for (uint32_t i = 0; i < page_graph_info.vertex_count; ++i) {
          uint32_t buf[2];
          read(fileno(labels_file), buf, 2*4);
          uint32_t v = buf[0];
          uint32_t label = buf[1];
          labels[v] = label;
        }
        fclose(labels_file);
      } else {
        std::cout << "No label file" << std::endl;
      }
  }

  DataGraph::DataGraph(const DataGraph *other)
    : page_graph_info(other->page_graph_info),
      labelled_graph(other->labelled_graph),
      labels(other->labels),
      data_graph(other->data_graph),
      known_labels(other->known_labels)
  {}

  DataGraph::DataGraph(DataGraph &&other)
    : page_graph_info(other.page_graph_info),
      rbi(other.rbi),
      forest_count(other.forest_count),
      labelled_graph(other.labelled_graph),
      labels(other.labels),
      data_graph(other.data_graph),
      graph_in_memory(other.graph_in_memory),
      known_labels(other.known_labels)
  {
    other.labels = nullptr;
    other.data_graph = nullptr;
    other.graph_in_memory = nullptr;
    other.page_graph_info.vertex_count = 0;
    other.page_graph_info.edge_count = 0;
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

  const std::vector<uint32_t> &DataGraph::get_upper_bounds(uint32_t v) const {
    return rbi.upper_bounds[v];
  }

  const std::vector<uint32_t> &DataGraph::get_lower_bounds(uint32_t v) const {
    return rbi.lower_bounds[v];
  }

  const adjlist &DataGraph::get_adj(uint32_t v) const {
    return data_graph[v-1];
  }

  const SmallGraph &DataGraph::get_vgs(unsigned vgsi) const {
    return rbi.vgs[vgsi];
  }

  const SmallGraph &DataGraph::get_pattern() const {
    return rbi.query_graph;
  }

  const std::vector<std::vector<uint32_t>> &DataGraph::get_qs(unsigned vgsi) const {
    return rbi.qs[vgsi];
  }

  uint32_t DataGraph::label(uint32_t dv) const {
    return labels[dv];
  }

  uint32_t DataGraph::vmap_at(unsigned vgsi, uint32_t v, unsigned qsi) const {
    return rbi.vmap[vgsi][v][qsi];
  }

  const std::vector<uint32_t> &DataGraph::get_qo(uint32_t vgsi) const {
    return rbi.qo_book[vgsi];
  }

  uint32_t DataGraph::get_vgs_count() const{
    return rbi.vgs.size();
  }
}
