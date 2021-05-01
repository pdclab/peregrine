#include <thread>
#include <atomic>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <execution>

#include <cassert>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "utils.hh"
#include "DataConverter.hh"

namespace Peregrine
{
  namespace DataConverter
  {
    const unsigned nthreads = std::thread::hardware_concurrency();
    bool is_directory(const std::string &path)
    {
      struct stat statbuf;
      return (stat(path.c_str(), &statbuf) == 0) && S_ISDIR(statbuf.st_mode);
    }
    
    bool file_exists(const std::string &path)
    {
      struct stat statbuf;
      return (stat(path.c_str(), &statbuf) == 0);
    }
    
    struct adjlist {
      adjlist() : length(0), ptr(nullptr) {}
      adjlist(uint32_t l, uint32_t *p) : length(l), ptr(p) {}
      std::atomic<uint32_t> length;
      uint32_t *ptr;
    };
    
    
    /**
     * task_size gives the approx. number of bytes to process per thread.
     * This is adjusted to the nearest line.
     */
    void calculate_degree_map(unsigned thread_id, const char *graph, size_t file_size, size_t task_size, uint64_t &num_edges, std::vector<uint32_t> &degree_maps, uint32_t &max_vid)
    {
      // estimate number of vertices based on file_size
      uint32_t vertex_count_estimate = 4096;
      if (file_size / (10lu*1024*1024*1024) != 0) // 10's of GB
      {
        vertex_count_estimate = 100'000'000;
      }
      else if (file_size / (1024*1024*1024) != 0) // GBs
      {
        vertex_count_estimate = 10'000'000;
      }
      else if (file_size / (100*1024*1024) != 0) // 100's of MBs
      {
        vertex_count_estimate = 1'000'000;
      }
      else if (file_size / (1024*1024) != 0) // MBs
      {
        vertex_count_estimate = 100'000;
      }

      std::vector<uint32_t> degree_map(vertex_count_estimate);
      uint64_t edges = 0;
      uint32_t max = 0;
    
      size_t start = thread_id * task_size;
      size_t end = start + task_size;

      if (thread_id == nthreads - 1)
      {
        end = file_size;
      }
    
      // don't start from middle of a line in the file
      if (thread_id != 0)
      {
        while (start <= file_size && graph[start-1] != '\n') ++start;
      }
    
      uint32_t u, v;
      size_t cursor = start;
      while (cursor < end)
      {
        // find the end of this edge
        size_t line_end = cursor;
        while (line_end < file_size && graph[line_end] != '\n') ++line_end;
    
        if (line_end >= file_size) break;
    
        // read the vertices from it
        std::string line(&graph[cursor], &graph[line_end]);
    
        cursor = line_end+1;
    
        if (line[0] == '#') continue;
    
        std::istringstream iss(line);
        iss >> u >> v;
    
        max = std::max(max, std::max(u, v));
    
        // make sure we have space in the degree map
        if (max >= degree_map.size()) degree_map.resize(max+1);
    
        // ignore loops
        if (u != v)
        {
          // count edges you've seen
          edges += 1;
          degree_map[u] += 1;
          degree_map[v] += 1;
        }
      }
    
      // return thread-local degree map
      num_edges = edges;
      max_vid = max;
      degree_maps.swap(degree_map);
    }
    
    /**
     * task_size gives the approx. number of bytes to process per thread.
     * This is adjusted to the nearest line.
     */
    void populate_graph(unsigned thread_id, const char *edges, size_t file_size, size_t task_size, adjlist *graph, const std::vector<uint32_t> &ids_rev_map)
    {
      size_t start = thread_id * task_size;
      size_t end = start + task_size;

      if (thread_id == nthreads - 1)
      {
        end = file_size;
      }
    
      // don't start from middle of a line in the file
      if (thread_id != 0)
      {
        while (start <= file_size && edges[start-1] != '\n') ++start;
      }
    
      uint32_t u, v;
      size_t cursor = start;
      while (cursor < end)
      {
        // find the end of this edge
        size_t line_end = cursor;
        while (line_end < file_size && edges[line_end] != '\n') ++line_end;
    
        if (line_end >= file_size) break;
    
        // read the vertices from it
        std::string line(&edges[cursor], &edges[line_end]);
    
        cursor = line_end+1;
    
        if (line[0] == '#') continue;
    
        std::istringstream iss(line);
        iss >> u >> v;
    
        // ignore loops
        if (u != v)
        {
          uint32_t uu = ids_rev_map[u];
          uint32_t vv = ids_rev_map[v];
          uint32_t u_pos = graph[uu].length++;
          // have to live with the 1-based vertex ID's forever...
          graph[uu].ptr[u_pos] = vv+1;
          uint32_t v_pos = graph[vv].length++;
          graph[vv].ptr[v_pos] = uu+1;
        }
      }
    }
    
    void write_local(unsigned thread_id, const adjlist *graph, uint32_t num_vertices, uint32_t task_size, const std::string &out_dir)
    {
    
      uint32_t start = thread_id * task_size;
      uint32_t end = std::min(start + task_size, num_vertices);
    
      if (thread_id == nthreads - 1)
      {
        end = num_vertices;
      }
    
      std::string out_path(out_dir + "/data_" + std::to_string(thread_id) + ".bin");
      std::ofstream output(out_path.c_str(), std::ios::binary | std::ios::trunc);
    
      for (uint32_t v = start; v < end; ++v)
      {
        uint32_t deg = graph[v].length.load(std::memory_order_relaxed);
        //output << reinterpret_cast<const char *>(v) << deg;
        //output.write(reinterpret_cast<const char *>(v), sizeof(v));
        output.write(reinterpret_cast<const char *>(&deg), sizeof(deg));
        output.write(reinterpret_cast<const char *>(graph[v].ptr), deg*sizeof(uint32_t));
      }
    }
    
    void convert_data(const std::string &edge_file, const std::string &label_file, const std::string &out_dir)
    {
      auto t1 = utils::get_timestamp();
    
      struct stat edge_file_stat;
      int edge_fd = open(edge_file.c_str(), O_RDONLY, 0);
      assert(edge_fd != -1);
    
      assert(fstat(edge_fd, &edge_file_stat) == 0);
    
      size_t edge_file_size = edge_file_stat.st_size;
      char *graph_data = static_cast<char*>(mmap(NULL,
            edge_file_size,
            PROT_READ,
            MAP_PRIVATE | MAP_POPULATE,
            edge_fd,
            0));
    
      auto t2 = utils::get_timestamp();
      utils::Log{} << "Read edge-list in " << (t2-t1)/1e6 << "s" << "\n";
    
    
      // calculate degree maps
      std::vector<uint64_t> edge_counts(nthreads);
      std::vector<std::vector<uint32_t>> degree_maps(nthreads);
      std::vector<uint32_t> max_vids(nthreads);
      auto t3 = utils::get_timestamp();
      {
        size_t task_size = edge_file_size / nthreads;
        std::vector<std::thread> pool;
        for (unsigned i = 0; i < nthreads; ++i)
        {
          pool.emplace_back(calculate_degree_map, i, graph_data, edge_file_size,
              task_size,
              std::ref(edge_counts[i]),
              std::ref(degree_maps[i]),
              std::ref(max_vids[i]));
        }
    
        for (auto &th : pool)
        {
          th.join();
        }
      }
    
      uint32_t max_vid = *std::max_element(max_vids.cbegin(), max_vids.cend());
    
      // find biggest per-thread degree map to use as the final degree map
      std::vector<uint32_t> &degree_map = *std::max_element(degree_maps.begin(),
          degree_maps.end(),
          [](auto &v1, auto &v2) { return v1.size() < v2.size(); });
    
      degree_map.resize(max_vid+1);
    
      // sum up degrees
      {
        std::for_each(std::execution::par_unseq, degree_map.begin(), degree_map.end(),
            [&degree_maps, &degree_map](uint32_t &total)
            {
              uint32_t v = &total - &degree_map[0];
              uint32_t prev_total = total;
              total = 0;
              for (auto &map : degree_maps)
              {
                total += map[v];
              }
              total += prev_total;
            });
      }
    
      uint64_t num_edges = std::accumulate(edge_counts.cbegin(), edge_counts.cend(), 0);
    
      auto t4 = utils::get_timestamp();
      utils::Log{} << "Calculated degree map in " << (t4-t3)/1e6 << "s" << "\n";
    
      auto t5 = utils::get_timestamp();
    
      // initialize map from new ID to original ID
      std::vector<uint32_t> ids_map(degree_map.size());
      std::iota(ids_map.begin(), ids_map.end(), 0);
    
      // sort by degree
      std::sort(std::execution::par_unseq, ids_map.begin(), ids_map.end(),
          [&degree_map](uint32_t u, uint32_t v)
          {
            return degree_map[u] > degree_map[v];
          });
    
      // remove degree 0 vertices
      ids_map.erase(std::remove_if(ids_map.begin(), ids_map.end(),
            [&degree_map](uint32_t v) { return degree_map[v] == 0; }),
          ids_map.end());
      uint32_t num_vertices = ids_map.size();
    
    
      // map from original ID to new ID
      std::vector<uint32_t> ids_rev_map(degree_map.size(), static_cast<uint32_t>(-1));
      {
        std::for_each(std::execution::par_unseq, ids_map.cbegin(), ids_map.cend(),
          [&ids_rev_map, &ids_map](const uint32_t &true_v)
          {
            uint32_t v = &true_v - &ids_map[0];
            ids_rev_map[true_v] = v;
          });
      }
    
      auto t6 = utils::get_timestamp();
      utils::Log{} << "Sorted vertices in " << (t6-t5)/1e6 << "s" << "\n";
    
      // initialize adjacency lists in memory
      auto t7 = utils::get_timestamp();
      adjlist *graph = new adjlist[num_vertices];
      std::for_each(std::execution::par_unseq, ids_map.cbegin(), ids_map.cend(),
          [&graph, &ids_rev_map, &degree_map](const uint32_t true_v)
          {
            uint32_t new_v = ids_rev_map[true_v];
            uint32_t deg = degree_map[true_v];
            graph[new_v].ptr = new uint32_t[deg];
            graph[new_v].length = 0;
          });
    
      // clear unused memory
      degree_maps.clear();
    
      // populate the in-memory graph
      {
        size_t task_size = edge_file_size / nthreads;
        std::vector<std::thread> pool;
        for (unsigned i = 0; i < nthreads; ++i)
        {
          pool.emplace_back(populate_graph, i, graph_data, edge_file_size, task_size, graph, std::cref(ids_rev_map));
        }
    
        for (auto &th : pool)
        {
          th.join();
        }
    
        // don't need the edges anymore
        munmap(graph_data, edge_file_size);
      }
      auto t8 = utils::get_timestamp();
      utils::Log{} << "Created in-memory graph in " << (t8-t7)/1e6 << "s" << "\n";
    
      // sort the in-memory adjacency lists
      auto t9 = utils::get_timestamp();
      std::for_each(std::execution::par_unseq, ids_map.cbegin(), ids_map.cend(),
          [&graph, &ids_rev_map](uint32_t true_v)
          {
            uint32_t v = ids_rev_map[true_v];
            std::sort(std::execution::unseq, graph[v].ptr, graph[v].ptr + graph[v].length);
          });
      auto t10 = utils::get_timestamp();
      utils::Log{} << "Sorted adjacency lists in " << (t10-t9)/1e6 << "s" << "\n";
    
      // write thread-local files
      auto t11 = utils::get_timestamp();
      {
        uint32_t task_size =  num_vertices / nthreads;
        std::vector<std::thread> pool;
        for (unsigned i = 0; i < nthreads; ++i)
        {
          pool.emplace_back(write_local, i, graph, num_vertices, task_size, std::ref(out_dir));
        }
    
        for (auto &th : pool)
        {
          th.join();
        }
      }
      auto t12 = utils::get_timestamp();
      utils::Log{} << "Wrote data to disk in " << (t12-t11)/1e6 << "s" << "\n";
    
      auto t13 = utils::get_timestamp();

      // clean up in memory graph
      std::for_each(std::execution::unseq, ids_map.cbegin(), ids_map.cend(),
          [&graph, &ids_rev_map](const uint32_t true_v)
          {
            uint32_t new_v = ids_rev_map.at(true_v);
            delete[] graph[new_v].ptr;
          });
      delete[] graph;

      // concatenate thread local files
      {
        std::string output_path(out_dir + "/data.bin");
        {
          std::ofstream output(output_path.c_str(), std::ios::binary | std::ios::trunc);
          // write header first:
          // - number of vertices in the graph
          // - number of edges in the graph
          output.write(reinterpret_cast<const char *>(&num_vertices), sizeof(num_vertices));
          output.write(reinterpret_cast<const char *>(&num_edges), sizeof(num_edges));
        }
    
        // now merge the thread-local files, deleting as you go
        for (unsigned i = 0; i < nthreads; ++i)
        {
          std::string thread_local_path(out_dir + "/data_" + std::to_string(i) + ".bin");
          {
            std::ifstream input(thread_local_path.c_str(), std::ios::binary);
            std::ofstream output(output_path.c_str(), std::ios::binary | std::ios::app);
            output << input.rdbuf();
          }
          std::remove(thread_local_path.c_str());
        }
      }
      auto t14 = utils::get_timestamp();
      utils::Log{} << "Concatenated output in " << (t14-t13)/1e6 << "s" << "\n";
    
      // associate labels if they exist
      if (file_exists(label_file))
      {
        auto t15 = utils::get_timestamp();
        {
          std::ifstream ifile(label_file.c_str());
    
          std::string output_labels(out_dir + "/labels.bin");
          std::ofstream ofile(output_labels.c_str(), std::ios::binary | std::ios::trunc);
    
          std::string line;
          while (std::getline(ifile, line))
          {
            // easy to predict properly
            if (line[0] == '#') continue;
            std::istringstream iss(line);
    
            uint32_t u;
            uint32_t new_u_label[2];
            iss >> u >> new_u_label[1];
            new_u_label[0] = ids_rev_map[u]+1;
    
            if (new_u_label[0] != 0) // ids_rev_map[u] != -1
            {
              ofile.write(reinterpret_cast<char *>(&new_u_label[0]), 2*sizeof(uint32_t));
            }
          }
        }
        auto t16 = utils::get_timestamp();
        utils::Log{} << "Converted labels to binary format in " << (t16-t15)/1e6 << "s" << "\n";
      }
    
      // write ID map
      auto t17 = utils::get_timestamp();
      {
        std::string output_ids(out_dir + "/ids.bin");
        std::ofstream file(output_ids.c_str(), std::ios::binary | std::ios::trunc);
        file.write(reinterpret_cast<char *>(&ids_map[0]), ids_map.size() * sizeof(uint32_t));
      }
      auto t18 = utils::get_timestamp();
      utils::Log{} << "ID map written in " << (t18-t17)/1e6 << "s" << "\n";
    
      utils::Log{} << "Finished in " << (t18-t1)/1e6 << "s" << "\n";
    }
  }
}
