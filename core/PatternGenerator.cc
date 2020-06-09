#include "PatternGenerator.hh"

#include <unordered_set>
#include <iostream>

namespace Peregrine
{
  std::vector<std::vector<std::pair<uint32_t, uint32_t>>>
  get_elists(uint32_t size);

  SmallGraph PatternGenerator::clique(uint32_t size) {
    std::vector<std::pair<uint32_t, uint32_t>> edge_list;
    for (uint32_t u = 1; u <= size; ++u) {
      for (uint32_t v = 1; v <= size; ++v) {
        if (u < v) {
          edge_list.push_back({u, v});
        }
      }
    }

    return SmallGraph(edge_list);
  }

  SmallGraph PatternGenerator::star(uint32_t size) {
    if (size < 3) {
      std::cerr << "Please pass in size > 2, we don't do pattern-matching on single edges/vertices" << std::endl;
    }

    std::vector<std::pair<uint32_t, uint32_t>> edge_list;
    for (uint32_t u = 2; u <= size; ++u) {
      edge_list.emplace_back(1, u);
    }

    return SmallGraph(edge_list);
  }

  /**
   * Undirected edge-induced extension:
   *  - generate all extensions containing "known" labels as fully-labelled
   *    patterns
   *  - generate all extensions NOT containing "known" labels as
   *    partially-labelled patterns
   * where "known" means existing in the unextended pattern.
   * For example, the 3-star
   *    1 a 2 3
   *    2 b 1
   *    3 c 1
   * will be extended to the following 4-stars:
   *    1 a 2 3 4
   *    2 b 1
   *    3 c 1
   *    4 a 1
   *
   *    1 a 2 3 4
   *    2 b 1
   *    3 c 1
   *    4 b 1
   *
   *    1 a 2 3 4
   *    2 b 1
   *    3 c 1
   *    4 c 1
   *
   *    1 a 2 3 4
   *    2 b 1
   *    3 c 1
   *    4 - 1
   *
   * Notice that only one (the last) is partially-labelled, with 3 fully-labelled
   * ones where the added vertex takes one of the pre-existing 3 labels.
   */
  std::vector<SmallGraph> PatternGenerator::extend(const std::vector<SmallGraph> &from, bool vertex_based, bool overwrite_anti_edges)
  {
    std::vector<SmallGraph> result;
    std::unordered_set<uint32_t> label_set;
    for (const auto &p : from)
    {
      // if labels
      label_set.insert(p.labels.begin(), p.labels.end());
    }

    // partially-labelled
    label_set.insert(static_cast<uint32_t>(-1));

    if (vertex_based)
    {
      for (const auto &p : from)
      {

        // indsets do not guarantee uniqueness in the result set!
        // consider an unlabelled 4chain:
        // the partial order on the centre points means there is no order on the
        // end points, hence they are in different indsets. However, the end
        // points are also symmetric.
        // Have to merge indsets which have symmetric parents maybe?
        const AnalyzedPattern rbi(p);
        const auto &indsets = rbi.indsets;

        uint32_t start_idx = result.size();

        const uint32_t new_v = p.num_vertices()+1;

        const auto &adj = p.true_adj_list;

        std::vector<uint32_t> perms(indsets.size() - 1, 0);
        uint32_t idx = indsets.size()-1; // begin at the end
        while (true)
        {
          if (idx == indsets.size() - 1)
          {
            // generate all the extensions

            auto new_adj(adj);

            for (uint32_t k = 0; k < perms.size(); ++k)
            {
              // add an edge from new_v to the first perms[k] elements in kth indset
              for (uint32_t j = 0; j < perms[k]; ++j)
              {
                new_adj[new_v].push_back(indsets[k][j]);
                new_adj[indsets[k][j]].push_back(new_v);
              }
            }

            // if an edge has been added from another indset, then make an
            // extension with no edge from the last indset
            if (std::any_of(perms.begin(), perms.end(), [](uint32_t n) { return n != 0; }))
            {
              if (p.labelling == Graph::LABELLED)
              {
                for (uint32_t label : label_set)
                {
                  std::vector<uint32_t> labels(p.labels);
                  labels.push_back(label);
                  result.emplace_back(new_adj, labels);
                }
              }
              else
              {
                result.emplace_back(new_adj);
              }
            }

            // add an edge from new_v to the all elements in last indset,
            // finishing an extension with each edge
            for (uint32_t j = 0; j < indsets.back().size(); ++j)
            {
              new_adj[new_v].push_back(indsets.back()[j]);
              new_adj[indsets.back()[j]].push_back(new_v);

              // finished one extension
              if (p.labelling == Graph::LABELLED)
              {
                for (uint32_t label : label_set)
                {
                  std::vector<uint32_t> labels(p.labels);
                  labels.push_back(label);
                  result.emplace_back(new_adj, labels);
                }
              }
              else
              {
                result.emplace_back(new_adj);
              }
            }

            idx -= 1;
            perms[idx] += 1;
          }
          else if (perms[idx] <= indsets[idx].size())
          {
            idx += 1;
          }
          else if (idx != 0)
          {
            perms[idx] = 0;
            idx -= 1;
            perms[idx] += 1;
          }
          else
          {
            break;
          }
        }


        // add anti-vertices and anti-edges back into all extensions of this pattern
        const auto &anti_vertices = rbi.anti_vertices;
        std::vector<std::vector<SmallGraph>::iterator> to_delete;
        for (uint32_t i = start_idx; i < result.size(); ++i)
        {
          auto &rp = result[i];

          for (const auto &[u, nbrs] : p.anti_adj_list)
          {
            for (uint32_t v : nbrs)
            {
              if (u > v || p.is_anti_vertex(u) || p.is_anti_vertex(v)) continue;

              if (!utils::search(rp.get_neighbours(u), v)) // hasn't been overwritten
              {
                rp.add_anti_edge(u, v);
              }
              else if (!overwrite_anti_edges)
              {
                // should be maintaining all anti-edges: delete this pattern
                to_delete.push_back(result.begin() + i);
              }
            }
          }

          uint32_t new_av = rp.num_vertices() + 1;
          for (uint32_t av : anti_vertices)
          {
            for (uint32_t regular_vert : p.get_anti_neighbours(av))
            {
              rp.add_anti_edge(new_av, regular_vert);
            }
            new_av += 1;
          }
        }

        for (auto it : to_delete) result.erase(it);
      }
    }
    else // edge-based
    {
      for (const auto &p : from)
      {
        const AnalyzedPattern rbi(p);
        const auto &indsets = rbi.indsets;

        uint32_t start_idx = result.size();

        const auto &adj = p.true_adj_list;

        // add edge between first pair of disconnected vertices in an indset
        for (const auto &indset : indsets)
        {
          for (uint32_t i = 0; i < indset.size(); ++i)
          {
            bool found = false;
            uint32_t u = indset[i];
            for (uint32_t missing = i+1; missing < indset.size(); ++missing)
            {
              uint32_t v = indset[missing];
              if (!utils::search(adj.at(u), v))
              {
                // add edge here
                auto new_adj(adj);
                new_adj[u].push_back(v);
                new_adj[v].push_back(u);
                result.emplace_back(new_adj, p.labels);

                // XXX: ugly. It's making sure this isn't marked labelled falsely
                // better to get rid of all the labels.resize() stuff in Graph.hh
                if (p.labelling != Graph::LABELLED) result.back().labelling = Graph::UNLABELLED;

                // break out of nested loop
                found = true;
                break;
              }
            }
            if (found) break;
          }
        }

        // TODO: document this properly somewhere
        // add edge between disconnected vertices in different indsets if:
        //  * src vertex's indset id < dst vertex's indset id
        //  * src vertex is the first in its indset to not be connected to dst
        //  * dst vertex is the first in its indset to not be connected to src
        for (uint32_t i = 0; i < indsets.size(); ++i)
        {
          const auto &indset = indsets[i];
          // for every other indset
          for (uint32_t other = i+1; other < indsets.size(); ++other)
          {
            bool found = false;
            const auto &other_indset = indsets[other];
            // check for missing edges between a vertex in this indset and one in the other
            for (uint32_t u : indset)
            {
              for (uint32_t v : other_indset)
              {
                if (!utils::search(adj.at(u), v))
                {
                  // add the edge
                  auto new_adj(adj);
                  new_adj[u].push_back(v);
                  new_adj[v].push_back(u);
                  result.emplace_back(new_adj, p.labels);
                  // XXX: ugly. It's making sure this isn't marked labelled falsely
                  // better to just get rid of all the labels.resize() stuff in Graph::*
                  if (p.labelling != Graph::LABELLED) result.back().labelling = Graph::UNLABELLED;

                  // break out of nested loop
                  found = true;
                  break;
                }
              }
              if (found) break;
            }
          }
        }

        // add new vertex

        if (p.labelling != Graph::LABELLED)
        {
          // easy case: Only one valid extension per symmetric vertex set
          uint32_t new_vid = p.num_vertices() + 1;
          for (const auto &indset : indsets)
          {
            uint32_t src = indset.front();
            auto new_adj(adj);
            new_adj[src].push_back(new_vid);
            new_adj[new_vid].push_back(src);
            result.emplace_back(new_adj);
          }
        }
        else if (p.labelling == Graph::LABELLED)
        {
          // hard case: If n distinct labels, then will get n labelled extensions
          // and 1 partially-labelled extension per symmetric vertex set

          //std::unordered_set<uint32_t> label_set(p.labels.begin(), p.labels.end());

          // partially-labelled
          //label_set.insert(static_cast<uint32_t>(-1));

          uint32_t new_vid = p.num_vertices() + 1;
          for (const auto &indset : indsets)
          {
            uint32_t src = indset.front();
            auto new_adj(adj);
            new_adj[src].push_back(new_vid);
            new_adj[new_vid].push_back(src);
            // extend with every possible label
            for (uint32_t l : label_set)
            {
              std::vector<uint32_t> labels(p.labels);
              labels.push_back(l);
              result.emplace_back(new_adj, labels);
            }
          }
        }
        else
        {
          // XXX: or can I?
          std::cerr
            << "WARN: cannot extend partially labelled patterns with new partially labelled vertices!"
            << "\nOnly adding edges between existing vertices."
            << std::endl;
        }

        // add anti-vertices back into all extensions of this pattern
        const auto &anti_vertices = rbi.anti_vertices;
        std::vector<std::vector<SmallGraph>::iterator> to_delete;
        for (uint32_t i = start_idx; i < result.size(); ++i)
        {
          auto &rp = result[i];

          for (const auto &[u, nbrs] : p.anti_adj_list)
          {
            for (uint32_t v : nbrs)
            {
              if (u > v || p.is_anti_vertex(u) || p.is_anti_vertex(v)) continue;

              if (!utils::search(rp.get_neighbours(u), v)) // hasn't been overwritten
              {
                rp.add_anti_edge(u, v);
              }
              else if (!overwrite_anti_edges)
              {
                // should be maintaining all anti-edges: delete this pattern
                to_delete.push_back(result.begin() + i);
              }
            }
          }

          uint32_t new_av = rp.num_vertices() + 1;
          for (uint32_t av : anti_vertices)
          {
            for (uint32_t regular_vert : p.get_anti_neighbours(av))
            {
              rp.add_anti_edge(new_av, regular_vert);
            }
            new_av += 1;
          }
        }

        for (auto it : to_delete) result.erase(it);

      }
    }


    std::unordered_set<uint64_t> hashes;
    std::erase_if(result, [&hashes](const SmallGraph &p)
    {
      // erase if hash exists, else insert hash
      uint64_t h = p.bliss_hash();

      if (!hashes.contains(h))
      {
        hashes.insert(h);
        return false;
      }
      else
      {
        return true;
      }
    });

    return result;
  }

  std::vector<SmallGraph> PatternGenerator::all(uint32_t size, bool vertex_based, bool anti_edges) {
    if (size > 9) {
      std::cerr << "Currently don't support generating all patterns of size greater than 9" << std::endl;
    }

    std::vector<SmallGraph> result;
    if (vertex_based) {
      const auto &elists = get_elists(size);
      for (const auto &el : elists) {
        result.emplace_back(el);
      }
    } else if (size > 2) {
      const auto &elists1 = get_elists(size);
      for (const auto &el : elists1) {
        if (el.size() == size) {
          result.emplace_back(el);
        }
      }

      const auto &elists2 = get_elists(size);
      for (const auto &el : elists2) {
        if (el.size() == size) {
          result.emplace_back(el);
        }
      }
    } else {
      result.push_back(star(3));
    }

    if (anti_edges)
    {
      for (auto &p : result)
      {
        uint32_t nv = p.num_vertices();
        for (uint32_t u = 1; u < nv; ++u)
        {
          for (uint32_t v = u + 1; v <= nv; ++v)
          {
            if (!utils::search(p.get_neighbours(u), v))
            {
              p.add_anti_edge(u, v);
            }
          }
        }
      }
    }

    return result;
  }

  bool PatternGenerator::is_connected_pattern(Edges edge_list){
    uint32_t psize = edge_list.size();
    std::vector<std::vector<uint32_t>> graph(psize+1);
    for (auto &edge : edge_list)
    {
        graph[edge.first].push_back(edge.second);
        graph[edge.second].push_back(edge.first);
    }
    std::vector<bool> visited(psize + 1, true);

    for (uint32_t i = 1;i<=psize;i++)
    {
        visited[i] = false;
    }
    std::queue<uint32_t> q;
    q.push(1);
    visited[1] = true;
    while (!q.empty())
    {
        uint32_t current_node = q.front();
        q.pop();
        for (auto vertex_id : graph[current_node])
        {
            if (!visited[vertex_id])
            {
                q.push(vertex_id);
                visited[vertex_id] = true;
            }
        }
    }
    for (size_t i = 1; i <= psize; i++)
    {
        if (!visited[i])
            return false;
    }

    return true;
  }
}
