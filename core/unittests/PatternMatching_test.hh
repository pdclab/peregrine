#include <iostream>

#include "../Options.hh"
#include "../utils.hh"
#include "../Graph.hh"
#include "../PatternGenerator.hh"
#include "../PatternMatching.hh"

using namespace Peregrine;
// assumes patterns are sorted in ascending order by number of edges
// for each pattern, calculate the vertex-based count
std::vector<uint64_t> convert_counts(std::vector<SmallGraph> &ps, const std::vector<uint64_t> &edge_based) {
  std::vector<uint64_t> vbased(edge_based.size());

  for (int32_t i = ps.size()-1; i >= 0; --i) {
    uint64_t count = edge_based[i];
    for (uint32_t j = i+1; j < ps.size(); ++j) {
      uint32_t n = Peregrine::num_mappings(ps[j], ps[i]);
      uint64_t inc = n * vbased[j];
      std::cout << "mapping " << ps[i].to_string() << " into " << ps[j].to_string() << ": " << n << std::endl;
      count -= inc;
    }
    vbased[i] = count;
  }

  return vbased;
}

static const SmallGraph s3 = PatternGenerator::star(3);
static const SmallGraph s4 = PatternGenerator::star(4);
static const SmallGraph cl3 = PatternGenerator::clique(3);
static const SmallGraph cl4 = PatternGenerator::clique(4);
static const SmallGraph c4({
    std::make_pair(1, 2),
    std::make_pair(1, 4),
    std::make_pair(2, 3)
});
static const SmallGraph twe({
    std::make_pair(1, 2),
    std::make_pair(1, 3),
    std::make_pair(1, 4),
    std::make_pair(2, 3)
});
static const SmallGraph cy4({
    std::make_pair(1, 2),
    std::make_pair(1, 4),
    std::make_pair(2, 3),
    std::make_pair(3, 4)
});
static const SmallGraph swd({
    std::make_pair(1, 2),
    std::make_pair(1, 3),
    std::make_pair(1, 4),
    std::make_pair(2, 3),
    std::make_pair(3, 4)
});

static const SmallGraph p6({
    std::make_pair(1, 4),
    std::make_pair(1, 5),
    std::make_pair(2, 3),
    std::make_pair(2, 5),
    std::make_pair(3, 5),
    std::make_pair(3, 6),
    std::make_pair(4, 5),
    std::make_pair(4, 6),
    std::make_pair(5, 6),
});

static const SmallGraph base({
    std::make_pair(3, 5),
    std::make_pair(1, 2),
    std::make_pair(1, 3),
    std::make_pair(2, 4)
});

static const SmallGraph labelling1 = SmallGraph(base)
      .set_label(1, 0)
      .set_label(2, 0)
      .set_label(3, 0)
      .set_label(4, 0)
      .set_label(5, 0);

static const SmallGraph labelling2 = SmallGraph(base)
      .set_label(1, 1)
      .set_label(2, 1)
      .set_label(3, 1)
      .set_label(4, 1)
      .set_label(5, 1);

// edge-induced counts on citeseer
std::unordered_map<SmallGraph, uint64_t> unlabelled_ground_truth_edge({
  std::make_pair(s3, 26878),  // [1-3][2-3]                     3-star
  std::make_pair(cl3, 1166),  // [1-2][2-3][1-3]                triangle
  std::make_pair(c4, 185589), // [1-2][1-4][2-3]                4-chain
  std::make_pair(s4, 250950), // [1-2][1-3][1-4]                4-star
  std::make_pair(twe, 34760), // [1-2][1-3][1-4][2-3]           twe
  std::make_pair(cy4, 6059),  // [1-2][1-4][2-3][3-4]           4-cycle
  std::make_pair(swd, 3730),  // [1-2][1-3][1-4][2-3][3-4]      4-chordal-cycle
  std::make_pair(cl4, 255),   // [1-2][1-3][1-4][2-3][2-4][3-4] 4-clique
  std::make_pair(p6, 81870),  // p6 from the paper
  // problem patterns from 4-fsm
  std::make_pair(labelling1, 54548),
  std::make_pair(labelling2, 1123126),
});


// vertex-induced counts on citeseer
std::unordered_map<SmallGraph, uint64_t> unlabelled_ground_truth_vtx({
   std::make_pair(s3, 23380),  // [1-3][2-3]                     3-star
   std::make_pair(cl3, 1166),  // [1-2][2-3][1-3]                triangle
   std::make_pair(c4, 111153), // [1-2][1-4][2-3]                4-chain
   std::make_pair(s4, 222630), // [1-2][1-3][1-4]                4-star
   std::make_pair(twe, 22900), // [1-2][1-3][1-4][2-3]           twe
   std::make_pair(cy4, 3094),  // [1-2][1-4][2-3][3-4]           4-cycle
   std::make_pair(swd, 2200),  // [1-2][1-3][1-4][2-3][3-4]      4-chordal-cycle
   std::make_pair(cl4, 255),   // [1-2][1-3][1-4][2-3][2-4][3-4] 4-clique
});


template <Graph::Labelling L = Graph::UNLABELLED, typename F>
void matchPattern(DataGraph &dg,
          const SmallGraph &p,
          F &process)
{
  std::vector<std::vector<uint32_t>> cands(p.num_vertices()+2, std::vector<uint32_t>{});
  for (auto &cand : cands) cand.reserve(10000);

  dg.set_rbi(p);
  bool has_anti_edges = dg.rbi.has_anti_edges();
  bool has_anti_vertices = !dg.rbi.anti_vertices.empty();

  uint32_t num_vertices = dg.get_vertex_count();
  uint32_t vgs_count = dg.get_vgs_count();

  if (has_anti_edges)
  {
    if (has_anti_vertices)
    {
      for (uint32_t vgsi = 0; vgsi < vgs_count; ++vgsi)
      {
        Peregrine::Matcher<true, UNSTOPPABLE, decltype(process)> m(std::stop_token(), dg.rbi, &dg, vgsi, cands, process);
        for (uint32_t v = 1; v <= num_vertices; ++v)
        {
          m.template map_into<L, true>(v);
        }
      }
    }
    else
    {
      for (uint32_t vgsi = 0; vgsi < vgs_count; ++vgsi)
      {
        Peregrine::Matcher<false, UNSTOPPABLE, decltype(process)> m(std::stop_token(), dg.rbi, &dg, vgsi, cands, process);
        for (uint32_t v = 1; v <= num_vertices; ++v)
        {
          m.template map_into<L, true>(v);
        }
      }
    }
  }
  else
  {
    if (has_anti_vertices)
    {
      for (uint32_t vgsi = 0; vgsi < vgs_count; ++vgsi)
      {
        Peregrine::Matcher<true, UNSTOPPABLE, decltype(process)> m(std::stop_token(), dg.rbi, &dg, vgsi, cands, process);
        for (uint32_t v = 1; v <= num_vertices; ++v)
        {
          m.template map_into<L, false>(v);
        }
      }
    }
    else
    {
      for (uint32_t vgsi = 0; vgsi < vgs_count; ++vgsi)
      {
        Peregrine::Matcher<false, UNSTOPPABLE, decltype(process)> m(std::stop_token(), dg.rbi, &dg, vgsi, cands, process);
        for (uint32_t v = 1; v <= num_vertices; ++v)
        {
          m.template map_into<L, false>(v);
        }
      }
    }
  }

}

template <Graph::Labelling L = Graph::UNLABELLED>
uint64_t countPattern(DataGraph &dg,
          const SmallGraph &p)
{
  std::vector<std::vector<uint32_t>> cands(p.num_vertices()+2, std::vector<uint32_t>{});
  for (auto &cand : cands) cand.reserve(10000);

  dg.set_rbi(p);
  bool has_anti_vertices = !dg.rbi.anti_vertices.empty();

  uint32_t num_vertices = dg.get_vertex_count();
  uint32_t vgs_count = dg.get_vgs_count();

  uint64_t c = 0;

  if (has_anti_vertices)
  {
    for (uint32_t vgsi = 0; vgsi < vgs_count; ++vgsi)
    {
      Peregrine::Counter<true> m(dg.rbi, &dg, vgsi, cands);
      for (uint32_t v = 1; v <= num_vertices; ++v)
      {
        c += m.template map_into<L>(v);
      }
    }
  }
  else
  {
    for (uint32_t vgsi = 0; vgsi < vgs_count; ++vgsi)
    {
      Peregrine::Counter<false> m(dg.rbi, &dg, vgsi, cands);
      for (uint32_t v = 1; v <= num_vertices; ++v)
      {
        c += m.template map_into<L>(v);
      }
    }
  }

  return c;
}

template <typename DataGraphT>
void matchcount(DataGraphT &&data_graph,
          const std::vector<SmallGraph> &patterns,
          std::unordered_map<SmallGraph, uint64_t> &truth)
{
  DataGraph dg(data_graph);
  for (uint32_t i = 0; i < patterns.size(); ++i)
  {
    uint64_t count = 0;
    const auto process = [&count](const Peregrine::CompleteMatch &) -> void { count += 1; };
    const auto &p = patterns[i];
    std::cout << p.to_string() << std::flush;
    Graph::Labelling L = p.get_labelling();
    auto t1 = utils::get_timestamp();
    switch (L)
    {
      case Graph::LABELLED:
        matchPattern<Graph::LABELLED>(dg, p, process);
        break;
      default:
        matchPattern(dg, p, process);
    }
    auto t2 = utils::get_timestamp();

    CHECK_EQUAL(truth[p], count);
    std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
  }
}

template <typename DataGraphT>
void count(DataGraphT &&data_graph,
          const std::vector<SmallGraph> &patterns,
          std::unordered_map<SmallGraph, uint64_t> &truth)
{
  DataGraph dg(data_graph);
  for (uint32_t i = 0; i < patterns.size(); ++i)
  {
    uint64_t count = 0;
    const auto &p = patterns[i];
    std::cout << p.to_string() << std::flush;
    Graph::Labelling L = p.get_labelling();
    auto t1 = utils::get_timestamp();
    switch (L)
    {
      case Graph::LABELLED:
        count += countPattern<Graph::LABELLED>(dg, p);
        break;
      default:
        count += countPattern(dg, p);
    }
    auto t2 = utils::get_timestamp();

    CHECK_EQUAL(truth[p], count);
    std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
  }
}


SUITE(PatternMatcherTests)
{
  TEST(AntiVertices)
  {
    DataGraph dgs[] = {
      SmallGraph(cl3), 
      SmallGraph(twe),
      SmallGraph(swd),
      SmallGraph(cl4)
    };

    SmallGraph p(cl3);
    p.add_anti_edge(1, 4);
    
    {
      uint64_t truth[] = {3, 2, 2, 0};
      uint64_t count = 0;
      const auto process = [&count](auto &&) { count += 1; };

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        matchPattern(dgs[i], p, process);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    p.add_anti_edge(2, 4);
    {
      uint64_t truth[] = {3, 3, 4, 0};
      uint64_t count = 0;
      const auto process = [&count](auto &&) { count += 1; };

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        matchPattern(dgs[i], p, process);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    p.add_anti_edge(3, 4);
    {
      uint64_t truth[] = {1, 1, 2, 0};
      uint64_t count = 0;
      const auto process = [&count](auto &&) { count += 1; };

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        matchPattern(dgs[i], p, process);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    p.add_anti_edge(1, 5);
    {
      uint64_t truth[] = {3, 2, 2, 0};
      uint64_t count = 0;
      const auto process = [&count](auto &&) { count += 1; };

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        matchPattern(dgs[i], p, process);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    p.remove_edge(2, 4);
    p.remove_edge(3, 4);
    {
      // automorphic anti-vertices are redundant: it's the same as if just one
      // were there
      uint64_t truth[] = {3, 2, 2, 0};
      uint64_t count = 0;
      const auto process = [&count](auto &&) { count += 1; };

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        matchPattern(dgs[i], p, process);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    // different vertex ID's: anti-vertex as ID 1
    SmallGraph pinv;
    pinv.add_anti_edge(1, 2)
      .add_edge(2, 3)
      .add_edge(2, 4)
      .add_edge(3, 4);
    {
      uint64_t count = 0;
      const auto process = [&count](auto &&) { count += 1; };
      {
        std::cout << pinv.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        CHECK_THROW(matchPattern(dgs[0], pinv, process), std::runtime_error);
        auto t2 = utils::get_timestamp();

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }


    SmallGraph p2(s3); // different base topology
    p2.add_anti_edge(1, 4);
    {
      uint64_t truth[] = {3, 2, 2, 0};
      uint64_t count = 0;
      const auto process = [&count](auto &&) { count += 1; };

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p2.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        matchPattern(dgs[i], p2, process);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    p2.add_anti_edge(2, 3); // anti-edge as well! to check interplay
    {
      uint64_t truth[] = {0, 0, 0, 0};
      uint64_t count = 0;
      const auto process = [&count](auto &&) { count += 1; };

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p2.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        matchPattern(dgs[i], p2, process);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    SmallGraph dgs2p[] = {
      SmallGraph(c4),
      SmallGraph(c4).add_edge(1, 5),
      SmallGraph(c4).add_edge(4, 5),
      SmallGraph(c4).add_edge(1, 5).add_edge(2, 5),
      SmallGraph(c4).add_edge(1, 5).add_edge(4, 5)
    };

    DataGraph dgs2[] = {
      dgs2p[0],
      dgs2p[1],
      dgs2p[2],
      dgs2p[3],
      dgs2p[4]
    };

    SmallGraph p3(c4); // different base topology
    p3.add_anti_edge(4, 5);
    {
      uint64_t truth[] = {2, 4, 2, 6, 2};
      uint64_t count = 0;
      const auto process = [&count](auto &&) { count += 1; };

      for (uint32_t i = 0; i < 5; ++i)
      {
        std::cout << p3.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        matchPattern(dgs2[i], p3, process);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    // labelled anti-vertices
    SmallGraph p4(cl3);
    p4.add_anti_edge(1, 4)
      .set_label(1, 1)
      .set_label(2, 2)
      .set_label(3, 3)
      .set_label(4, 1);

    DataGraph dgs3[] = {
      SmallGraph(cl3, {1, 2, 3}),
      SmallGraph(twe, {1, 2, 3, 4}),
      SmallGraph(twe, {1, 2, 3, 1}),
      SmallGraph(swd, {1, 2, 3, 4}),
      SmallGraph(swd, {1, 2, 3, 1}),
      SmallGraph(cl4, {1, 2, 3, 4}),
      SmallGraph(cl4, {1, 2, 3, 1})
    };

    {
      uint64_t truth[] = {1, 1, 0, 1, 0, 1, 0};
      uint64_t count = 0;
      const auto process = [&count](const Peregrine::CompleteMatch &&) { count += 1; };

      for (uint32_t i = 0; i < 7; ++i)
      {
        std::cout << p4.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        matchPattern<Graph::LABELLED>(dgs3[i], p4, process);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    p4.set_label(3, 1);

    DataGraph dgs4[] = {
      SmallGraph(cl3, {1, 2, 1}),
      SmallGraph(twe, {1, 2, 1, 4}),
      SmallGraph(twe, {1, 2, 1, 1}),
      SmallGraph(swd, {1, 2, 1, 4}),
      SmallGraph(swd, {1, 2, 1, 1}),
      SmallGraph(cl4, {1, 2, 1, 4}),
      SmallGraph(cl4, {1, 2, 1, 1})
    };

    {
      uint64_t truth[] = {2, 2, 1, 2, 0, 2, 0};
      uint64_t count = 0;
      const auto process = [&count](const Peregrine::CompleteMatch &&) { count += 1; };

      for (uint32_t i = 0; i < 7; ++i)
      {
        std::cout << p4.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        matchPattern<Graph::LABELLED>(dgs4[i], p4, process);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }


  }

  TEST(AntiEdges)
  {
    // mico is too slow single-threaded
    DataGraph data_graph("data/citeseer");
    SmallGraph p(s4);

    // filter
    const auto filter1 = [&](uint64_t &c, uint32_t x, uint32_t y)
    {
      return [&c, &data_graph, x, y](const std::vector<uint32_t> &mapping)
      {
        assert(mapping.size() == 4);
        assert(x < mapping.size());
        assert(y < mapping.size());
        uint32_t v1 = mapping[x]; // vertex 2
        uint32_t v2 = mapping[y]; // vertex 3

        const adjlist &adj1 = data_graph.get_adj(v1);

        if (!std::binary_search(adj1.ptr, adj1.ptr + adj1.length, v2))
        {
          c += 1;
        }
      };
    };

    const auto filter2 = [&](uint64_t &c, uint32_t x, uint32_t y, uint32_t z)
    {
      return [&c, &data_graph, x, y, z](auto &&mapping)
        {
          uint32_t v1 = mapping[x]; // vertex 2
          uint32_t v2 = mapping[y]; // vertex 2
          uint32_t v3 = mapping[z]; // vertex 4

          const adjlist &adj1 = data_graph.get_adj(v1);

          if (!std::binary_search(adj1.ptr, adj1.ptr + adj1.length, v2)
              && !std::binary_search(adj1.ptr, adj1.ptr + adj1.length, v3))
          {
            c += 1;
          }
        };
    };

    const auto filter3 = [&](uint64_t &c)
    {
      return [&c, &data_graph](auto &&mapping)
      {
        uint32_t v1 = mapping[1]; // vertex 2
        uint32_t v2 = mapping[2]; // vertex 3
        uint32_t v3 = mapping[3]; // vertex 4

        const adjlist &adj1 = data_graph.get_adj(v1);
        const adjlist &adj2 = data_graph.get_adj(v2);

        if (!std::binary_search(adj1.ptr, adj1.ptr + adj1.length, v2)
            && !std::binary_search(adj1.ptr, adj1.ptr + adj1.length, v3)
            && !std::binary_search(adj2.ptr, adj2.ptr + adj2.length, v3))
        {
          c += 1;
        }
      };
    };

    for (uint32_t i = 2; i < 4; ++i)
    {
      for (uint32_t j = i+1; j <= 4; ++j)
      {
        SmallGraph antip(p);
        antip.add_anti_edge(i, j);
        std::cout << antip.to_string() << std::flush;

        uint64_t truth = 0;
        const auto filter = filter1(truth, i-1, j-1);
        matchPattern(data_graph, p, filter);

        CHECK(AnalyzedPattern(antip).has_anti_edges());

        uint64_t count = 0;
        const auto process = [&count](auto &&) { count += 1; };
        auto t1 = utils::get_timestamp();
        matchPattern(data_graph, antip, process);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth, count);
        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    for (uint32_t i = 2; i <= 4; ++i)
    {
      for (uint32_t j = 2; j <= 4; ++j)
      {
        if (j == i) continue;
        for (uint32_t k = j+1; k <= 4; ++k)
        {
          if (k == i) continue;
          SmallGraph antip(p);
          antip.add_anti_edge(i, j);
          antip.add_anti_edge(i, k);
          CHECK(AnalyzedPattern(antip).has_anti_edges());
          CHECK(!AnalyzedPattern(p).has_anti_edges());

          std::cout << antip.to_string() << std::flush;

          uint64_t truth = 0;
          const auto filter = filter2(truth, i-1, j-1, k-1);
          matchPattern(data_graph, p, filter);

          uint64_t count = 0;
          const auto process = [&count](auto &&) { count += 1; };
          auto t1 = utils::get_timestamp();
          matchPattern(data_graph, antip, process);
          auto t2 = utils::get_timestamp();

          CHECK_EQUAL(truth, count);
          std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
        }
      }
    }

    {
     SmallGraph antip(p);
     antip.add_anti_edge(2, 3);
     antip.add_anti_edge(2, 4);
     antip.add_anti_edge(3, 4);
     CHECK(!AnalyzedPattern(antip).has_anti_edges()); // vertex induced!
     CHECK(!AnalyzedPattern(p).has_anti_edges());

     std::cout << antip.to_string() << std::flush;

     uint64_t truth = 0;
     const auto filter = filter3(truth);
     matchPattern(data_graph, p, filter);

     uint64_t count = 0;
     const auto process = [&count](auto &&) { count += 1; };
     auto t1 = utils::get_timestamp();
     matchPattern(data_graph, antip, process);
     auto t2 = utils::get_timestamp();

     CHECK_EQUAL(truth, count);
     std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
    }
  }

  TEST(Unlabelled3MotifsNoAnti)
  {
    std::string data_graph= "data/citeseer";
    std::vector<SmallGraph> patterns = PatternGenerator::all(3,
        PatternGenerator::VERTEX_BASED,
        PatternGenerator::EXCLUDE_ANTI_EDGES);

    matchcount(data_graph, patterns, unlabelled_ground_truth_edge);

    std::cout << "DONE 3MOTIFS EDGE" << std::endl;
  }

  TEST(Unlabelled4MotifsNoAnti)
  {
    std::string data_graph= "data/citeseer";
    std::vector<SmallGraph> patterns = PatternGenerator::all(4,
        PatternGenerator::VERTEX_BASED,
        PatternGenerator::EXCLUDE_ANTI_EDGES);

    matchcount(data_graph, patterns, unlabelled_ground_truth_edge);
    std::cout << "DONE 4MOTIFS EDGE" << std::endl;
  }

  TEST(Unlabelled3Motifs)
  {
    std::string data_graph= "data/citeseer";
    std::vector<SmallGraph> patterns = PatternGenerator::all(3,
        PatternGenerator::VERTEX_BASED,
        PatternGenerator::INCLUDE_ANTI_EDGES);

    matchcount(data_graph, patterns, unlabelled_ground_truth_vtx);

    std::cout << "DONE 3MOTIFS" << std::endl;
  }

  TEST(Unlabelled4Motifs)
  {
    std::string data_graph= "data/citeseer";
    std::vector<SmallGraph> patterns = PatternGenerator::all(4,
        PatternGenerator::VERTEX_BASED,
        PatternGenerator::INCLUDE_ANTI_EDGES);

    matchcount(data_graph, patterns, unlabelled_ground_truth_vtx);
    std::cout << "DONE 4MOTIFS" << std::endl;
  }

  TEST(ToughPatterns)
  {
    std::string data_graph= "data/citeseer";

    std::vector<SmallGraph> patterns = {p6, labelling1, labelling2};

    matchcount(data_graph, patterns, unlabelled_ground_truth_edge);

    std::cout << "DONE TOUGH PATTERNS" << std::endl;
  }
}


SUITE(PatternCounterTests)
{
  TEST(AntiVertices)
  {
    DataGraph dgs[] = {
      SmallGraph(cl3), 
      SmallGraph(twe),
      SmallGraph(swd),
      SmallGraph(cl4)
    };

    SmallGraph p(cl3);
    p.add_anti_edge(1, 4);
    
    {
      uint64_t truth[] = {3, 2, 2, 0};

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        uint64_t count = countPattern(dgs[i], p);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    p.add_anti_edge(2, 4);
    {
      uint64_t truth[] = {3, 3, 4, 0};

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        uint64_t count = countPattern(dgs[i], p);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);
        count = 0;

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    p.add_anti_edge(3, 4);
    {
      uint64_t truth[] = {1, 1, 2, 0};

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        uint64_t count = countPattern(dgs[i], p);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    p.add_anti_edge(1, 5);
    {
      uint64_t truth[] = {3, 2, 2, 0};

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        uint64_t count = countPattern(dgs[i], p);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    p.remove_edge(2, 4);
    p.remove_edge(3, 4);
    {
      // automorphic anti-vertices are redundant: it's the same as if just one
      // were there
      uint64_t truth[] = {3, 2, 2, 0};

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        uint64_t count = countPattern(dgs[i], p);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }


    SmallGraph p2(s3); // different base topology
    p2.add_anti_edge(1, 4);
    {
      uint64_t truth[] = {3, 2, 2, 0};

      for (uint32_t i = 0; i < 4; ++i)
      {
        std::cout << p2.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        uint64_t count = countPattern(dgs[i], p2);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }

    SmallGraph dgs2p[] = {
      SmallGraph(c4),
      SmallGraph(c4).add_edge(1, 5),
      SmallGraph(c4).add_edge(4, 5),
      SmallGraph(c4).add_edge(1, 5).add_edge(2, 5),
      SmallGraph(c4).add_edge(1, 5).add_edge(4, 5)
    };

    DataGraph dgs2[] = {
      dgs2p[0],
      dgs2p[1],
      dgs2p[2],
      dgs2p[3],
      dgs2p[4]
    };

    SmallGraph p3(c4); // different base topology
    p3.add_anti_edge(4, 5);
    {
      uint64_t truth[] = {2, 4, 2, 6, 2};

      for (uint32_t i = 0; i < 5; ++i)
      {
        std::cout << p3.to_string() << std::flush;
        auto t1 = utils::get_timestamp();
        uint64_t count = countPattern(dgs2[i], p3);
        auto t2 = utils::get_timestamp();

        CHECK_EQUAL(truth[i], count);

        std::cout << "\t✓ " << (t2-t1)/1e6 << "s" << std::endl;
      }
    }
  }

  TEST(Unlabelled3MotifsNoAnti)
  {
    std::string data_graph= "data/citeseer";
    std::vector<SmallGraph> patterns = PatternGenerator::all(3,
        PatternGenerator::VERTEX_BASED,
        PatternGenerator::EXCLUDE_ANTI_EDGES);

    count(data_graph, patterns, unlabelled_ground_truth_edge);

    std::cout << "DONE 3MOTIFS EDGE" << std::endl;
  }

  TEST(Unlabelled4MotifsNoAnti)
  {
    std::string data_graph= "data/citeseer";
    std::vector<SmallGraph> patterns = PatternGenerator::all(4,
        PatternGenerator::VERTEX_BASED,
        PatternGenerator::EXCLUDE_ANTI_EDGES);

    count(data_graph, patterns, unlabelled_ground_truth_edge);

    std::cout << "DONE 4MOTIFS EDGE" << std::endl;
  }

  TEST(Unlabelled3Motifs)
  {
    std::string data_graph= "data/citeseer";
    std::vector<SmallGraph> patterns = PatternGenerator::all(3,
        PatternGenerator::VERTEX_BASED,
        PatternGenerator::INCLUDE_ANTI_EDGES);

    count(data_graph, patterns, unlabelled_ground_truth_vtx);

    std::cout << "DONE 3MOTIFS" << std::endl;
  }

  TEST(Unlabelled4Motifs)
  {
    std::string data_graph= "data/citeseer";
    std::vector<SmallGraph> patterns = PatternGenerator::all(4,
        PatternGenerator::VERTEX_BASED,
        PatternGenerator::INCLUDE_ANTI_EDGES);

    count(data_graph, patterns, unlabelled_ground_truth_vtx);

    std::cout << "DONE 4MOTIFS" << std::endl;
  }

  TEST(ToughPatterns)
  {
    std::string data_graph= "data/citeseer";
    DataGraph dg2(p6);
    std::vector<SmallGraph> patterns = {p6, labelling1, labelling2};

    count(data_graph, patterns, unlabelled_ground_truth_edge);
    CHECK_EQUAL(1, countPattern(dg2, {p6}));
    std::cout << "DONE TOUGH PATTERNS" << std::endl;
  }
}
