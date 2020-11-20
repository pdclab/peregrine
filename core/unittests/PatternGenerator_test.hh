#include <iostream>

#include "../Graph.hh"
#include "../PatternGenerator.hh"

using namespace Peregrine;
SUITE(PatternGeneratorTests)
{
  SUITE(VertexExtension)
  {
    TEST(StarUnlabelled)
    {
      auto s = PatternGenerator::star(3);
      auto res = PatternGenerator::extend({s}, PatternGenerator::VERTEX_BASED);
      // 0 distinct labels (n), 2 indsets with 1 (m1) and 2 (m2) elements,
      // (n+1)*((m1+1)*(m2+1) - 1) = 1*2*3 - 1 = 5
      // the minus 1 is for the vertex extension with no added edges that we discount
      CHECK_EQUAL(5, res.size());
    }

    TEST(StarSameLabels)
    {
      auto s = PatternGenerator::star(3);
      s.set_label(1, 1)
       .set_label(2, 1)
       .set_label(3, 1);
      auto res = PatternGenerator::extend({s}, PatternGenerator::VERTEX_BASED);
      // 1 distinct label (n), 2 indsets with 1 (m1) and 2 (m2) elements,
      // (n+1)*((m1+1)*(m2+1) - 1) = 2(2*3 - 1) = 10
      CHECK_EQUAL(10, res.size());
    }

    TEST(StarDifferentEndLabels)
    {
      auto s = PatternGenerator::star(3);
      s.set_label(1, 1)
       .set_label(2, 2)
       .set_label(3, 1);
      auto res = PatternGenerator::extend({s}, PatternGenerator::VERTEX_BASED);
      // 2 distinct labels (n), 3 indsets with 1 (m1), 1 (m2), and 1 (m3) elements,
      // (n+1)*((m1+1)*(m2+1)*(m3+1) - 1) = 3(2*2*2 - 1) = 21
      CHECK_EQUAL(21, res.size());
    }
  }

  SUITE(ChainEdgeExtension)
  {
    TEST(ChainUnlabelled)
    {
      auto s = SmallGraph(std::vector<std::pair<uint32_t,uint32_t>>{{1,2},{2,3},{3,4}});
      auto res = PatternGenerator::extend({s}, PatternGenerator::EDGE_BASED);
      // 0 labels (n), 2 indsets (m), and 2 asymmetric pairs of disconnected vertices (d)
      // (n+1)*m + d = 2 + 2 = 4
      CHECK_EQUAL(4, res.size());
    }
  }

  SUITE(StarEdgeExtension)
  {
    TEST(Unlabelled)
    {
      auto s = PatternGenerator::star(3);
      auto res = PatternGenerator::extend({s}, PatternGenerator::EDGE_BASED);
      // 0 labels (n), 2 indsets (m), and 1 pair of disconnected vertices (d)
      // (n+1)*m + d = 1*2 + 1 = 3
      CHECK_EQUAL(3, res.size()); // triangle, 4-chain, 4-star
    }

    SUITE(Labelled)
    {
      TEST(SameLabels)
      {
        auto s = PatternGenerator::star(3);
        s.set_label(1, 1)
         .set_label(2, 1)
         .set_label(3, 1);
        auto res = PatternGenerator::extend({s}, PatternGenerator::EDGE_BASED);
        // fully-labelled triangle, 4-chain, and 4-star, then
        // a partially-labelled 4-chain and partially-labelled 4-star
        CHECK_EQUAL(5, res.size());
      }

      TEST(DifferentCentreLabel)
      {
        auto s = PatternGenerator::star(3);
        s.set_label(1, 2)
         .set_label(2, 1)
         .set_label(3, 1);
        auto res = PatternGenerator::extend({s}, PatternGenerator::EDGE_BASED);
        // 2 distinct labels (n), 2 indsets (m), 1 pair of disconnected vertices (d)
        // so (n+1)*m + d = 3*2+1 = 7
        CHECK_EQUAL(7, res.size());
      }

      TEST(DifferentEndLabels)
      {
        auto s = PatternGenerator::star(3);
        s.set_label(1, 1)
         .set_label(2, 1)
         .set_label(3, 2);
        auto res = PatternGenerator::extend({s}, PatternGenerator::EDGE_BASED);
        // fully-labelled triangle, 4 4-chains, and 2 4-stars, then
        // 2 partially-labelled 4-chains and 1 partially-labelled 4-star
        // another way: 2 distinct labels (n), 3 indsets (m), 1 pair of
        // disconnected vertices (d) so (n+1)*m + d = 3*3+1 = 10
        CHECK_EQUAL(10, res.size());
      }

      TEST(AllDifferentLabels)
      {
        auto s = PatternGenerator::star(3);
        s.set_label(1, 1)
         .set_label(2, 2)
         .set_label(3, 3);
        auto res = PatternGenerator::extend({s}, PatternGenerator::EDGE_BASED);
        // fully-labelled triangle, 6 4-chain, and 3 4-star, then
        // 2 partially-labelled 4-chain and 1 partially-labelled 4-star
        // another way: 3 labels, 3 indsets, 1 pair of disconnected vertices,
        // so (n+1)*m + d = 4*3 + 1 = 13
        CHECK_EQUAL(13, res.size());
      }
    }
  }

  SUITE(LabelledMultiplePatternEdgeExtension)
  {
    TEST(StarsIso)
    {
      auto s1 = PatternGenerator::star(3);
      auto s2 = PatternGenerator::star(3);
      s1.set_label(1, 1)
        .set_label(2, 2)
        .set_label(3, 1);
      s2.set_label(1, 1)
        .set_label(2, 1)
        .set_label(3, 2);
      auto res = PatternGenerator::extend({s1, s2}, PatternGenerator::EDGE_BASED);
      // s1 and s2 are isomorphic, so their extensions will also be isomorphic
      CHECK_EQUAL(10, res.size());
    }

    TEST(StarsOverlap)
    {
      auto s1 = PatternGenerator::star(3);
      auto s2 = PatternGenerator::star(3);
      s1.set_label(1, 1)
        .set_label(2, 2)
        .set_label(3, 1);
      s2.set_label(1, 2)
        .set_label(2, 1)
        .set_label(3, 1);
      auto res = PatternGenerator::extend({s1, s2}, PatternGenerator::EDGE_BASED);
      // s1 and s2 are not isomorphic, but at least two extensions of s1 are also
      // extensions of s2:
      //  * the 4-chain with one centre vertex labeled 2 and
      //    the other labeled 1
      //  * the triangle with one vertex labeled 2 and the others labeled 1
      CHECK_EQUAL(15, res.size());
    }
  }

  TEST(AntiVerticesEdgeBased)
  {
    auto s1 = PatternGenerator::star(3);
    s1.add_anti_edge(1, 4);

    auto res1 = PatternGenerator::extend({s1}, PatternGenerator::EDGE_BASED,
        PatternGenerator::MAINTAIN_ANTI_EDGES);

    auto res2 = PatternGenerator::extend({s1}, PatternGenerator::EDGE_BASED,
        PatternGenerator::OVERWRITE_ANTI_EDGES);

    CHECK(res1 == res2);
    for (const auto &p : res1)
    {
      if (p.num_vertices() == 4)
      {
        CHECK(p.get_anti_neighbours(1) == std::vector<uint32_t>{5});
        CHECK(p.get_anti_neighbours(5) == std::vector<uint32_t>{1});
        CHECK(p.get_anti_neighbours(4).empty());
        CHECK(!p.get_neighbours(4).empty());
        CHECK(p.num_vertices() < 5);
      }
      else
      {
        CHECK(p.get_anti_neighbours(1) == std::vector<uint32_t>{4});
        CHECK(p.get_anti_neighbours(4) == std::vector<uint32_t>{1});
        CHECK(p.num_vertices() < 4);
      }
    }

    s1.add_anti_edge(2, 4);
    auto res6 = PatternGenerator::extend({s1}, PatternGenerator::EDGE_BASED);
    for (const auto &p : res6)
    {
      if (p.num_vertices() == 4)
      {
        CHECK((p.get_anti_neighbours(5) == std::vector<uint32_t>{1, 2}));
        CHECK(p.get_anti_neighbours(1) == std::vector<uint32_t>{5});
        CHECK(p.get_anti_neighbours(2) == std::vector<uint32_t>{5});
        CHECK(p.get_anti_neighbours(4).empty());
        CHECK(!p.get_neighbours(4).empty());
        CHECK(p.num_vertices() < 5);
      }
      else
      {
        CHECK((p.get_anti_neighbours(4) == std::vector<uint32_t>{1, 2}));
        CHECK(p.get_anti_neighbours(1) == std::vector<uint32_t>{4});
        CHECK(p.get_anti_neighbours(2) == std::vector<uint32_t>{4});
        CHECK(p.num_vertices() < 4);
      }
    }
  }

  TEST(AntiVerticesVertBased)
  {
    auto s1 = PatternGenerator::star(3);
    s1.add_anti_edge(1, 4);

    auto res3 = PatternGenerator::extend({s1}, PatternGenerator::VERTEX_BASED,
        PatternGenerator::MAINTAIN_ANTI_EDGES);

    auto res4 = PatternGenerator::extend({s1}, PatternGenerator::VERTEX_BASED,
        PatternGenerator::OVERWRITE_ANTI_EDGES);

    // equal true edges, but overwrite should add an anti-edge between 2,3 in res4
    CHECK(res3 == res4);
    for (const auto &p : res3)
    {
      CHECK(p.get_anti_neighbours(1) == std::vector<uint32_t>{5});
      // 1 is still the centre
      CHECK((p.get_neighbours(1).size() > 2) || (p.get_neighbours(1) == std::vector<uint32_t>{2, 3}));
      CHECK(p.get_anti_neighbours(5) == std::vector<uint32_t>{1});
      CHECK(p.get_anti_neighbours(4).empty());
      CHECK(!p.get_neighbours(4).empty());
      CHECK(p.num_vertices() < 5);
    }

    s1.add_anti_edge(2, 4);
    auto res5 = PatternGenerator::extend({s1}, PatternGenerator::VERTEX_BASED, PatternGenerator::MAINTAIN_ANTI_EDGES);
    for (const auto &p : res5)
    {
      CHECK((p.get_anti_neighbours(5) == std::vector<uint32_t>{1, 2}));
      CHECK(p.get_anti_neighbours(1) == std::vector<uint32_t>{5});
      CHECK(p.get_anti_neighbours(2) == std::vector<uint32_t>{5});
      CHECK(p.get_anti_neighbours(4).empty());
      CHECK(!p.get_neighbours(4).empty());
      CHECK(p.num_vertices() == 4);
    }

    auto res6 = PatternGenerator::extend({s1}, PatternGenerator::VERTEX_BASED, PatternGenerator::OVERWRITE_ANTI_EDGES);
    for (const auto &p : res6)
    {
      CHECK((p.get_anti_neighbours(5) == std::vector<uint32_t>{1, 2}));
      CHECK(p.get_anti_neighbours(1).size() >= 1);
      CHECK(utils::search(p.get_anti_neighbours(1), 5u));
      CHECK(p.get_anti_neighbours(2).size() > 1);
      CHECK(utils::search(p.get_anti_neighbours(2), 5u));
      CHECK(p.get_anti_neighbours(4).size() < 4);
      CHECK(!utils::search(p.get_anti_neighbours(4), 5u));
      CHECK(!p.get_neighbours(4).empty());
      CHECK(p.num_vertices() == 4);
    }
  }

  TEST(AntiEdgesVertBased)
  {
    auto s1 = PatternGenerator::star(3);
    s1.add_anti_edge(2, 3);

    auto res1 = PatternGenerator::extend({s1}, PatternGenerator::VERTEX_BASED,
        PatternGenerator::OVERWRITE_ANTI_EDGES);

    auto res2 = PatternGenerator::extend({s1}, PatternGenerator::VERTEX_BASED,
        PatternGenerator::MAINTAIN_ANTI_EDGES);

    auto is_vinduced = [](const SmallGraph &p) {
        uint32_t m = p.num_anti_edges() + p.num_true_edges();
        uint32_t n = p.num_vertices();
        return m == (n*(n-1))/2;
      };
    
    // vertex extension doesn't add edges between existing vertices and s1 is
    // vinduced, so there should be no difference between
    // overwriting/maintaining
    CHECK(res1 == res2);
    for (const auto &p : res1)
    {
      CHECK(is_vinduced(p));
    }
  }

  TEST(AntiEdgesEdgeBased)
  {
    auto s1 = PatternGenerator::star(3);
    s1.add_anti_edge(2, 3);

    auto res3 = PatternGenerator::extend({s1}, PatternGenerator::EDGE_BASED,
        PatternGenerator::OVERWRITE_ANTI_EDGES);

    auto res4 = PatternGenerator::extend({s1}, PatternGenerator::EDGE_BASED,
        PatternGenerator::MAINTAIN_ANTI_EDGES);

    // res4 should be missing the triangle
    CHECK_EQUAL(res4.size() + 1, res3.size());
    for (const auto &p : res4)
    {
      CHECK_EQUAL(4, p.num_vertices());
      CHECK(p.get_anti_neighbours(2) == std::vector<uint32_t>{3});
      CHECK(p.get_anti_neighbours(3) == std::vector<uint32_t>{2});
    }

    for (const auto &p : res3)
    {
      if (p.num_vertices() == 4)
      {
        CHECK(p.get_anti_neighbours(2) == std::vector<uint32_t>{3});
        CHECK(p.get_anti_neighbours(3) == std::vector<uint32_t>{2});
      }
      else
      {
        CHECK(p.get_anti_neighbours(2).empty());
        CHECK(p.get_anti_neighbours(3).empty());
        CHECK((p.get_neighbours(2) == std::vector<uint32_t>{1, 3}));
        CHECK((p.get_neighbours(3) == std::vector<uint32_t>{1, 2}));
      }
    }

  }
}
