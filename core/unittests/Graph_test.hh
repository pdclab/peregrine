#include "../Graph.hh"
#include <iostream>
#include <vector>
#include "../utils.hh"

using namespace Peregrine;
SUITE(GraphPatternAnalysis)
{
  TEST(AntiVertices)
  {
    SmallGraph s(std::vector<std::pair<uint32_t,uint32_t>>{{1,2}, {1,3}, {2, 3}});
    AnalyzedPattern r1(s);

    s.add_anti_edge(1, 4);
    s.add_anti_edge(2, 4);
    s.add_anti_edge(3, 4);

    AnalyzedPattern r2(s);

    // shouldn't affect core pattern
    CHECK_EQUAL(r1.vgs.size(), r2.vgs.size());

    CHECK_EQUAL(r1.vgs.front().bliss_hash(), r2.vgs.front().bliss_hash());
  }

  TEST(AntiEdges)
  {
    SmallGraph s(std::vector<std::pair<uint32_t,uint32_t>>{{1,2}, {1,3}, {1, 4}});
    AnalyzedPattern r1(s);
    CHECK(!r1.has_anti_edges());

    s.add_anti_edge(2, 5);
    AnalyzedPattern r2(s);
    CHECK(!r2.has_anti_edges()); // only anti-edges are for anti-vertices

    s.add_edge(1, 5);
    AnalyzedPattern r3(s);
    CHECK(r3.has_anti_edges()); // no longer has an anti-vertex

    s.add_anti_edge(2, 3);
    s.add_anti_edge(2, 4);
    s.add_anti_edge(3, 4);
    s.add_anti_edge(3, 5);
    s.add_anti_edge(4, 5);

    AnalyzedPattern r4(s);
    CHECK(!r4.has_anti_edges()); // vertex induced
  }

  TEST(Conditions)
  {
    SmallGraph s1(std::vector<std::pair<uint32_t,uint32_t>>{{1,2}, {1,4}, {2, 3}});
    SmallGraph s2(std::vector<std::pair<uint32_t,uint32_t>>{{1,2}, {2,3}, {3, 4}});

    AnalyzedPattern r1(s1);
    AnalyzedPattern r2(s2);

    CHECK_EQUAL(r1.conditions.size(), r2.conditions.size());

    CHECK_EQUAL(r1.conditions[0].first,  1);
    CHECK_EQUAL(r1.conditions[0].second, 2);

    CHECK_EQUAL(r2.conditions[0].first,  2);
    CHECK_EQUAL(r2.conditions[0].second, 3);
  }

  TEST(IndSetsUnlabelled)
  {
    // 4 chain, should have indsets {1, 3} and {2, 4}
    SmallGraph s(std::vector<std::pair<uint32_t,uint32_t>>{{1,2}, {2,3}, {3,4}});
    AnalyzedPattern r(s);
    CHECK_EQUAL(3, r.indsets.size());
  }

  TEST(IndSetsLabelled)
  {
    // 3-star with differently labeled endpoints
    SmallGraph s(std::vector<std::pair<uint32_t,uint32_t>>{{1,2}, {1,3}});
    s.set_label(1, 1)
     .set_label(2, 2)
     .set_label(3, 3);
    AnalyzedPattern r(s);

    CHECK_EQUAL(3, r.indsets.size());
    CHECK_EQUAL(1, r.indsets[0].size());
    CHECK_EQUAL(1, r.indsets[1].size());
    CHECK_EQUAL(1, r.indsets[2].size());
  }

  TEST(NonRedOrder)
  {
    // chordal square
    SmallGraph s(std::vector<std::pair<uint32_t,uint32_t>> {
        {1,2}, {1,3}, {1,4}, {2,3}, {3,4}
    });
    AnalyzedPattern r(s);
    CHECK((r.order_groups.front() == std::vector<uint32_t>{2,4}));

    // with labels
    s.set_label(1, 1)
     .set_label(2, 1)
     .set_label(3, 3)
     .set_label(4, 3);
    AnalyzedPattern r2(s);
    CHECK((r2.order_groups.front() == std::vector<uint32_t>{2}));
    CHECK((r2.order_groups.back() == std::vector<uint32_t>{4}));

    // 5-star + edge
    SmallGraph s2(std::vector<std::pair<uint32_t,uint32_t>> {
        {1,2}, {1,3}, {1,4}, {1,5}, {2,3}
    });

    AnalyzedPattern r3(s2);
    CHECK((r3.order_groups.front() == std::vector<uint32_t>{3}));
    CHECK((r3.order_groups.back() == std::vector<uint32_t>{4,5}));

    // 4-star
    SmallGraph s3(std::vector<std::pair<uint32_t,uint32_t>> {
        {1,2}, {1,3}, {1,4}
    });
    AnalyzedPattern r4(s3);
    CHECK((r4.order_groups.front() == std::vector<uint32_t>{2, 4, 3}));
    CHECK(r4.sibling_groups.front() == 2);
  }

  TEST(GraphBuilder)
  {
    SmallGraph p;

    CHECK_EQUAL(0, p.num_vertices());
    CHECK_EQUAL(0, p.num_true_edges());
    CHECK_EQUAL(0, p.num_anti_edges());
    CHECK_EQUAL(Graph::UNLABELLED, p.get_labelling());

    p.add_edge(1, 4);
    CHECK_EQUAL(2, p.num_vertices());
    CHECK_EQUAL(1, p.num_true_edges());
    CHECK_EQUAL(0, p.num_anti_edges());
    CHECK_EQUAL(Graph::UNLABELLED, p.get_labelling());

    p.set_label(4, 'a');
    CHECK_EQUAL(2, p.num_vertices());
    CHECK_EQUAL(Graph::LABELLED, p.get_labelling());
    CHECK_EQUAL('a', p.label(4));
    CHECK_EQUAL(0, p.label(1));

    // Since the vertex 4 exists, we assume there will be vertices 2 and 3 at some point
    p.set_label(2, 'b');
    CHECK_EQUAL(2, p.num_vertices());
    CHECK_EQUAL(1, p.num_true_edges());
    CHECK_EQUAL(0, p.num_anti_edges());
    CHECK_EQUAL('a', p.label(4));
    CHECK_EQUAL('b', p.label(2));

    p.add_edge(2, 1);
    CHECK_EQUAL(3, p.num_vertices());
    CHECK_EQUAL(2, p.num_true_edges());
    CHECK_EQUAL(0, p.num_anti_edges());
    CHECK_EQUAL('a', p.label(4));
    CHECK_EQUAL('b', p.label(2));
    CHECK_EQUAL(0, p.label(1));

    p.set_label(1, static_cast<uint32_t>(-1));
    CHECK_EQUAL(Graph::PARTIALLY_LABELLED, p.get_labelling());
    CHECK_EQUAL('a', p.label(4));
    CHECK_EQUAL('b', p.label(2));
    CHECK_EQUAL(-1, p.label(1));
  }

  TEST(Automorphisms)
  {
    SmallGraph p;
    p.add_edge(1, 2)
     .add_edge(1, 4)
     .add_edge(2, 3)
     .add_edge(3, 4)
     .set_label(1, 0)
     .set_label(2, 0)
     .set_label(3, 4)
     .set_label(4, 4);

    SmallGraph pr; // rotated clockwise
    pr.add_edge(1, 2)
      .add_edge(1, 4)
      .add_edge(2, 3)
      .add_edge(3, 4)
      .set_label(1, 4)
      .set_label(2, 0)
      .set_label(3, 0)
      .set_label(4, 4);

    SmallGraph prr; // rotated clockwise twice
    prr.add_edge(1, 2)
      .add_edge(1, 4)
      .add_edge(2, 3)
      .add_edge(3, 4)
      .set_label(1, 4)
      .set_label(2, 4)
      .set_label(3, 0)
      .set_label(4, 0);

    CHECK_EQUAL(p, pr);
    CHECK_EQUAL(p, prr);
    CHECK_EQUAL(p.bliss_hash(), pr.bliss_hash());
    CHECK_EQUAL(p.bliss_hash(), prr.bliss_hash());

    AnalyzedPattern ap(p);
    AnalyzedPattern apr(pr);
    AnalyzedPattern aprr(prr);

    CHECK_EQUAL(ap.conditions.size(), apr.conditions.size());
    CHECK_EQUAL(ap.conditions.size(), aprr.conditions.size());
  }
}
