#include "../Peregrine.hh"

std::vector<uint64_t> vtx_counts = {
  111153, // [2-3][1-2][1-4](3~4)(2~4)(1~3)
  222630, // [1-2][1-3][1-4](3~4)(2~3)(2~4)
  22900,  // [2-3][1-2][1-3][1-4](3~4)(2~4)
  3094,   // [3-4][2-3][1-2][1-4](2~4)(1~3)
  2200,   // [3-4][2-3][1-2][1-3][1-4](2~4)
  255,    // [3-4][2-3][2-4][1-2][1-3][1-4]
};

std::vector<uint64_t> edge_counts = {
  185589, // [2-3][1-2][1-4]
  250950, // [1-2][1-3][1-4]
  34760,  // [2-3][1-2][1-3][1-4]
  6059,   // [3-4][2-3][1-2][1-4]
  3730,   // [3-4][2-3][1-2][1-3][1-4]
  255,    // [3-4][2-3][2-4][1-2][1-3][1-4]
};

TEST(CountingIntegration)
{
  // the pattern matching unit tests check correctness for counting, this is
  // just to make sure nothing gets messed up at the boundary of the Counter

  const std::string data_graph_name("data/citeseer");
  size_t nthreads = std::thread::hardware_concurrency();

  std::cout << "RUNNING COUNTING INTEGRATION TEST" << std::endl;

  auto t1 = utils::get_timestamp();

  uint32_t k = 4;
  {
    // include anti-edges
    std::vector<Peregrine::SmallGraph> patterns = Peregrine::PatternGenerator::all(k,
        Peregrine::PatternGenerator::VERTEX_BASED,
        Peregrine::PatternGenerator::INCLUDE_ANTI_EDGES);

    auto res = Peregrine::count(data_graph_name, patterns, nthreads);
    int i = 0;
    for (auto &[p, count] : res)
    {
      CHECK_EQUAL(vtx_counts[i++], count);
    }
  }

  {
    // exclude anti-edges
    std::vector<Peregrine::SmallGraph> patterns = Peregrine::PatternGenerator::all(k,
        Peregrine::PatternGenerator::VERTEX_BASED,
        Peregrine::PatternGenerator::EXCLUDE_ANTI_EDGES);

    auto res = Peregrine::count(data_graph_name, patterns, nthreads);
    int i = 0;
    for (auto &[p, count] : res)
    {
      CHECK_EQUAL(edge_counts[i++], count);
    }
  }

  auto t2 = utils::get_timestamp();

  std::cout << "COUNTING INTEGRATION TEST FINISHED IN " << (t2-t1)/1e6 << "s" << std::endl;
}
