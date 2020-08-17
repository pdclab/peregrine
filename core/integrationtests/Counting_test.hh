#include "../Peregrine.hh"

TEST(CountingIntegration)
{
  // the pattern matching unit tests check correctness for counting, this is
  // just to make sure nothing gets messed up at the boundary of the Counter
  const std::string data_graph_name("data/citeseer");
  size_t nthreads = std::thread::hardware_concurrency();

  std::cout << "RUNNING COUNTING INTEGRATION TEST" << std::endl;

  auto t1 = utils::get_timestamp();

  for (unsigned k : {4, 5, 6})
  {
    {
      std::unordered_map<Peregrine::SmallGraph, uint64_t> truth;
      {
        std::ifstream truth_file("core/integrationtests/truth/" + std::to_string(k) + "m.txt");
        CHECK(truth_file.is_open());

        // Using the graph string representation may be expensive, but the test is
        // fast _enough_ and this makes the truth file human readable.
        std::string id;
        uint64_t count;
        while (truth_file >> id >> count)
        {
          Peregrine::SmallGraph p = from_string(id);
          truth[p] = count;
        }
      }

      std::cout << k << "-MOTIFS VERTEX-INDUCED" << std::flush;
      auto tt1 = utils::get_timestamp();
      // include anti-edges
      std::vector<Peregrine::SmallGraph> patterns = Peregrine::PatternGenerator::all(k,
          Peregrine::PatternGenerator::VERTEX_BASED,
          Peregrine::PatternGenerator::INCLUDE_ANTI_EDGES);

      auto res = Peregrine::count(data_graph_name, patterns, nthreads);
      for (auto &[p, count] : res)
      {
        CHECK_EQUAL(truth.at(p), count);
      }
      auto tt2 = utils::get_timestamp();
      std::cout << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
    }


    {
      std::unordered_map<Peregrine::SmallGraph, uint64_t> truth;
      {
        std::ifstream truth_file("core/integrationtests/truth/" + std::to_string(k) + "m-edge.txt");
        CHECK(truth_file.is_open());

        // Using the graph string representation may be expensive, but the test is
        // fast _enough_ and this makes the truth file human readable.
        std::string id;
        uint64_t count;
        while (truth_file >> id >> count)
        {
          Peregrine::SmallGraph p = from_string(id);
          truth[p] = count;
        }
      }

      std::cout << k << "-MOTIFS EDGE-INDUCED" << std::flush;
      auto tt1 = utils::get_timestamp();

      // exclude anti-edges
      std::vector<Peregrine::SmallGraph> patterns = Peregrine::PatternGenerator::all(k,
          Peregrine::PatternGenerator::VERTEX_BASED,
          Peregrine::PatternGenerator::EXCLUDE_ANTI_EDGES);

      auto res = Peregrine::count(data_graph_name, patterns, nthreads);
      for (auto &[p, count] : res)
      {
        CHECK_EQUAL(truth.at(p), count);
      }

      auto tt2 = utils::get_timestamp();
      std::cout << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
    }
  }

  auto t2 = utils::get_timestamp();

  std::cout << "COUNTING INTEGRATION TEST FINISHED IN " << (t2-t1)/1e6 << "s" << std::endl;
}
