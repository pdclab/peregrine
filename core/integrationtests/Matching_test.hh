#include "../Peregrine.hh"
#include "../../apps/Domain.hh"


TEST(MatchingIntegration)
{
  const std::string data_graph_name("data/citeseer");
  uint32_t k = 4;
  uint64_t threshold = 1;
  size_t nthreads = std::thread::hardware_concurrency();

  std::cout << "RUNNING MATCHING INTEGRATION TEST" << std::endl;

  std::vector<std::unordered_map<Peregrine::SmallGraph, uint64_t>> truth(k-1);
  {
    std::ifstream truth_file("core/integrationtests/truth/cs-supports.txt");
    CHECK(truth_file.is_open());

    // Using the graph string representation may be expensive, but the test is
    // fast _enough_ and this makes the truth file human readable.
    std::string id;
    uint64_t supp;
    while (truth_file >> id >> supp)
    {
      Peregrine::SmallGraph p = from_string_labelled(id);
      truth[p.num_true_edges()-2][p] = supp;
    }
  }

  auto t1 = utils::get_timestamp();

  std::vector<uint64_t> supports;
  std::vector<Peregrine::SmallGraph> freq_patterns;

  const auto view = [](auto &&v) { return v.get_support(); };

  // initial discovery
  {
    auto tt1 = utils::get_timestamp();
    const auto process = [](auto &&a, auto &&cm) {
      uint32_t merge = cm.pattern[1] == cm.pattern[2] ? 1 : 2;
      a.map(cm.pattern, std::make_pair(cm.mapping, merge));
    };

    std::vector<Peregrine::SmallGraph> patterns = {Peregrine::PatternGenerator::star(3)};
    patterns.front().set_labelling(Peregrine::Graph::DISCOVER_LABELS);
    auto psupps = Peregrine::match<Peregrine::Pattern, DiscoveryDomain<2>, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE>(data_graph_name, patterns, nthreads, process, view);
    for (const auto &[p, supp] : psupps)
    {
      if (supp >= threshold)
      {
        freq_patterns.push_back(p);
        supports.push_back(supp);
        REQUIRE CHECK_EQUAL(supp, truth[0].at(p));
      }
    }
    REQUIRE CHECK_EQUAL(truth[0].size(), freq_patterns.size());

    auto tt2 = utils::get_timestamp();
    std::cout << "Step 1\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
  }

  std::vector<Peregrine::SmallGraph> patterns = Peregrine::PatternGenerator::extend(freq_patterns, Peregrine::PatternGenerator::EDGE_BASED);

  const auto process = [](auto &&a, auto &&cm) {
    a.map(cm.pattern, cm.mapping);
  };

  uint32_t step = 2;
  while (step < k && !patterns.empty())
  {
    freq_patterns.clear();
    supports.clear();

    auto tt1 = utils::get_timestamp();
    auto psupps = Peregrine::match<Peregrine::Pattern, Domain, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE>(data_graph_name, patterns, nthreads, process, view);

    for (const auto &[p, supp] : psupps)
    {
      if (supp >= threshold)
      {
        freq_patterns.push_back(p);
        supports.push_back(supp);
        REQUIRE CHECK_EQUAL(supp, truth[step-1].at(p));
      }
    }
    REQUIRE CHECK_EQUAL(truth[step-1].size(), freq_patterns.size());

    auto tt2 = utils::get_timestamp();
    std::cout << "Step " << step << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;

    patterns = Peregrine::PatternGenerator::extend(freq_patterns, Peregrine::PatternGenerator::EDGE_BASED);

    step += 1;
  }
  auto t2 = utils::get_timestamp();

  std::cout << "MATCHING INTEGRATION TEST FINISHED IN " << (t2-t1)/1e6 << "s" << std::endl;
}

TEST(EdgeLabelDiscoveryIntegration)
{
  const std::string data_graph_name("data/citeseer");
  uint32_t k = 4;
  uint64_t threshold = 1;
  size_t nthreads = std::thread::hardware_concurrency();

  std::cout << "RUNNING MATCHING INTEGRATION TEST" << std::endl;

  std::vector<std::unordered_map<Peregrine::SmallGraph, uint64_t>> truth(k-1);
  {
    std::ifstream truth_file("core/integrationtests/truth/cs-supports.txt");
    CHECK(truth_file.is_open());

    // Using the graph string representation may be expensive, but the test is
    // fast _enough_ and this makes the truth file human readable.
    std::string id;
    uint64_t supp;
    while (truth_file >> id >> supp)
    {
      Peregrine::SmallGraph p = from_string_labelled(id);
      truth[p.num_true_edges()-2][p] = supp;
    }
  }

  auto t1 = utils::get_timestamp();

  std::vector<uint64_t> supports;
  std::vector<Peregrine::SmallGraph> freq_patterns;

  const auto view = [](auto &&v) { return v.get_support(); };

  // initial discovery
  {
    auto tt1 = utils::get_timestamp();
    const auto process = [](auto &&a, auto &&cm) {
      uint32_t merge = cm.pattern[0] == cm.pattern[1] ? 0 : 1;
      a.map(cm.pattern, std::make_pair(cm.mapping, merge));
    };

    std::vector<Peregrine::SmallGraph> patterns = {Peregrine::PatternGenerator::star(2)};
    patterns.front().set_labelling(Peregrine::Graph::DISCOVER_LABELS);
    auto psupps = Peregrine::match<Peregrine::Pattern, DiscoveryDomain<1>, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE>(data_graph_name, patterns, nthreads, process, view);
    for (const auto &[p, supp] : psupps)
    {
      if (supp >= threshold)
      {
        freq_patterns.push_back(p);
      }
    }

    auto tt2 = utils::get_timestamp();
    std::cout << "Step 0\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
  }

  std::vector<Peregrine::SmallGraph> patterns = Peregrine::PatternGenerator::extend(freq_patterns, Peregrine::PatternGenerator::EDGE_BASED);

  const auto process = [](auto &&a, auto &&cm) {
    a.map(cm.pattern, cm.mapping);
  };

  uint32_t step = 1;
  while (step < k && !patterns.empty())
  {
    freq_patterns.clear();
    supports.clear();

    auto tt1 = utils::get_timestamp();
    auto psupps = Peregrine::match<Peregrine::Pattern, Domain, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE>(data_graph_name, patterns, nthreads, process, view);

    for (const auto &[p, supp] : psupps)
    {
      if (supp >= threshold)
      {
        freq_patterns.push_back(p);
        supports.push_back(supp);
        REQUIRE CHECK_EQUAL(supp, truth[step-1].at(p));
      }
    }
    REQUIRE CHECK_EQUAL(truth[step-1].size(), freq_patterns.size());

    auto tt2 = utils::get_timestamp();
    std::cout << "Step " << step << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;

    patterns = Peregrine::PatternGenerator::extend(freq_patterns, Peregrine::PatternGenerator::EDGE_BASED);

    step += 1;
  }
  auto t2 = utils::get_timestamp();

  std::cout << "MATCHING INTEGRATION TEST FINISHED IN " << (t2-t1)/1e6 << "s" << std::endl;
}
