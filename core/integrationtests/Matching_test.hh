#include "../Peregrine.hh"
#include "../../apps/Domain.hh"

auto from_string(const std::string &str)
{
  Peregrine::SmallGraph p;
  {
    unsigned c = 0;
    while (c < str.size())
    {
      if (str[c] == '[')
      {
        unsigned edge_start = c;
        while (str[c] != ']' && c < str.size()) c += 1;
        // edges are in format [u,lu-v,lv] where each is a one digit integer
        assert(c-edge_start == 8);

        uint32_t u = str[edge_start+1] - '0';
        uint32_t lu = str[edge_start+3] - '0';
        uint32_t v = str[edge_start+5] - '0';
        uint32_t lv = str[edge_start+7] - '0';

        p.add_edge(u, v);
        p.set_label(u, lu);
        p.set_label(v, lv);
      }
      c += 1;
    }
  }

  return p;
}


TEST(MatchingIntegration)
{
  const std::string data_graph_name("data/citeseer");
  uint32_t k = 4;
  uint64_t threshold = 1;
  size_t nthreads = std::thread::hardware_concurrency();

  std::cout << "RUNNING MATCHING INTEGRATION TEST" << std::endl;

  std::unordered_map<Peregrine::SmallGraph, uint64_t> truth;
  {
    std::ifstream truth_file("core/integrationtests/cs-supports.txt");
    CHECK(truth_file.is_open());

    // Using the graph string representation may be expensive, but the test is
    // fast _enough_ and this makes the truth file human readable.
    std::string id;
    uint64_t supp;
    while (truth_file >> id >> supp)
    {
      truth[from_string(id)] = supp;
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
    auto psupps = Peregrine::match<Peregrine::Pattern, DiscoveryDomain, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE>(data_graph_name, patterns, nthreads, process, view);
    for (const auto &[p, supp] : psupps)
    {
      if (supp >= threshold)
      {
        freq_patterns.push_back(p);
        supports.push_back(supp);
        REQUIRE CHECK_EQUAL(supp, truth.at(p));
      }
    }
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
        REQUIRE CHECK_EQUAL(supp, truth.at(p));
      }
    }
    auto tt2 = utils::get_timestamp();
    std::cout << "Step " << step << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;

    patterns = Peregrine::PatternGenerator::extend(freq_patterns, Peregrine::PatternGenerator::EDGE_BASED);
    step += 1;
  }
  auto t2 = utils::get_timestamp();

  std::cout << "MATCHING INTEGRATION TEST FINISHED IN " << (t2-t1)/1e6 << "s" << std::endl;
}
