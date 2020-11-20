#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "Peregrine.hh"

#include "Domain.hh"

int main(int argc, char *argv[])
{
  if (argc < 4)
  {
    std::cerr << "USAGE: " << argv[0] << " <data graph> <max size> <support threshold> [vertex-induced] [# threads]" << std::endl;
    return -1;
  }

  const std::string data_graph_name(argv[1]);
  uint32_t k = std::stoi(argv[2]);

  uint64_t threshold = std::stoi(argv[3]);

  size_t nthreads = std::thread::hardware_concurrency();
  bool extension_strategy = Peregrine::PatternGenerator::EDGE_BASED;

  uint32_t step = 1;

  // decide whether user provided # threads or extension strategy
  if (argc == 5)
  {
    std::string arg(argv[4]);
    if (arg.starts_with("v")) // asking for vertex-induced
    {
      extension_strategy = Peregrine::PatternGenerator::VERTEX_BASED;
      step = 2;
    }
    else if (!arg.starts_with("e")) // not asking for edge-induced
    {
      nthreads = std::stoi(arg);
    }
  }
  else if (argc == 6)
  {
    for (std::string arg : {argv[4], argv[5]})
    {
      if (arg.starts_with("v")) // asking for vertex-induced
      {
        extension_strategy = Peregrine::PatternGenerator::VERTEX_BASED;
        step = 2;
      }
      else if (!arg.starts_with("e")) // not asking for edge-induced
      {
        nthreads = std::stoi(arg);
      }
    }
  }


  const auto view = [](auto &&v) { return v.get_support(); };

  std::vector<uint64_t> supports;
  std::vector<Peregrine::SmallGraph> freq_patterns;

  std::cout << k << "-FSM with threshold " << threshold << std::endl;

  Peregrine::DataGraph dg(data_graph_name);

  // initial discovery
  auto t1 = utils::get_timestamp();
  {
    const auto process = [](auto &&a, auto &&cm) {
      uint32_t merge = cm.pattern[0] == cm.pattern[1] ? 0 : 1;
      a.map(cm.pattern, std::make_pair(cm.mapping, merge));
    };

    std::vector<Peregrine::SmallGraph> patterns = {Peregrine::PatternGenerator::star(2)};
    patterns.front().set_labelling(Peregrine::Graph::DISCOVER_LABELS);
    auto psupps = Peregrine::match<Peregrine::Pattern, DiscoveryDomain<1>, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE>(dg, patterns, nthreads, process, view);
    for (const auto &[p, supp] : psupps)
    {
      if (supp >= threshold)
      {
        freq_patterns.push_back(p);
        supports.push_back(supp);
      }
    }
  }

  std::vector<Peregrine::SmallGraph> patterns = Peregrine::PatternGenerator::extend(freq_patterns, extension_strategy);

  const auto process = [](auto &&a, auto &&cm) {
    a.map(cm.pattern, cm.mapping);
  };

  while (step < k && !patterns.empty())
  {
    freq_patterns.clear();
    supports.clear();
    auto psupps = Peregrine::match<Peregrine::Pattern, Domain, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE>(dg, patterns, nthreads, process, view);

    for (const auto &[p, supp] : psupps)
    {
      if (supp >= threshold)
      {
        freq_patterns.push_back(p);
        supports.push_back(supp);
      }
    }

    patterns = Peregrine::PatternGenerator::extend(freq_patterns, extension_strategy);
    step += 1;
  }
  auto t2 = utils::get_timestamp();

  std::cout << freq_patterns.size() << " frequent patterns: " << std::endl;
  for (uint32_t i = 0; i < freq_patterns.size(); ++i)
  {
    std::cout << freq_patterns[i].to_string() << ": " << supports[i] << std::endl;
  }

  std::cout << "finished in " << (t2-t1)/1e6 << "s" << std::endl;
  return 0;
}
