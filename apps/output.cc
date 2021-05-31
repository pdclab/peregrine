#include <sys/types.h>
#include <sys/stat.h>
#include <filesystem>
#include <unistd.h>

#include "Peregrine.hh"


int main(int argc, char *argv[])
{
  if (argc < 5)
  {
    std::cerr << "USAGE: " << argv[0] << " <data graph> <pattern | #-motifs | #-clique> <output directory> <bin | csv> [# threads]" << std::endl;
    return -1;
  }

  const std::string data_graph_name(argv[1]);
  const std::string pattern_name(argv[2]);
  const std::string output_dir(argv[3]);
  const std::string fmt_string(argv[4]);
  size_t nthreads = argc < 6 ? 1 : std::stoi(argv[5]);

  std::vector<Peregrine::SmallGraph> patterns;
  if (auto end = pattern_name.rfind("motifs"); end != std::string::npos)
  {
    auto k = std::stoul(pattern_name.substr(0, end-1));
    patterns = Peregrine::PatternGenerator::all(k,
        Peregrine::PatternGenerator::VERTEX_BASED,
        Peregrine::PatternGenerator::INCLUDE_ANTI_EDGES);
  }
  else if (auto end = pattern_name.rfind("clique"); end != std::string::npos)
  {
    auto k = std::stoul(pattern_name.substr(0, end-1));
    patterns.emplace_back(Peregrine::PatternGenerator::clique(k));
  }
  else
  {
    patterns.emplace_back(pattern_name);
  }


  Peregrine::DataGraph G(data_graph_name);
  std::vector<std::tuple<Peregrine::SmallGraph, uint64_t, std::filesystem::path>> result;

  if (fmt_string == "csv")
  {
    result = Peregrine::output<Peregrine::CSV>(G, patterns, nthreads, output_dir);
  }
  else if (fmt_string == "bin")
  {
    result = Peregrine::output<Peregrine::BIN>(G, patterns, nthreads, output_dir);
  }
  else
  {
    std::cerr << "Unrecognized output format: '" << fmt_string << "'\n"
      << "Currently supported output formats are 'csv' and 'bin'\n";
    return -1;
  }

  for (const auto &[p, c, d] : result)
  {
    std::cout << p << ": " << c << " matches stored in " << d << std::endl;
  }

  return 0;
}
