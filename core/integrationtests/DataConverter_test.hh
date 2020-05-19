#include <filesystem>

#include "../DataConverter.hh"
#include "../Peregrine.hh"


TEST(DataConverterIntegration)
{
  // # vertices in data graph
  const uint32_t dg_size = 6;
  const uint32_t nthreads = std::thread::hardware_concurrency();

  // directory for writing edge list, label file, and converted graph
  std::filesystem::path dir = std::filesystem::temp_directory_path() /
    ("DataConverter_test_" + std::to_string(utils::get_timestamp()));

  std::cout << "RUNNING DATA CONVERTER INTEGRATION TEST" << std::endl;

  auto t1 = utils::get_timestamp();
  {
    std::filesystem::create_directory(dir);

    // generate edge list and labels
    std::vector<std::pair<uint32_t, uint32_t>> edges;
    for (uint32_t i = 1; i < dg_size; ++i)
    {
      for (uint32_t j = i+1; j <= dg_size; ++j)
      {
        edges.emplace_back(i, j);
      }
    }
    // we want a dg_size-clique
    CHECK(Peregrine::AnalyzedPattern(edges).is_clique());

    std::vector<uint32_t> labels(dg_size+1, 0);

    // write the edges and labels to disk
    {
      std::ofstream edge_file(dir / "edges.txt");
      // add comments to the file
      edge_file << "# hello\n";
      edge_file << "# useless comment\n";
      edge_file << "# n 3 m 10\n";
      for (const auto &[u, v] : edges) edge_file << u << " " << v << "\n";
    }

    {
      std::ofstream label_file(dir / "labels.txt");
      // add comments to the file
      label_file << "# hello\n";
      label_file << "# useless comment\n";
      for (uint32_t i = 1; i <= dg_size; ++i) label_file << i << " " << labels[i] << "\n";
    }


    // convert the data
    Peregrine::DataConverter::convert_data(dir / "edges.txt", dir / "labels.txt", dir);

    // make sure number of edges/vertices is correct
    Peregrine::DataGraph g(dir);
    CHECK_EQUAL(dg_size, g.get_vertex_count());
    CHECK_EQUAL(edges.size(), g.get_edge_count());

    // test some matching to make sure edges and labels were written ok
    for (uint32_t k = 3; k <= dg_size; ++k)
    {
      auto &&kclique = Peregrine::PatternGenerator::clique(k);
      for (uint32_t i = 1; i <= k; ++i)
      {
        kclique.set_label(i, 0);
      }

      std::cout << k << "-clique" << std::flush;
      auto &&res = Peregrine::count(dir, {kclique}, nthreads);
      for (const auto &[_, count] : res)
      {
        CHECK_EQUAL(Peregrine::binom(dg_size, k), count);
      }
      std::cout << "\t✓" << std::endl;

      // set a label that didn't exist before
      kclique.set_label(1, 1);
      res = Peregrine::count(dir, {kclique}, nthreads);
      for (const auto &[_, count] : res)
      {
        CHECK_EQUAL(0, count);
      }
    }

    std::filesystem::remove_all(dir);
  }
  auto t2 = utils::get_timestamp();
  std::cout << "DATA CONVERTER INTEGRATION TEST FINISHED IN " << (t2-t1)/1e6 << "s" << std::endl;
}
