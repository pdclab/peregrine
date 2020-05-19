#include "../Peregrine.hh"

TEST(EarlyTerminationIntegration)
{
  const std::string data_graph_name("data/citeseer");
  size_t nthreads = std::thread::hardware_concurrency();

  std::cout << "RUNNING EARLY TERMINATION INTEGRATION TEST" << std::endl;

  auto t1 = utils::get_timestamp();

  {
    // maximum clique in citeseer
    uint32_t k = 6;
    Peregrine::SmallGraph kclique = Peregrine::PatternGenerator::clique(k);

    auto process = [](auto &&a, auto &&cm) {
      a.map(cm.pattern, true);
      a.stop();
    };

    std::vector<std::pair<Peregrine::SmallGraph, bool>> results = Peregrine::match<Peregrine::Pattern, bool, Peregrine::AT_THE_END, Peregrine::STOPPABLE>(data_graph_name, {kclique}, nthreads, process);

    CHECK(results.front().second);
    CHECK(Peregrine::Context::task_ctr < Peregrine::DataGraph(data_graph_name).get_vertex_count());
    std::cout << k << "-clique exists" << "\t✓" << std::endl;
  }

  {
    // larger than maximum clique in citeseer
    uint32_t k = 7;
    Peregrine::SmallGraph kclique = Peregrine::PatternGenerator::clique(k);

    auto process = [](auto &&a, auto &&cm) {
      a.map(cm.pattern, true);
      a.stop();
    };

    std::vector<std::pair<Peregrine::SmallGraph, bool>> results = Peregrine::match<Peregrine::Pattern, bool, Peregrine::ON_THE_FLY, Peregrine::STOPPABLE>(data_graph_name, {kclique}, nthreads, process);

    CHECK(!results.front().second);
    CHECK(Peregrine::Context::task_ctr >= Peregrine::DataGraph(data_graph_name).get_vertex_count());
    std::cout << k << "-clique does not exist" << "\t✓" << std::endl;
  }

  auto t2 = utils::get_timestamp();
  std::cout << "EARLY TERMINATION INTEGRATION TEST FINISHED IN " << (t2-t1)/1e6 << "s" << std::endl;
}

