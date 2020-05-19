#include "../Peregrine.hh"

TEST(OnTheFlyIntegration)
{
  const std::string data_graph_name("data/citeseer");
  size_t nthreads = std::thread::hardware_concurrency();

  std::cout << "RUNNING ON-THE-FLY AGGREGATOR INTEGRATION TEST" << std::endl;

  auto t1 = utils::get_timestamp();

  {
    // maximum clique in citeseer
    uint32_t k = 4;
    Peregrine::SmallGraph kclique = Peregrine::PatternGenerator::star(k);

    std::atomic_flag f;
    auto process = [&f](auto &&a, auto &&cm) {
      using namespace std::chrono_literals;
      if (!f.test_and_set()) // first thread only
      {
        uint64_t first_count = a.read_value(cm.pattern);
        std::this_thread::sleep_for(600ms); // wait for aggregator to run
        uint64_t second_count = a.read_value(cm.pattern);
        CHECK(first_count < second_count);
      }

      a.map(cm.pattern, 1);
    };

    auto result = Peregrine::match<Peregrine::Pattern, uint64_t, Peregrine::ON_THE_FLY, Peregrine::UNSTOPPABLE>(data_graph_name, {kclique}, nthreads, process);

    auto expected = Peregrine::count(data_graph_name, {kclique}, nthreads);

    CHECK_EQUAL(expected.front().second, result.front().second);
  }

  auto t2 = utils::get_timestamp();
  std::cout << "ON-THE-FLY AGGREGATOR INTEGRATION TEST FINISHED IN " << (t2-t1)/1e6 << "s" << std::endl;
}


