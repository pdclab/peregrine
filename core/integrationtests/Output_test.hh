#include "../Peregrine.hh"

SUITE(DiskOutputIntegration)
{
  // This seems like a unit test, but because it touches Peregrine::Context I've put it here
  TEST(OutputManager)
  {
    // check defaults
    const std::filesystem::path default_root = ".";

    Peregrine::SmallGraph p;
    p.add_edge(1, 2);
    p.add_edge(2, 3);
    Peregrine::Context::current_pattern = std::make_shared<Peregrine::AnalyzedPattern>(p);

    const std::filesystem::path default_path = default_root / p.to_string();
    CHECK_EQUAL(default_path, Peregrine::OutputManager<Peregrine::DISK>::get_path());

    // check set_root_directory
    const std::filesystem::path new_root = "./results";
    const std::filesystem::path new_path = new_root / p.to_string();

    Peregrine::OutputManager<Peregrine::DISK>::set_root_directory(new_root);
    CHECK_EQUAL(new_path, Peregrine::OutputManager<Peregrine::DISK>::get_path());
    CHECK(std::filesystem::is_directory(new_root));

    //
    // try setting up
    //
    unsigned id = 0;
    auto *bm = new Peregrine::OutputManager<Peregrine::DISK>;
    std::filesystem::path output_path = new_root / p.to_string() / std::to_string(id);
    bm->reset(id);
    CHECK(std::filesystem::exists(output_path));

    //
    // check outputs in BIN
    //
    Peregrine::DataGraph dg(p);
    Peregrine::Context::data_graph = &dg;

    bm->template output<Peregrine::BIN>({1, 2});
    bm->flush();
    CHECK(bm->persist());
    {
      // make sure we didn't write the entire page or anything
      CHECK_EQUAL(2*sizeof(uint32_t), std::filesystem::file_size(output_path));

      int fd = open(output_path.c_str(), O_RDONLY);
      CHECK(fd != -1);
      uint32_t read_buf[2];
      int bytes_read = read(fd, read_buf, 2*sizeof(uint32_t));
      CHECK_EQUAL(2*sizeof(uint32_t), bytes_read);
      CHECK_EQUAL(1, read_buf[0]);
      CHECK_EQUAL(2, read_buf[1]);
      close(fd);
    }

    //
    // check outputs in CSV
    //
    bm->reset(++id);
    output_path = new_root / p.to_string() / std::to_string(id);
    CHECK(std::filesystem::exists(output_path));
    bm->template output<Peregrine::CSV>({1, 2});
    bm->flush();
    CHECK(bm->persist());
    {

      // make sure we didn't write the entire page or anything
      // 4 bytes: u ',' v '\n'
      CHECK_EQUAL(4, std::filesystem::file_size(output_path));

      int fd = open(output_path.c_str(), O_RDONLY);
      CHECK(fd != -1);
      char read_buf[4];
      int bytes_read = read(fd, read_buf, 4);
      CHECK_EQUAL(4, bytes_read);
      CHECK_EQUAL('1', read_buf[0]);
      CHECK_EQUAL(',', read_buf[1]);
      CHECK_EQUAL('2', read_buf[2]);
      CHECK_EQUAL('\n', read_buf[3]);
      close(fd);
    }
    
    //
    // check for flush on reset
    //
    bm->reset(++id);
    output_path = new_root / p.to_string() / std::to_string(id);
    CHECK(std::filesystem::exists(output_path));
    bm->template output<Peregrine::CSV>({1, 2});
    bm->reset(++id);
    {
      CHECK(bm->persist());

      // make sure we didn't write the entire page or anything
      CHECK_EQUAL(4, std::filesystem::file_size(output_path));

      int fd = open(output_path.c_str(), O_RDONLY);
      CHECK(fd != -1);
      char read_buf[4];
      int bytes_read = read(fd, read_buf, 4);
      CHECK_EQUAL(4, bytes_read);
      CHECK_EQUAL('1', read_buf[0]);
      CHECK_EQUAL(',', read_buf[1]);
      CHECK_EQUAL('2', read_buf[2]);
      CHECK_EQUAL('\n', read_buf[3]);
      close(fd);
    }

    output_path = new_root / p.to_string() / std::to_string(id);
    CHECK(std::filesystem::exists(output_path));

    //
    // check that nothing is written when buffer is empty
    //
    bm->reset(++id);
    CHECK_EQUAL(std::filesystem::file_size(output_path), 0);

    // switch to current output_path
    output_path = new_root / p.to_string() / std::to_string(id);
    CHECK(std::filesystem::exists(output_path));

    //
    // check flush on full buffer for BIN
    //
    {
      const std::vector<uint32_t> vertices {1, 2, 3};
      if (PRG_OUTPUT_BLOCK_SIZE % 8 == 0)
      {
        size_t to_write = PRG_OUTPUT_BLOCK_SIZE/8 - 1;
        // write to_write many lines with 8 bytes each, then try to write one
        // with 12 bytes.
        // since PRG_OUTPUT_BLOCK_SIZE is a power of 2 greater than 8, it divides by 8 so
        // there will be exactly enough space for one more 8 byte entry in the buffer
        // and trying to write the 12 bytes will cause a flush. BIN doesn't
        // pack everything if it can't, so it will leave all 12 bytes in the buffer.
        for (size_t i = 0; i < to_write; ++i)
        {
          bm->template output<Peregrine::BIN>({1, 2});
        }

        CHECK(bm->persist());
        // shouldn't have flushed yet
        CHECK_EQUAL(0, std::filesystem::file_size(output_path));

        bm->template output<Peregrine::BIN>({1, 2, 3});
        CHECK(bm->persist());

        // should have caused a flush
        CHECK_EQUAL(PRG_OUTPUT_BLOCK_SIZE - 8, std::filesystem::file_size(output_path));

        // if we flush again the last entry should be written
        bm->flush();
        CHECK(bm->persist());
        CHECK_EQUAL(PRG_OUTPUT_BLOCK_SIZE+4, std::filesystem::file_size(output_path));
      }
      else
      {
        std::cerr << "WARNING: skipping flush BIN on full test because PRG_OUTPUT_BLOCK_SIZE % 8 != 0\n";
      }
    }

    //
    // check flush on full buffer for CSV
    //
    bm->reset(++id);
    output_path = new_root / p.to_string() / std::to_string(id);
    CHECK(std::filesystem::exists(output_path));
    {
      const std::vector<uint32_t> vertices {1, 2, 3};
      if (PRG_OUTPUT_BLOCK_SIZE % 4 == 0)
      {
        size_t to_write = PRG_OUTPUT_BLOCK_SIZE/4 - 1;
        // write to_write many lines with 4 bytes each, then try to write one
        // with 6 bytes.
        // since PRG_OUTPUT_BLOCK_SIZE is a power of 2, it divides by 4 so
        // there will be exactly enough space for one more line in the buffer
        // and writing the 6 bytes will cause a flush, leaving 2 bytes in the
        // buffer afterwards.
        for (size_t i = 0; i < to_write; ++i)
        {
          bm->template output<Peregrine::CSV>({1, 2});
        }

        CHECK(bm->persist());
        // shouldn't have flushed yet
        CHECK_EQUAL(0, std::filesystem::file_size(output_path));

        bm->template output<Peregrine::CSV>({1, 2, 3});
        
        CHECK(bm->persist());

        // should have caused a flush
        CHECK_EQUAL(PRG_OUTPUT_BLOCK_SIZE, std::filesystem::file_size(output_path));

        // if we flush again the remaining 2 bytes should be written
        bm->flush();
        CHECK(bm->persist());
        CHECK_EQUAL(PRG_OUTPUT_BLOCK_SIZE+2, std::filesystem::file_size(output_path));
      }
      else
      {
        std::cerr << "WARNING: skipping flush CSV on full test because PRG_OUTPUT_BLOCK_SIZE % 4 != 0\n";
      }
    }

    //
    // check for flush on destruction
    //
    bm->reset(++id);
    output_path = new_root / p.to_string() / std::to_string(id);
    bm->template output<Peregrine::CSV>({1, 2});
    CHECK(bm->persist());
    CHECK_EQUAL(std::filesystem::file_size(output_path), 0);
    {
      delete bm;
      int fd = open(output_path.c_str(), O_RDONLY);
      CHECK(fd != -1);
      CHECK_EQUAL(0, fsync(fd));
      char read_buf[4];
      int bytes_read = read(fd, read_buf, 4);
      CHECK_EQUAL(4, bytes_read);
      CHECK_EQUAL('1', read_buf[0]);
      CHECK_EQUAL(',', read_buf[1]);
      CHECK_EQUAL('2', read_buf[2]);
      CHECK_EQUAL('\n', read_buf[3]);
      close(fd);

      CHECK_EQUAL(4, std::filesystem::file_size(output_path));
    }

    // clean up
    std::filesystem::remove_all(new_root);
  }

  struct OutputIntegrationFixture
  {
    Peregrine::DataGraph dg;
    Peregrine::SmallGraph pattern;
    const size_t nthreads = std::thread::hardware_concurrency();
    std::filesystem::path dst;

    OutputIntegrationFixture()
      : dg("data/citeseer")
    {
      // triangle
      pattern.add_edge(1, 2);
      pattern.add_edge(1, 3);
      pattern.add_edge(2, 3);
      dst = std::filesystem::temp_directory_path();
    }

    ~OutputIntegrationFixture()
    {
    }
  };

  TEST_FIXTURE(OutputIntegrationFixture, SingleValueAggregator)
  {
    std::cout << "RUNNING SingleValueAggregator OUTPUT INTEGRATION TEST" << std::endl;

    const uint32_t truth = 1166;

    auto t1 = utils::get_timestamp();
    {
      {
        std::cout << "\tOUTPUT ALL TRIANGLES" << std::flush;
        auto tt1 = utils::get_timestamp();
        auto res = Peregrine::output<Peregrine::BIN>(dg, {pattern}, nthreads, dst);
        for (auto &[p, count, path] : res)
        {
          CHECK_EQUAL(truth, count);
          CHECK_EQUAL(dst / pattern.to_string(), path);
          size_t total_size = 0;
          for (unsigned i = 0; i < nthreads; ++i)
          {
            total_size += std::filesystem::file_size(path / std::to_string(i));
          }
          CHECK_EQUAL(truth*pattern.num_vertices()*sizeof(uint32_t), total_size);
        }
        auto tt2 = utils::get_timestamp();
        std::cout << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
      }

      {
        std::cout << "\tOUTPUT EVERY SECOND TRIANGLE" << std::flush;
        std::atomic<uint32_t> count(0);
        const auto process = [&count](auto &&a, auto &&cm) { if (count++ % 2 == 0) { a.template output<Peregrine::BIN>(cm.mapping); } a.map(cm.pattern, 1); };
        const auto view = [](uint64_t v) { return v; };
        auto tt1 = utils::get_timestamp();
        auto res = Peregrine::output<Peregrine::Pattern, uint64_t, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE>(dg, {pattern}, nthreads, process, view, dst);
        for (auto &[p, count, path] : res)
        {
          CHECK_EQUAL(truth, count);
          CHECK_EQUAL(dst / pattern.to_string(), path);
          size_t total_size = 0;
          for (unsigned i = 0; i < nthreads; ++i)
          {
            total_size += std::filesystem::file_size(path / std::to_string(i));
          }
          CHECK_EQUAL((truth/2)*pattern.num_vertices()*sizeof(uint32_t), total_size);
        }
        auto tt2 = utils::get_timestamp();
        std::cout << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
      }

      {
        std::cout << "\tOUTPUT USING match()" << std::flush;
        std::atomic<uint32_t> count(0);
        const auto process = [&count](auto &&a, auto &&cm) { if (count++ % 2 == 0) { a.template output<Peregrine::BIN>(cm.mapping); } a.map(cm.pattern, 1); };
        const auto view = [](uint64_t v) { return v; };
        auto tt1 = utils::get_timestamp();
        Peregrine::OutputManager<Peregrine::DISK>::set_root_directory(dst);
        auto res = Peregrine::match<Peregrine::Pattern, uint64_t, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE, Peregrine::DataGraph &, decltype(process), decltype(view), Peregrine::DISK>(dg, {pattern}, nthreads, std::move(process), view);
        for (auto &[p, count, path] : res)
        {
          CHECK_EQUAL(truth, count);
          CHECK_EQUAL(dst / pattern.to_string(), path);
          size_t total_size = 0;
          for (unsigned i = 0; i < nthreads; ++i)
          {
            total_size += std::filesystem::file_size(path / std::to_string(i));
          }
          CHECK_EQUAL((truth/2)*pattern.num_vertices()*sizeof(uint32_t), total_size);
        }
        auto tt2 = utils::get_timestamp();
        std::cout << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
      }
    }

    auto t2 = utils::get_timestamp();

    std::cout << "SingleValueAggregator OUTPUT INTEGRATION TEST FINISHED IN " << (t2-t1)/1e6 << "s" << std::endl;
  }

  TEST_FIXTURE(OutputIntegrationFixture, VectorAggregator)
  {
    std::cout << "RUNNING VectorAggregator OUTPUT INTEGRATION TEST" << std::endl;

    const uint32_t truth = 92;
    pattern.set_label(1, 1);
    pattern.set_label(2, 1);
    pattern.set_label(3, -1);

    auto t1 = utils::get_timestamp();
    {
      {
        std::cout << "\tOUTPUT ALL TRIANGLES" << std::flush;
        auto tt1 = utils::get_timestamp();
        auto res = Peregrine::output<Peregrine::BIN>(dg, {pattern}, nthreads, dst);
        size_t total_count = 0;
        size_t total_size = 0;
        for (auto &[p, count, path] : res)
        {
          total_count += count;
          CHECK_EQUAL(dst / pattern.to_string(), path);
        }

        for (unsigned i = 0; i < nthreads; ++i)
        {
          total_size += std::filesystem::file_size(dst / pattern.to_string() / std::to_string(i));
        }

        CHECK_EQUAL(truth*pattern.num_vertices()*sizeof(uint32_t), total_size);
        CHECK_EQUAL(truth, total_count);
        auto tt2 = utils::get_timestamp();
        std::cout << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
      }

      {
        std::cout << "\tOUTPUT EVERY SECOND TRIANGLE" << std::flush;
        std::atomic<uint32_t> count(0);
        const auto process = [&count](auto &&a, auto &&cm) { if (count++ % 2 == 0) { a.template output<Peregrine::BIN>(cm.mapping); } a.map(cm.pattern, 1); };
        const auto view = [](uint64_t v) { return v; };
        auto tt1 = utils::get_timestamp();
        auto res = Peregrine::output<Peregrine::Pattern, uint64_t, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE>(dg, {pattern}, nthreads, process, view, dst);

        size_t total_count = 0;
        size_t total_size = 0;
        for (auto &[p, count, path] : res)
        {
          total_count += count;
          CHECK_EQUAL(dst / pattern.to_string(), path);
        }

        for (unsigned i = 0; i < nthreads; ++i)
        {
          total_size += std::filesystem::file_size(dst / pattern.to_string() / std::to_string(i));
        }

        CHECK_EQUAL((truth/2)*pattern.num_vertices()*sizeof(uint32_t), total_size);
        CHECK_EQUAL(truth, total_count);
        auto tt2 = utils::get_timestamp();
        std::cout << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
      }

      {
        std::cout << "\tOUTPUT USING match()" << std::flush;
        std::atomic<uint32_t> count(0);
        const auto process = [&count](auto &&a, auto &&cm) { if (count++ % 2 == 0) { a.template output<Peregrine::BIN>(cm.mapping); } a.map(cm.pattern, 1); };
        const auto view = [](uint64_t v) { return v; };
        auto tt1 = utils::get_timestamp();
        Peregrine::OutputManager<Peregrine::DISK>::set_root_directory(dst);
        auto res = Peregrine::match<Peregrine::Pattern, uint64_t, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE, Peregrine::DataGraph &, decltype(process), decltype(view), Peregrine::DISK>(dg, {pattern}, nthreads, std::move(process), view);

        size_t total_count = 0;
        size_t total_size = 0;
        for (auto &[p, count, path] : res)
        {
          total_count += count;
          CHECK_EQUAL(dst / pattern.to_string(), path);
        }

        for (unsigned i = 0; i < nthreads; ++i)
        {
          total_size += std::filesystem::file_size(dst / pattern.to_string() / std::to_string(i));
        }

        CHECK_EQUAL((truth/2)*pattern.num_vertices()*sizeof(uint32_t), total_size);
        CHECK_EQUAL(truth, total_count);
        auto tt2 = utils::get_timestamp();
        std::cout << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
      }
    }

    auto t2 = utils::get_timestamp();

    std::cout << "VectorAggregator OUTPUT INTEGRATION TEST FINISHED IN " << (t2-t1)/1e6 << "s" << std::endl;
  }

  TEST_FIXTURE(OutputIntegrationFixture, GeneralPurposeAggregator)
  {
    std::cout << "RUNNING Aggregator OUTPUT INTEGRATION TEST" << std::endl;

    const uint32_t truth = 26878;
    // label discovery wedge
    pattern.remove_edge(2, 3);
    pattern.set_labelling(Peregrine::Graph::DISCOVER_LABELS);

    auto t1 = utils::get_timestamp();
    {
      {
        std::cout << "\tOUTPUT ALL TRIANGLES" << std::flush;
        auto tt1 = utils::get_timestamp();
        auto res = Peregrine::output<Peregrine::BIN>(dg, {pattern}, nthreads, dst);
        size_t total_count = 0;
        size_t total_size = 0;
        for (auto &[p, count, path] : res)
        {
          total_count += count;
          CHECK_EQUAL(dst / pattern.to_string(), path);
        }

        for (unsigned i = 0; i < nthreads; ++i)
        {
          total_size += std::filesystem::file_size(dst / pattern.to_string() / std::to_string(i));
        }

        CHECK_EQUAL(truth, total_count);
        CHECK_EQUAL(truth*pattern.num_vertices()*sizeof(uint32_t), total_size);
        auto tt2 = utils::get_timestamp();
        std::cout << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
      }

      {
        std::cout << "\tOUTPUT EVERY SECOND TRIANGLE" << std::flush;
        std::atomic<uint32_t> count(0);
        const auto process = [&count](auto &&a, auto &&cm) { if (count++ % 2 == 0) { a.template output<Peregrine::BIN>(cm.mapping); } a.map(cm.pattern, 1); };
        const auto view = [](uint64_t v) { return v; };
        auto tt1 = utils::get_timestamp();
        auto res = Peregrine::output<Peregrine::Pattern, uint64_t, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE>(dg, {pattern}, nthreads, process, view, dst);

        size_t total_count = 0;
        size_t total_size = 0;
        for (auto &[p, count, path] : res)
        {
          total_count += count;
          CHECK_EQUAL(dst / pattern.to_string(), path);
        }

        for (unsigned i = 0; i < nthreads; ++i)
        {
          total_size += std::filesystem::file_size(dst / pattern.to_string() / std::to_string(i));
        }

        CHECK_EQUAL(truth, total_count);
        CHECK_EQUAL((truth/2)*pattern.num_vertices()*sizeof(uint32_t), total_size);
        auto tt2 = utils::get_timestamp();
        std::cout << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
      }

      {
        std::cout << "\tOUTPUT USING match()" << std::flush;
        std::atomic<uint32_t> count(0);
        const auto process = [&count](auto &&a, auto &&cm) { if (count++ % 2 == 0) { a.template output<Peregrine::BIN>(cm.mapping); } a.map(cm.pattern, 1); };
        const auto view = [](uint64_t v) { return v; };
        auto tt1 = utils::get_timestamp();
        Peregrine::OutputManager<Peregrine::DISK>::set_root_directory(dst);
        auto res = Peregrine::match<Peregrine::Pattern, uint64_t, Peregrine::AT_THE_END, Peregrine::UNSTOPPABLE, Peregrine::DataGraph &, decltype(process), decltype(view), Peregrine::DISK>(dg, {pattern}, nthreads, std::move(process), view);

        size_t total_count = 0;
        size_t total_size = 0;
        for (auto &[p, count, path] : res)
        {
          total_count += count;
          CHECK_EQUAL(dst / pattern.to_string(), path);
        }

        for (unsigned i = 0; i < nthreads; ++i)
        {
          total_size += std::filesystem::file_size(dst / pattern.to_string() / std::to_string(i));
        }

        CHECK_EQUAL(truth, total_count);
        CHECK_EQUAL((truth/2)*pattern.num_vertices()*sizeof(uint32_t), total_size);
        auto tt2 = utils::get_timestamp();
        std::cout << "\t✓ " << (tt2-tt1)/1e6 << "s" << std::endl;
      }
    }

    auto t2 = utils::get_timestamp();

    std::cout << "Aggregator OUTPUT INTEGRATION TEST FINISHED IN " << (t2-t1)/1e6 << "s" << std::endl;
  }
}

