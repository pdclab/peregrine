#ifndef OUTPUT_MANAGER_HH
#define OUTPUT_MANAGER_HH

#include "Options.hh"

#include <vector>
#include <charconv>
#include <system_error>
#include <sys/stat.h>
#include <fcntl.h>


#define PRG_OUTPUT_BLOCK_SIZE 16777216

namespace Peregrine
{

// Make it empty in the general case, where it's useless
template <OutputOption Output>
class OutputManager
{
  void reset() = delete;
  void output(const std::vector<uint32_t> &) = delete;
  void flush() = delete;
};

template <>
class OutputManager<DISK>
{
  public:
    OutputManager()
    {
      bufsz = 0;
      buf = new char[PRG_OUTPUT_BLOCK_SIZE];
    }

    ~OutputManager()
    {
      if (bufsz > 0)
      {
        flush();
      }
      assert(bufsz == 0);
      delete[] buf;

      if (fd != -1) close(fd);
    }

    static const std::filesystem::path get_path()
    {
      return output_root / Context::current_pattern->query_graph.to_string();
    }

    static void set_root_directory(const std::filesystem::path &new_root)
    {
      output_root = new_root;
      std::filesystem::create_directory(output_root);
    }

    void reset(uint32_t id)
    {
      if (bufsz > 0)
      {
        flush();
      }

      if (fd != -1) // the first time, fd won't be a valid file descriptor
      {
        close(fd);
      }

      std::filesystem::create_directory(get_path());
      fd = open((get_path() / std::to_string(id)).c_str(), O_WRONLY | O_CREAT | O_APPEND | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
      assert(fd != -1);
    }

    template <OutputFormat fmt>
    void output(const std::vector<uint32_t> &vertices)
    {
      if constexpr (fmt == BIN)
      {
        if ((PRG_OUTPUT_BLOCK_SIZE-bufsz) <= vertices.size()*4)
        {
          flush();
        }

        std::transform(std::execution::unseq,
            vertices.cbegin(), vertices.cend(), reinterpret_cast<uint32_t *>(&buf[bufsz]),
            [](uint32_t u) { return Context::data_graph->original_id(u); });
        bufsz += 4*vertices.size();
      }
      else if constexpr (fmt == CSV)
      {
        std::for_each(vertices.cbegin(), std::prev(vertices.cend()),
            [this](uint32_t v)
            {
              auto [ptr, ec] = std::to_chars(&buf[bufsz], &buf[PRG_OUTPUT_BLOCK_SIZE], Context::data_graph->original_id(v));
              if (ec == std::errc::value_too_large)
              {
                flush();
                assert(bufsz == 0);
                // try again: it should only fail if the ID needs more than PRG_OUTPUT_BLOCK_SIZE digits... unlikely
                auto [_ptr, ec] = std::to_chars(&buf[bufsz], &buf[PRG_OUTPUT_BLOCK_SIZE], Context::data_graph->original_id(v));
                ptr = _ptr;
                assert(ec == std::errc());
              }
              else if (ptr == &buf[PRG_OUTPUT_BLOCK_SIZE]) // ptr is past-the-end, so we aren't out of space at BLOCK_SIZE-1
              {
                bufsz += std::distance(&buf[bufsz], ptr); // write what you can before flushing
                flush(); // make room for comma
                assert(bufsz == 0);
                ptr = &buf[bufsz]; // reset ptr: it pointed past-the-end
              }

              *ptr++ = ',';
              bufsz += std::distance(&buf[bufsz], ptr);
            });

        // do the last vertex
        {
          uint32_t v = vertices.back();

          auto [ptr, ec] = std::to_chars(&buf[bufsz], &buf[PRG_OUTPUT_BLOCK_SIZE], Context::data_graph->original_id(v));
          if (ec == std::errc::value_too_large)
          {
            flush();
            assert(bufsz == 0);
            // try again: it should only fail if the ID needs more than PRG_OUTPUT_BLOCK_SIZE digits... unlikely
            auto [_ptr, ec] = std::to_chars(&buf[bufsz], &buf[PRG_OUTPUT_BLOCK_SIZE], Context::data_graph->original_id(v));
            ptr = _ptr;
            assert(ec == std::errc());
          }
          else if (ptr == &buf[PRG_OUTPUT_BLOCK_SIZE]) // ptr is past-the-end, so we aren't out of space at BLOCK_SIZE-1
          {
            flush(); // make room for newline
            assert(bufsz == 0);
            ptr = &buf[bufsz]; // reset ptr: it pointed past-the-end
          }

          *ptr++ = '\n';
          bufsz += std::distance(&buf[bufsz], ptr);
        }
      }
      else
      {
        // Never happens
        throw std::runtime_error("Unrecognized output format."); 
      }
    }

    bool persist()
    {
      int err = fsync(fd);
      if (err != 0)
      {
        std::cerr << "fsync failed: " << strerror(errno) << std::endl;
      }
      return err == 0;
    }

    void flush()
    {
      ssize_t bytes = 0;
      while (bytes < bufsz)
      {
        bytes += write(fd, reinterpret_cast<void *>(&buf[bytes]), bufsz - bytes);
        if (bytes == -1)
        {
          throw std::runtime_error("Match output failed: " + std::string(strerror(errno)));
        }
      }
      assert(bytes == bufsz);
      bufsz -= bytes;
      assert(bufsz == 0);
    }

  private:
    static inline std::filesystem::path output_root = ".";
    int fd = -1;
    char *buf;
    uint32_t bufsz;
};

} // namespace Peregrine

#endif
