#include <iostream>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>


#include "DataConverter.hh"

bool is_directory(const std::string &path)
{
  struct stat statbuf;
  return (stat(path.c_str(), &statbuf) == 0) && S_ISDIR(statbuf.st_mode);
}

bool file_exists(const std::string &path)
{
  struct stat statbuf;
  return (stat(path.c_str(), &statbuf) == 0);
}

int main(int argc, char *argv[])
{
  if (argc < 3)
  {
    std::cerr
      << "USAGE: " << argv[0]
      << " <input edge file> [input label file] <output directory>"
      << std::endl;
    return -1;
  }

  std::string edge_file(argv[1]);
  std::string label_file(argc == 4 ? argv[2] : "");
  std::string out_dir(argc == 4 ? argv[3] : argv[2]);

  if (!file_exists(edge_file))
  {
    std::cerr
      << "ERROR: " << edge_file << " could not be opened."
      << std::endl;
    return -1;
  }
  else if (argc == 4 && !file_exists(label_file))
  {
    std::cerr
      << "ERROR: " << label_file << " could not be opened."
      << std::endl;
    return -1;
  }
  else if (!is_directory(out_dir))
  {
    std::cerr
      << "ERROR: " << out_dir << " is not a valid directory."
      << std::endl;
    return -1;
  }

  Peregrine::DataConverter::convert_data(edge_file, label_file, out_dir);

  return 0;
}
