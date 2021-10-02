#include <iostream>
#include <filesystem>
#include <string>

#include "DataConverter.hh"

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

  if (!std::filesystem::exists(edge_file))
  {
    std::cerr
      << "ERROR: " << edge_file << " could not be opened."
      << std::endl;
    return -1;
  }
  else if (argc == 4 && !std::filesystem::exists(label_file))
  {
    std::cerr
      << "ERROR: " << label_file << " could not be opened."
      << std::endl;
    return -1;
  }
  else if (!std::filesystem::is_directory(out_dir))
  {
    if (std::filesystem::exists(out_dir))
    {
      std::cerr
        << "ERROR: " << out_dir << " already exists and is not a directory."
        << std::endl;
      return -1;
    }
    else if (!std::filesystem::create_directory(out_dir))
    {
      std::cerr
        << "ERROR: could not create directory " << out_dir << "."
        << std::endl;
      return -1;
    }
  }

  Peregrine::DataConverter::convert_data(edge_file, label_file, out_dir);

  return 0;
}
