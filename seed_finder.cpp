#include <iostream>
#include <string>

#include "include/seed_finder.hpp"

void help(const char* call) {
  std::cout << R"(
Attempt seed extraction from selex data.

Usage: )" << call
            << R"( [OPTION]... [-b bg_fasta] <sig_fasta>

-b bg_fasta  Background fasta file. Required for now.
sig_fasta    Signal fasta file.
-a           Gap any location, not just middle.
-p <val>     p value to use. default 0.01.
-h           Print this and terminate. Overrides all other options.

)" << std::endl;
  exit(0);
}

const constexpr uint16_t max_gap = MAX_GAP;

int main(int argc, char const* argv[]) {
  std::string bg_path = "";
  std::string sig_path = "";
  bool middle_gap_only = true;
  double p = 0.01;
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "-h") {
      help(argv[0]);
    } else if (arg == "-a") {
      middle_gap_only = false;
    } else if (arg == "-b") {
      bg_path = argv[++i];
    } else if (arg == "-p") {
      p = std::stod(argv[++i]);
    } else {
      sig_path = arg;
    }
  }
  if (bg_path.size() == 0 || sig_path.size() == 0) {
    std::cerr << "Input fasta files are required." << std::endl;
    help(argv[0]);
  }
  
  if (middle_gap_only) {
    sf::seed_finder<true, max_gap> sf(bg_path.c_str(), sig_path.c_str(), p);
    sf.find_seeds();
    for (auto g : sf.get_seeds()) {
      std::cout << g.to_string() << std::endl;
    }
  } else {
    sf::seed_finder<false, max_gap> sf(bg_path.c_str(), sig_path.c_str(), p);
    sf.find_seeds();
    for (auto g : sf.get_seeds()) {
      std::cout << g.to_string() << std::endl;
    }
  }

  return 0;
}
