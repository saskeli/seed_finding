#include "include/seed_finder.hpp"

#include <unistd.h>

#include <iostream>
#include <string>

uint64_t available_gigs() {
  uint64_t mem = sysconf(_SC_PHYS_PAGES);
  mem *= sysconf(_SC_PAGE_SIZE);
  mem /= 2;
  mem /= 1000;
  mem /= 1000;
  uint64_t gigs = mem / 1000;
  if (mem - gigs * 1000 < (gigs + 1) * 1000 - mem) {
    ++gigs;
  }
  return gigs;
}

void help(const char* call, uint64_t gigs, uint64_t max_k, double p,
          double log_fold, double p_ext) {
  std::cerr << R"(
Attempt seed extraction from read data.

With max gap size = )"
            << MAX_GAP << R"(

Usage: )" << call
            << R"( [OPTION]... [-b bg_fasta] <sig_fasta>

-b bg_fasta  Background fasta file. Required for now.
sig_fasta    Signal fasta file.
-a           Gap any location, not just middle.
-p <val>     p value to use for signal to background count comparison. ()"
            << p << R"().
-pext <val>  p value to use for extension when background counts are 0. ()"
            << p_ext << R"().
-lf <val>    Discard all mers with log fold change smaller than this ()"
            << log_fold << R"().
-h           Print this and terminate. Overrides all other options.
-mk <val>    Maximum mer length. ()"
            << max_k << R"().
-mem <val>   Memory limit for lookup tables (ish). ()"
            << gigs << " GB).\n\n"
            << std::endl;
  exit(0);
}

const constexpr uint16_t max_gap = MAX_GAP;

int main(int argc, char const* argv[]) {
  std::string bg_path = "";
  std::string sig_path = "";
  bool middle_gap_only = true;
  double p = 0.01;
  double p_ext = 0.01;
  double log_fold = 0.5;
  uint8_t max_k = 20;
  double mem_limit = available_gigs();
  bool print_help = false;
  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "-h") {
      print_help = true;
    } else if (arg == "-a") {
      middle_gap_only = false;
    } else if (arg == "-b") {
      bg_path = argv[++i];
    } else if (arg == "-p") {
      p = std::stod(argv[++i]);
    } else if (arg == "-pext") {
      p_ext = std::stod(argv[++i]);
    } else if (arg == "-lf") {
      log_fold = std::stod(argv[++i]);
    } else if (arg == "-mk") {
      max_k = std::stoi(argv[++i]);
    } else if (arg == "-mem") {
      mem_limit = std::stod(argv[++i]);
    } else {
      sig_path = arg;
    }
  }
  if (print_help) {
    help(argv[0], mem_limit, max_k, p, log_fold, p_ext);
    exit(0);
  }
  if (bg_path.size() == 0 || sig_path.size() == 0) {
    std::cerr << "Input fasta files are required." << std::endl;
    help(argv[0], mem_limit, max_k, p, log_fold, p_ext);
    exit(1);
  }

  mem_limit *= 1000;
  mem_limit *= 1000;
  mem_limit *= 1000;

  if (middle_gap_only) {
    sf::seed_finder<true, max_gap> sf(sig_path.c_str(), bg_path.c_str(), p,
                                      log_fold, max_k, mem_limit, p_ext);
    sf.find_seeds();
    for (auto seed : sf.get_seeds()) {
      std::cout << seed.g.to_string() << "\t" << seed.p << "\t("
                << seed.sig_count - 1 << ", " << seed.bg_count - 1 << ")"
                << std::endl;
    }
  } else {
    sf::seed_finder<false, max_gap> sf(sig_path.c_str(), bg_path.c_str(), p,
                                       log_fold, max_k, mem_limit, p_ext);
    sf.find_seeds();
    for (auto seed : sf.get_seeds()) {
      std::cout << seed.g.to_string() << "\t" << seed.p << "\t("
                << seed.sig_count - 1 << ", " << seed.bg_count - 1 << ")"
                << std::endl;
    }
  }

  return 0;
}
