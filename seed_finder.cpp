#include "include/seed_finder.hpp"
#include "include/seed_clusterer.hpp"

#ifndef MAX_GAP
#define MAX_GAP 5
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
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
          double log_fold, double p_ext, uint64_t threads) {
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
-mk <val>    Maximum mer length in [6, 24] range. ()"
            << max_k << R"().
-t <val>     Total number of threads to use. ()"
            << threads << R"()
-s           Disable smoothing of counted mers.
-mem <val>   Memory limit for lookup tables (ish). ()"
            << gigs << " GB).\n\n"
            << std::endl;
  exit(0);
}

const constexpr uint16_t max_gap = MAX_GAP;

void filter_seeds(auto& seeds, auto callback) {
  std::cerr << "Filtering.." << std::endl;
  uint64_t discarded_seeds = 0;
  auto it = seeds.begin();
  for (; it != seeds.end(); ++it) {
    auto iit = it;
    ++iit;
    bool is_partial = false;
    for (; iit != seeds.end(); ++iit) {
      auto a = (*it).g;
      auto b = (*iit).g;
      if (a.aligns_to(b)) {
        is_partial = true;
        ++discarded_seeds;
        break;
      }
    }
    if (not is_partial) {
      auto seed = *it;
      callback(seed);
    }
  }
  std::cerr << discarded_seeds << " seeds discarded in filtering" << std::endl;
}

int main(int argc, char const* argv[]) {
  std::string bg_path = "";
  std::string sig_path = "";
  bool middle_gap_only = true;
  double p = 0.01;
  double p_ext = 0.01;
  double log_fold = 0.5;
#ifdef _OPENMP
  uint64_t threads = omp_get_max_threads();
#else
  uint64_t threads = 1;
#endif
  uint8_t max_k = 20;
  double mem_limit = available_gigs();
  bool print_help = false;
  bool enable_smoothing = true;
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
    } else if (arg == "-t") {
      threads = std::stoull(argv[++i]);
    } else if (arg == "-s") {
      enable_smoothing = false;
    } else if (arg.starts_with("-")) {
      std::cerr << "Invalid parameter \"" << arg << "\"." << std::endl;
      print_help = true;
    } else {
      sig_path = arg;
    }
  }
  if (print_help) {
    help(argv[0], mem_limit, max_k, p, log_fold, p_ext, threads);
    exit(0);
  }
  if (bg_path.size() == 0 || sig_path.size() == 0) {
    std::cerr << "Input fasta files are required." << std::endl;
    help(argv[0], mem_limit, max_k, p, log_fold, p_ext, threads);
    exit(1);
  }
  if (max_k < 6 || max_k > 24) {
    std::cerr << "invalide value for maximum k: " << max_k
              << ", should be in [6, 24] range." << std::endl;
    help(argv[0], mem_limit, max_k, p, log_fold, p_ext, threads);
    exit(1);
  }
#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif

  mem_limit *= 1000;
  mem_limit *= 1000;
  mem_limit *= 1000;
  if (enable_smoothing) {
    if (middle_gap_only) {
      sf::seed_finder<true, max_gap, true, false> sf(
          sig_path, bg_path, p, log_fold, max_k, mem_limit, p_ext);
      sf.find_seeds();
      sf::seed_clusterer<true, max_gap, decltype(sf.get_seeds())> sc(sf.get_seeds(), sig_path, bg_path);
      for (size_t i = 0; i < 4; ++i) {
        if (sc.size() < 10) {
          break;
        }
        sc.output_cluster();
      }
    } else {
      sf::seed_finder<false, max_gap, true, false> sf(
          sig_path, bg_path, p, log_fold, max_k, mem_limit, p_ext);
      sf.find_seeds();
      sf::seed_clusterer<false, max_gap, decltype(sf.get_seeds())> sc(sf.get_seeds(), sig_path, bg_path);
      for (size_t i = 0; i < 4; ++i) {
        if (sc.size() < 10) {
          break;
        }
        sc.output_cluster();
      }
    }
  } else {
    if (middle_gap_only) {
      sf::seed_finder<true, max_gap, false, false> sf(
          sig_path, bg_path, p, log_fold, max_k, mem_limit, p_ext);
      sf.find_seeds();
      sf::seed_clusterer<true, max_gap, decltype(sf.get_seeds())> sc(sf.get_seeds(), sig_path, bg_path);
      for (size_t i = 0; i < 4; ++i) {
        if (sc.size() < 10) {
          break;
        }
        sc.output_cluster();
      }
    } else {
      sf::seed_finder<false, max_gap, false, false> sf(
          sig_path, bg_path, p, log_fold, max_k, mem_limit, p_ext);
      sf.find_seeds();
      sf::seed_clusterer<false, max_gap, decltype(sf.get_seeds())> sc(sf.get_seeds(), sig_path, bg_path);
      for (size_t i = 0; i < 4; ++i) {
        if (sc.size() < 10) {
          break;
        }
        sc.output_cluster();
      }
    }
  }

  return 0;
}
