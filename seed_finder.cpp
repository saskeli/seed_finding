#include "include/seed_finder.hpp"

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
-mk <val>    Maximum mer length. ()"
            << max_k << R"().
-t <val>     Total number of threads to use. ()"
            << threads << R"(
-s           Disable smoothing of counted mers.
-f           Disable mer filtering.
--ouput_mers Ouput found siginificant mers to standard out.
-mem <val>   Memory limit for lookup tables (ish). ()"
            << gigs << " GB).\n\n"
            << std::endl;
  exit(0);
}

const constexpr uint16_t max_gap = MAX_GAP;

void filter_seeds(auto& seeds) {
  std::cerr << "Filtering.." << std::endl;
  uint64_t discarded_seeds = 0;
  size_t trg = 0;
  for (size_t i = 0; i < seeds.size(); ++i) {
    bool keep = true;
    for (size_t j = i + 1; j < seeds.size(); ++j) {
      if (seeds[i].g.aligns_to(seeds[j].g)) {
        keep = false;
        ++discarded_seeds;
        break;
      }
    }
    if (keep) {
      seeds[trg++] = seeds[i];
    }
  }
  seeds.resize(trg);
  std::cerr << discarded_seeds << " seeds discarded in filtering" << std::endl;
}

template <bool middle_gap_only, bool do_smoothing, bool do_filtering>
void run(const char* sig_path, const char* bg_path, double p, double log_fold,
         uint8_t max_k, double memory_limit, double p_ext, bool output_mers) {
  sf::seed_finder<true, max_gap, do_smoothing, do_filtering> sf(
      sig_path, bg_path, p, log_fold, max_k, memory_limit, p_ext);
  sf.find_seeds();
  auto seeds = sf.get_seeds();
  if constexpr (do_filtering) {
    filter_seeds(seeds);
  }
  if (output_mers) {
    for (auto seed : seeds) {
      std::cout << seed.g.to_string() << "\t" << seed.p << "\t("
                << seed.sig_count - 1 << ", " << seed.bg_count - 1 << ")"
                << std::endl;
    }
  }
  auto alignment = sf.align();
  for (size_t i = 0; i < alignment.size(); ++i) {
    if (i >= 4) {
      break;
    }
    std::string fname("align/al_");
    fname.append(std::to_string(i));
    fname.append(".fa");
    std::cerr << "Dumping " << alignment[i].size() << " sequences to " << fname << std::endl;
    std::ofstream outf(fname);
    size_t seq_c = 0;
    for (auto s : alignment[i]) {
      outf << "> seq " << seq_c << "\n";
      outf << s << "\n";
    }
  }
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
  bool enable_filtering = true;
  bool output_mers = false;
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
    } else if (arg == "-f") {
      enable_filtering = false;
    } else if (arg == "--output_mers") {
      output_mers = true;
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
#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif

  mem_limit *= 1000;
  mem_limit *= 1000;
  mem_limit *= 1000;

  if (enable_filtering) {
    if (enable_smoothing) {
      if (middle_gap_only) {
        run<true, true, true>(sig_path.c_str(), bg_path.c_str(), p, log_fold,
                              max_k, mem_limit, p_ext, output_mers);
      } else {
        run<false, true, true>(sig_path.c_str(), bg_path.c_str(), p, log_fold,
                               max_k, mem_limit, p_ext, output_mers);
      }
    } else {
      if (middle_gap_only) {
        run<true, false, true>(sig_path.c_str(), bg_path.c_str(), p, log_fold,
                               max_k, mem_limit, p_ext, output_mers);
      } else {
        run<false, false, true>(sig_path.c_str(), bg_path.c_str(), p, log_fold,
                                max_k, mem_limit, p_ext, output_mers);
      }
    }
  } else {
    if (enable_smoothing) {
      if (middle_gap_only) {
        run<true, true, false>(sig_path.c_str(), bg_path.c_str(), p, log_fold,
                               max_k, mem_limit, p_ext, output_mers);
      } else {
        run<false, true, false>(sig_path.c_str(), bg_path.c_str(), p, log_fold,
                                max_k, mem_limit, p_ext, output_mers);
      }
    } else {
      if (middle_gap_only) {
        run<true, false, false>(sig_path.c_str(), bg_path.c_str(), p, log_fold,
                                max_k, mem_limit, p_ext, output_mers);
      } else {
        run<false, false, false>(sig_path.c_str(), bg_path.c_str(), p, log_fold,
                                 max_k, mem_limit, p_ext, output_mers);
      }
    }
  }
  return 0;
}
