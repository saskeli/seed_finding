#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif
#include <unistd.h>

#include <args.hxx>

#include "include/args.hpp"
#include "include/dot_writer.hpp"
#include "include/seed_clusterer.hpp"
#include "include/seed_finder.hpp"
#include "include/version.hpp"

#ifndef MAX_GAP
#define MAX_GAP 5
#endif

#define SF_STRINGIFY_(X) #X
#define SF_STRINGIFY(X) SF_STRINGIFY_(X)


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

// We use uint16_t since it is easier to output correctly than uint8_t (i.e. not
// as a character).
const constexpr uint16_t max_gap = MAX_GAP;

void filter_seeds(auto &seeds, auto callback) {
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

int main(int argc, char const *argv[]) {
  std::string bg_path = "";
  std::string sig_path = "";
  std::string prefix = "";
  std::string dot_output = "";
  bool middle_gap_only = true;
  bool should_output_all_matches = false;
  double p = 0.0001;
  double p_ext = 0.01;
  double log_fold = 2;
  double h1_weight = 1;
  size_t print_lim = 20;
  size_t lookup_k = 0;
  size_t max_aligns = 0;
#ifdef _OPENMP
  uint64_t threads = omp_get_max_threads();
#else
  uint64_t threads = 1;
#endif
  uint16_t max_k = 20;
  double mem_limit = available_gigs();
  bool enable_smoothing = true;

  {
    args::ArgumentParser parser(
        "Attempt seed extraction from read data.\n\nWith max gap size "
        "= " SF_STRINGIFY(MAX_GAP) ".");
    parser.SetArgumentSeparations(true, true, true, true);

    {
      args::HelpParams help_params{};
      help_params.addDefault = true;
      help_params.shortPrefix = "-";
      help_params.longPrefix = "--";
      help_params.shortSeparator = " ";
      help_params.longSeparator = " ";
      help_params.defaultString = " Default: ";
      parser.helpParams = help_params;
    }

    args::HelpFlag help_(parser, "help", "Display this help.", {'h', "help"});
    sf::args::version_flag version_(parser, "version",
                                    "Output the version number and exit.",
                                    {'V', "version"});
    args::CompletionFlag completion_(parser, {"complete"});
    sf::args::value_flag bg_path_(
        parser, "background_path", "Background FASTA file, required for now.",
        {'b', "background"}, bg_path, args::Options::Required);
    args::Positional<std::string> sig_path_(parser, "signal_path",
                                            "Signal FASTA file.", sig_path,
                                            args::Options::Required);
    sf::args::value_flag dot_path_(
        parser, "dot_path",
        "Compute huddinge graph and output to dot file path.", {"dot"},
        dot_output);
    args::Flag gap_any_(parser, "gap_any",
                        "Allow gaps at any location, not just in the middle.",
                        {'a', "gap-at-any-location"});
    sf::args::value_flag p_sig_(
        parser, "p_signal",
        "p value to use for signal to background comparison.",
        {'p', "p-signal"}, p);
    sf::args::value_flag p_ext_(
        parser, "p_extension",
        "p value to use for extension when background counts are zero.",
        {"p-ext"}, p_ext);
    sf::args::value_flag log_fold_(
        parser, "lf",
        "Discard all mers with log fold change smaller than this.", {"lf"},
        log_fold);
    sf::args::value_flag max_k_(
        parser, "mk", "Maximum mer length in [6, 24] range.", {"mk"}, max_k);
    sf::args::value_flag lookup_k_(
        parser, "lokup_k",
        "Limit for lookup table-based k-mer counting in [5, max_k] range.",
        {"lookup_k"}, lookup_k);
    sf::args::value_flag threads_(parser, "threads",
                                  "Number of threads to use.", {'t', "threads"},
                                  threads);
    sf::args::value_flag prefix_(parser, "output_prefix",
                                 "Prefix for alignment output, If not given, "
                                 "no alignments will be output.",
                                 {"pref"}, prefix);
    sf::args::value_flag print_lim_(
        parser, "max_s",
        "Maximum number of “best” seeds to output (0 -> all seeds).", {"max-s"},
        print_lim);
    sf::args::value_flag max_aligns_(
        parser, "max_a",
        "Maximum number of alignments to output. (0 -> max_s).", {"max-a"},
        max_aligns);
    args::Flag should_output_all_matches_(
        parser, "output_all_matches",
        "Output also matches that are substrings of other matches.",
        {"output-all-matches"});
    args::Flag disable_smoothing_(parser, "disable_smoothing",
                                  "Disable smoothing of counted mers.",
                                  {'s', "no-smoothing"});
    sf::args::value_flag mem_limit_(
        parser, "memory_limit",
        "Approximate memory limit for lookup tables in gigabytes.", {"mem"},
        mem_limit);
    sf::args::value_flag h1_weight_(
        parser, "h1_weight",
        "Relative impact of H1 neighbourhood enrichment on mer priority. (0 -> "
        "no impact, 1 -> 0.5 h1 neighbourhod 0.5 mer enrichment)",
        {"h1_weight"}, h1_weight);

    // Parse and check.
    try {
      parser.ParseCLI(argc, argv);
    } catch (args::Help const &) {
      std::cerr << parser;
      std::exit(0);
    } catch (args::Completion const &) {
      std::cerr << parser;
      std::exit(0);
    } catch (sf::args::output_version const &) {
      std::cout << "seed_finder " << sf::version << '\n';
      std::exit(0);
    } catch (std::runtime_error const &err) {
      std::cerr << err.what() << '\n';
      std::exit(1);
    }

    sig_path = args::get(sig_path_);
    if (gap_any_) middle_gap_only = true;
    if (should_output_all_matches_) should_output_all_matches = true;
    if (disable_smoothing_) enable_smoothing = false;

    if (print_lim == 0) {
      print_lim = ~print_lim;
    }
    if (max_aligns == 0) {
      max_aligns = print_lim;
    }
  }

  if (max_k < 6 || max_k > 24) {
    std::cerr << "invalid value for maximum k: " << max_k
              << ", should be in [6, 24] range." << std::endl;
    exit(1);
  }

#ifdef DEBUG
  {
    std::cerr << "bg_path:                   " << bg_path << '\n';
    std::cerr << "sig_path:                  " << sig_path << '\n';
    std::cerr << "prefix:                    " << prefix << '\n';
    std::cerr << "middle_gap_only:           " << middle_gap_only << '\n';
    std::cerr << "should_output_all_matches: " << should_output_all_matches
              << '\n';
    std::cerr << "p:                         " << p << '\n';
    std::cerr << "p_ext:                     " << p_ext << '\n';
    std::cerr << "log_fold:                  " << log_fold << '\n';
    std::cerr << "print_lim:                 " << print_lim << '\n';
    std::cerr << "max_aligns:                " << max_aligns << "\n";
    std::cerr << "threads:                   " << threads << '\n';
    std::cerr << "max_k:                     " << max_k << '\n';
    std::cerr << "lookup_k:                  " << lookup_k << '\n';
    std::cerr << "mem_limit:                 " << mem_limit << '\n';
    std::cerr << "enable_smoothing:          " << enable_smoothing << '\n';
  }
#endif

  mem_limit *= 1000;
  mem_limit *= 1000;
  mem_limit *= 1000;

  if (lookup_k == 0) {
    lookup_k = 5;
    // figure out how big lookup tables will fit in memory.
    if (middle_gap_only) {
      while (sf::gapmer_count<true, max_gap>::lookup_bytes(lookup_k) <
             mem_limit) {
        ++lookup_k;
      }
    } else {
      while (sf::gapmer_count<false, max_gap>::lookup_bytes(lookup_k) <
             mem_limit) {
        ++lookup_k;
      }
    }
    if (lookup_k > 10) {
      lookup_k = 10;
    }
  }
  if (lookup_k < 5 || lookup_k > max_k) {
    std::cerr << "Invalid lookup_k value (" << lookup_k << " not in [" << 5
              << ", " << max_k << "]).";
    exit(1);
  }

#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif

  if (prefix.length() > 0) {
    std::filesystem::create_directory(prefix);
  }
  if (enable_smoothing) {
    if (middle_gap_only) {
      sf::seed_finder<true, max_gap, true, false> finder(
          sig_path, bg_path, p, log_fold, max_k, mem_limit, p_ext, lookup_k);
      finder.find_seeds();
      if (dot_output.size() > 0) {
        sf::Dot_Writer::write_dot<decltype(finder.get_seeds()),
                                  decltype(finder)::G>(
            dot_output, finder.get_seeds(), max_k);
      }
      sf::seed_clusterer<true, max_gap, decltype(finder.get_seeds())> sc(
          finder.get_seeds(), sig_path, bg_path, p_ext, h1_weight, finder.x());
      std::cout << "Seed\tcounts\tp\tpriority" << std::endl;
      for (size_t i = 0; i < print_lim; ++i) {
        if (not sc.has_next()) {
          break;
        }
        if (max_aligns == i) {
          prefix = "";
        }
        sc.output_cluster(prefix, should_output_all_matches);
      }
    } else {
      sf::seed_finder<false, max_gap, true, false> finder(
          sig_path, bg_path, p, log_fold, max_k, mem_limit, p_ext, lookup_k);
      finder.find_seeds();
      if (dot_output.size() > 0) {
        sf::Dot_Writer::write_dot<decltype(finder.get_seeds()),
                                  decltype(finder)::G>(
            dot_output, finder.get_seeds(), max_k);
      }
      sf::seed_clusterer<false, max_gap, decltype(finder.get_seeds())> sc(
          finder.get_seeds(), sig_path, bg_path, p_ext, h1_weight, finder.x());
      std::cout << "Seed\tcounts\tp\tpriority" << std::endl;
      for (size_t i = 0; i < print_lim; ++i) {
        if (not sc.has_next()) {
          break;
        }
        if (max_aligns == i) {
          prefix = "";
        }
        sc.output_cluster(prefix, should_output_all_matches);
      }
    }
  } else {
    if (middle_gap_only) {
      sf::seed_finder<true, max_gap, false, false> finder(
          sig_path, bg_path, p, log_fold, max_k, mem_limit, p_ext, lookup_k);
      finder.find_seeds();
      if (dot_output.size() > 0) {
        sf::Dot_Writer::write_dot<decltype(finder.get_seeds()),
                                  decltype(finder)::G>(
            dot_output, finder.get_seeds(), max_k);
      }
      sf::seed_clusterer<true, max_gap, decltype(finder.get_seeds())> sc(
          finder.get_seeds(), sig_path, bg_path, p_ext, h1_weight, finder.x());
      std::cout << "Seed\tcounts\tp\tpriority" << std::endl;
      for (size_t i = 0; i < print_lim; ++i) {
        if (not sc.has_next()) {
          break;
        }
        if (max_aligns == i) {
          prefix = "";
        }
        sc.output_cluster(prefix, should_output_all_matches);
      }
    } else {
      sf::seed_finder<false, max_gap, false, false> finder(
          sig_path, bg_path, p, log_fold, max_k, mem_limit, p_ext, lookup_k);
      finder.find_seeds();
      if (dot_output.size() > 0) {
        sf::Dot_Writer::write_dot<decltype(finder.get_seeds()),
                                  decltype(finder)::G>(
            dot_output, finder.get_seeds(), max_k);
      }
      sf::seed_clusterer<false, max_gap, decltype(finder.get_seeds())> sc(
          finder.get_seeds(), sig_path, bg_path, p_ext, h1_weight, finder.x());
      std::cout << "Seed\tcounts\tp\tpriority" << std::endl;
      for (size_t i = 0; i < print_lim; ++i) {
        if (not sc.has_next()) {
          break;
        }
        if (max_aligns == i) {
          prefix = "";
        }
        sc.output_cluster(prefix, should_output_all_matches);
      }
    }
  }

  return 0;
}
