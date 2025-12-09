#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "pack_characters.hpp"
#include "reader_adapter.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif
#include <unistd.h>

#include <args.hxx>

#include "include/args.hpp"
#include "include/configuration.hpp"
#include "include/dot_writer.hpp"
#include "include/gapmer_count.hpp"
#include "include/reader_adapter.hpp"
#include "include/seed_finder.hpp"
#include "seed_clusterer.hpp"
#include "include/util.hpp"
#include "include/version.hpp"

#ifndef MAX_GAP
#define MAX_GAP 5
#endif

#define SF_STRINGIFY_(X) #X
#define SF_STRINGIFY(X) SF_STRINGIFY_(X)


namespace {

inline uint64_t available_gigs() {
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


void print_invocation(int argc, char const* argv[]) {
  std::cerr << "Invocation:";
  for (int i{}; i < argc; ++i) std::cerr << ' ' << argv[i];
  std::cerr << '\n';
}


struct configuration {
  std::string bg_path{};
  std::string sig_path{};
  std::string prefix{};
  std::string dot_output{};
  bool middle_gap_only{true};
  bool prune{false};
  bool should_output_all_matches{false};
  double p{0.0001};
  double p_ext{0.01};
  double log_fold{2};
  double h1_weight{1};
  size_t print_lim{20};
  size_t lookup_k{};
  size_t max_aligns{};
  int threads{1};
  uint16_t max_k{20};
  double mem_limit{};
  bool enable_smoothing{true};

  configuration()
      :
#ifdef _OPENMP
        threads{omp_get_max_threads()},
#endif
        mem_limit(available_gigs()) {
  }
};


std::ostream& operator<<(std::ostream& os, configuration const& conf) {
  os << "bg_path:                   " << conf.bg_path << '\n';
  os << "sig_path:                  " << conf.sig_path << '\n';
  os << "prefix:                    " << conf.prefix << '\n';
  os << "dot_output:                " << conf.dot_output << '\n';
  os << "middle_gap_only:           " << conf.middle_gap_only << '\n';
  os << "prune:                     " << conf.prune << '\n';
  os << "should_output_all_matches: " << conf.should_output_all_matches << '\n';
  os << "p:                         " << conf.p << '\n';
  os << "p_ext:                     " << conf.p_ext << '\n';
  os << "log_fold:                  " << conf.log_fold << '\n';
  os << "h1_weight:                 " << conf.h1_weight << '\n';
  os << "print_lim:                 " << conf.print_lim << '\n';
  os << "lookup_k:                  " << conf.lookup_k << '\n';
  os << "max_aligns:                " << conf.max_aligns << "\n";
  os << "threads:                   " << conf.threads << '\n';
  os << "max_k:                     " << conf.max_k << '\n';
  os << "mem_limit:                 " << conf.mem_limit << '\n';
  os << "enable_smoothing:          " << conf.enable_smoothing << '\n';

  return os;
}


configuration parse_command_line_arguments(int argc, char const* argv[]) {
  configuration retval;

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
        parser, "path", "Background FASTA file, required for now.",
        {'b', "background"}, retval.bg_path, args::Options::Required);
    args::Positional<std::string> sig_path_(
        parser, "signal_path", "Signal FASTA file.", retval.sig_path,
        args::Options::Required);
    sf::args::value_flag dot_path_(
        parser, "path", "Compute Huddinge graph and output to dot file path.",
        {"dot"}, retval.dot_output);
    args::Flag gap_any_(parser, "gap_any",
                        "Allow gaps at any location, not just in the middle.",
                        {'a', "gap-at-any-location"});
    sf::args::value_flag p_sig_(
        parser, "p_value",
        "p value to use for signal to background comparison.",
        {'p', "p-signal"}, retval.p);
    sf::args::value_flag p_ext_(
        parser, "p_value",
        "p value to use for extension when background counts are zero.",
        {"p-ext"}, retval.p_ext);
    sf::args::value_flag log_fold_(
        parser, "fold_change",
        "Discard all mers with log fold change smaller than this.", {"lf"},
        retval.log_fold);
    sf::args::value_flag max_k_(parser, "length",
                                "Maximum mer length in [6, 24] range.", {"mk"},
                                retval.max_k);
    sf::args::value_flag lookup_k_(
        parser, "length",
        "Limit for lookup table-based k-mer counting in [5, max_k] range. Use "
        "zero to calculate automatically from available memory.",
        {"lookup-k"}, retval.lookup_k);
    sf::args::value_flag threads_(parser, "threads",
                                  "Number of threads to use.", {'t', "threads"},
                                  retval.threads);
    sf::args::value_flag prefix_(parser, "output_prefix",
                                 "Prefix for alignment output, If not given, "
                                 "no alignments will be output.",
                                 {"pref"}, retval.prefix);
    sf::args::value_flag print_lim_(
        parser, "max_s",
        "Maximum number of “best” seeds to output (0 -> all seeds).", {"max-s"},
        retval.print_lim);
    sf::args::value_flag max_aligns_(
        parser, "max_a",
        "Maximum number of alignments to output. (0 -> max_s).", {"max-a"},
        retval.max_aligns);
    args::Flag should_output_all_matches_(
        parser, "output_all_matches",
        "Output also matches that are substrings of other matches.",
        {"output-all-matches"});
    args::Flag disable_smoothing_(parser, "disable_smoothing",
                                  "Disable smoothing of counted mers.",
                                  {'s', "no-smoothing"});
    args::Flag enable_pruning_(
        parser, "enable pruning",
        "Enable pruning of extendable mers for partial counting.", {"pruning"});
    sf::args::value_flag mem_limit_(
        parser, "gigabytes",
        "Approximate memory limit for lookup tables in gigabytes.", {"mem"},
        retval.mem_limit);
    sf::args::value_flag h1_weight_(
        parser, "weight",
        "Relative impact of H1 neighbourhood enrichment on mer priority. (0 -> "
        "no impact, 1 -> 0.5 h1 neighbourhod 0.5 mer enrichment)",
        {"h1-weight"}, retval.h1_weight);

    // Parse and check.
    try {
      parser.ParseCLI(argc, argv);
    } catch (args::Help const&) {
      std::cerr << parser;
      std::exit(0);
    } catch (args::Completion const&) {
      std::cerr << parser;
      std::exit(0);
    } catch (sf::args::output_version const&) {
      std::cout << "seed_finder " << sf::version << '\n';
      std::exit(0);
    } catch (std::runtime_error const& err) {
      std::cerr << err.what() << '\n';
      std::exit(1);
    }

    retval.sig_path = args::get(sig_path_);
    if (gap_any_) retval.middle_gap_only = false;
    if (should_output_all_matches_) retval.should_output_all_matches = true;
    if (disable_smoothing_) retval.enable_smoothing = false;
    if (enable_pruning_) retval.prune = true;
  }

  if (retval.print_lim == 0) {
    retval.print_lim = ~retval.print_lim;
  }

  if (retval.max_aligns == 0) {
    retval.max_aligns = retval.print_lim;
  }

  if (retval.max_k < 6 || retval.max_k > 24) {
    std::cerr << "invalid value for maximum k: " << retval.max_k
              << ", should be in [6, 24] range." << std::endl;
    exit(1);
  }

  return retval;
}


struct reader_adapter_delegate : public sf::reader_adapter_delegate {
  bool should_report_errors_for_path(sf::reader_adapter&,
                                     std::string_view path) override {
    return true;
  }

  void found_first_read_with_unexpected_character(
      sf::reader_adapter&, std::string_view path,
      std::uint64_t lineno) override {
    std::cerr << "WARNING: Skipping reads with unexpected characters in "
              << path << "; first one on line " << lineno << ".\n";
  }

  void found_total_reads_with_unexpected_characters(
      sf::reader_adapter&, std::string_view path,
      std::uint64_t count) override {
    std::cerr << "WARNING: Skipped " << count << " reads in " << path << ".\n";
  }
};
}  // namespace


int main(int argc, char const* argv[]) {
  print_invocation(argc, argv);
  configuration conf{parse_command_line_arguments(argc, argv)};

#ifdef DEBUG
  std::cerr << conf;
#endif

  conf.mem_limit *= 1000;
  conf.mem_limit *= 1000;
  conf.mem_limit *= 1000;

  // Calculate lookup_k from the available memory if set to zero.
  if (conf.lookup_k == 0) {
    conf.lookup_k = 5;

    // Figure out how big lookup tables will fit in memory.
    sf::call_with_constant(
        conf.middle_gap_only, [&](auto const middle_gap_only) {
          while (sf::gapmer_count<middle_gap_only, max_gap>::lookup_bytes(
                     conf.lookup_k) < conf.mem_limit &&
                 conf.lookup_k <= 10) {
            ++conf.lookup_k;
          }
        });
  }

  if (conf.lookup_k < 5 || conf.lookup_k > conf.max_k) {
    std::cerr << "Invalid lookup_k value (" << conf.lookup_k << " not in [" << 5
              << ", " << conf.max_k << "]).";
    exit(1);
  }

#ifdef _OPENMP
  omp_set_num_threads(conf.threads);
#endif

  sf::packed_read_vector signal_reads;
  sf::packed_read_vector background_reads;

  // Process the input reads.
  {
    reader_adapter_delegate delegate;
    sf::reader_adapter_type reader{delegate};

    auto const process_path{
        [&](std::string const& path, sf::packed_read_vector& dst) {
          sf::reader_adapter_guard guard{reader};
          reader.read_from_path(path);
          while (reader.retrieve_next_read())
            dst.emplace_back(reader.read_buffer(), reader.read_length());
        }};

    process_path(conf.sig_path, signal_reads);
    process_path(conf.bg_path, background_reads);
  }

  // Prepare to output the alignments if needed.
  if (not conf.prefix.empty()) {
    std::filesystem::create_directory(conf.prefix);
  }

  // Run the algorithm. The parameters are std::bool_constants and hence can be
  // used as non-type template parameters.
  auto const run([&](auto const middle_gap_only,
                     auto const enable_smoothing) -> void {
    typedef sf::seed_finder<middle_gap_only, max_gap, enable_smoothing, false>
        seed_finder_type;
    typedef typename seed_finder_type::gapmer_type gapmer_type;

    seed_finder_type finder(signal_reads, background_reads, conf.p,
                            conf.log_fold, conf.max_k, conf.mem_limit,
                            conf.p_ext, conf.lookup_k, conf.prune);

    finder.find_seeds();
    if (!conf.dot_output.empty()) {
      sf::Dot_Writer::write_dot<gapmer_type>(conf.dot_output,
                                             finder.get_seeds(), conf.max_k);
    }

    auto sc{make_seed_clusterer<middle_gap_only, max_gap>(
        finder.get_seeds(), signal_reads, background_reads, conf.p_ext,
        conf.h1_weight, finder.x())};
    std::cout << "Seed\tcounts\tp\tpriority" << std::endl;
    for (size_t i = 0; i < conf.print_lim; ++i) {
      if (not sc.has_next()) {
        break;
      }
      if (conf.max_aligns == i) {
        conf.prefix = "";
      }
      sc.output_cluster(conf.prefix, conf.should_output_all_matches);
    }
  });

  // Convert runtime parameters to constants.
  sf::call_with_constant(conf.middle_gap_only, [&](auto const middle_gap_only) {
    sf::call_with_constant(conf.enable_smoothing,
                           [&](auto const enable_smoothing) {
                             run(middle_gap_only, enable_smoothing);
                           });
  });

  return 0;
}
