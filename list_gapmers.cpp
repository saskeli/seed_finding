/*
 * Copyright (c) 2026 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <omp.h>

#include <args.hxx>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>

#include "include/args.hpp"
#include "include/configuration.hpp"
#include "include/count_base.hpp"
#include "include/gapmer.hpp"
#include "include/gapmer_count.hpp"
#include "include/packed_read.hpp"
#include "include/partial_count.hpp"
#include "include/reader_adapter.hpp"
#include "include/version.hpp"

namespace {

typedef sf::gapmer<true, 15> gapmer_type;


enum class action_type { list, count };


enum class counting_method_type { map, gapmer_count, partial_count };


struct configuration {
  std::string input_path;
  std::uint16_t kk{};
  action_type action{action_type::list};
  counting_method_type counting_method{counting_method_type::map};
};


class count_base : public sf::count_base<gapmer_type> {
  virtual uint64_t offset(uint8_t gap_start, uint8_t gap_length) const override;
  virtual inline void increment_signal_count_gapped(
      gapmer_type gg, uint64_t off, counting_context const&) override;

  virtual void increment_background_count(gapmer_type gg,
                                          counting_context const&) override;
  virtual void increment_background_count_gapped(
      gapmer_type gg, uint64_t off, counting_context const&) override;
};


template <typename t_count>
class count_wrapper; // Fwd.

template <typename t_count>
void output_all(count_wrapper<t_count> const&); // Fwd.

template <typename t_count>
class count_wrapper final : public count_base {
  friend void output_all<>(count_wrapper const&);

  typedef t_count count_type;
  typedef std::set<gapmer_type> gapmer_set;

  gapmer_set m_seen_gapmers;
  count_type m_count;

  virtual uint64_t offset(uint8_t gap_start, uint8_t gap_length) const override;
  virtual void increment_signal_count(gapmer_type gg,
                                      counting_context const&) override;
  virtual void increment_signal_count_gapped(gapmer_type gg, uint64_t off,
                                             counting_context const&) override;

 public:
  template <typename... t_args>
  explicit count_wrapper(t_args&&... args)
      : m_count(std::forward<t_args>(args)...) {}
};


typedef count_wrapper<sf::gapmer_count<gapmer_type>> gapmer_count;
typedef count_wrapper<sf::partial_count<gapmer_type, true>> partial_count;


class count_stl_map final : public count_base {
  friend void output_all(count_stl_map const&);

  typedef std::map<gapmer_type, std::uint32_t> count_map;

  count_map m_counts;

  virtual void increment_signal_count(gapmer_type gg,
                                      counting_context const&) override;
};


class list_gapmers final : public count_base {
  virtual void increment_signal_count(gapmer_type gg,
                                      counting_context const&) override;
};


struct reader_adapter_delegate : public sf::reader_adapter_delegate {
  bool should_report_errors_for_path(sf::reader_adapter&,
                                     std::string_view path) override;
  void found_first_read_with_unexpected_character(
      sf::reader_adapter&, std::string_view path,
      std::uint64_t lineno) override;
  void found_total_reads_with_unexpected_characters(
      sf::reader_adapter&, std::string_view path, std::uint64_t count) override;
};


void parse_arguments(int const argc, char const* const* const argv,
                     configuration& conf) {
  std::unordered_map<std::string, counting_method_type> const counting_methods{
      {"stl-map", counting_method_type::map},
      {"gapmer-count", counting_method_type::gapmer_count},
      {"partial-count", counting_method_type::partial_count}};

  args::ArgumentParser parser("List gapmers of given length in input FASTA/Q.");

  args::HelpFlag help_(parser, "help", "Display this help.", {'h', "help"});
  sf::args::version_flag version_(parser, "version",
                                  "Output the version number and exit.",
                                  {'V', "version"});
  args::CompletionFlag completion_(parser, {"complete"});

  args::MapFlag<std::string, counting_method_type> count_(
      parser, "count", "Count gapmers instead of listing", {"count"},
      counting_methods, counting_method_type::map, args::Options::Single);

  args::Positional kk_(parser, "lookup_k", "Number of defined characters",
                       conf.kk, args::Options::Required);
  args::Positional input_path_(parser, "input_path", "Input FASTA/Q file",
                               conf.input_path, args::Options::Required);

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

  conf.kk = args::get(kk_);
  conf.input_path = args::get(input_path_);

  if (count_) {
    conf.action = action_type::count;
    conf.counting_method = args::get(count_);
  }
}


void read_input(std::string const& path, sf::packed_read_vector& dst) {
  reader_adapter_delegate delegate;
  sf::reader_adapter_type reader{delegate};

  {
    sf::reader_adapter_guard guard{reader};
    reader.read_from_path(path);
    while (reader.retrieve_next_read())
      dst.emplace_back(reader.read_buffer(), reader.read_length());
  }
}


uint64_t count_base::offset(uint8_t gap_start, uint8_t gap_length) const {
  return 0;
}


void count_base::increment_signal_count_gapped(gapmer_type gg, uint64_t off,
                                               counting_context const& ctx) {
  increment_signal_count(gg, ctx);
}


void count_base::increment_background_count(gapmer_type gg,
                                            counting_context const&) {
  throw std::runtime_error("Should not be called");
}


void count_base::increment_background_count_gapped(gapmer_type gg, uint64_t off,
                                                   counting_context const&) {
  throw std::runtime_error("Should not be called");
}


void list_gapmers::increment_signal_count(gapmer_type gg,
                                          counting_context const& ctx) {
  auto const read_idx{ctx.read_index / 2};
  auto const is_reverse_complement{ctx.read_index % 2};
  std::cout << read_idx << '\t' << is_reverse_complement << '\t'
            << ctx.lhs_position << '\t' << gg << '\n';
}


void count_stl_map::increment_signal_count(gapmer_type gg,
                                           counting_context const& ctx) {
  ++m_counts[gg];
}


template <typename t_count>
uint64_t count_wrapper<t_count>::offset(uint8_t gap_start,
                                        uint8_t gap_length) const {
  return m_count.offset(gap_start, gap_length);
}


template <typename t_count>
void count_wrapper<t_count>::increment_signal_count(
    gapmer_type gg, counting_context const& ctx) {
  m_seen_gapmers.insert(gg);
  m_count.increment_signal_count(gg, ctx);
}


template <typename t_count>
void count_wrapper<t_count>::increment_signal_count_gapped(
    gapmer_type gg, uint64_t off, counting_context const& ctx) {
  m_seen_gapmers.insert(gg);
  m_count.increment_signal_count_gapped(gg, off, ctx);
}


bool reader_adapter_delegate::should_report_errors_for_path(
    sf::reader_adapter&, std::string_view path) {
  return true;
}


void reader_adapter_delegate::found_first_read_with_unexpected_character(
    sf::reader_adapter&, std::string_view path, std::uint64_t lineno) {
  std::cerr << "WARNING: Skipping reads with unexpected characters in " << path
            << "; first one on line " << lineno << ".\n";
}


void reader_adapter_delegate::found_total_reads_with_unexpected_characters(
    sf::reader_adapter&, std::string_view path, std::uint64_t count) {
  std::cerr << "WARNING: Skipped " << count << " reads in " << path << ".\n";
}


void output_all(count_stl_map const& map) {
  for (auto const& kv : map.m_counts) {
    std::cout << kv.first << '\t' << kv.second << '\n';
  }
}


template <typename t_count>
void output_all(count_wrapper<t_count> const& map) {
  for (auto const& gg : map.m_seen_gapmers) {
    auto const cc{map.m_count.count(gg)};
    std::cout << gg << '\t' << cc.signal_count << '\n';
  }
}


void do_list(sf::packed_read_vector const& reads, configuration const& conf) {
  std::cout << "read_index\tis_reverse_complement\tlhs_position\tsequence\n";
  list_gapmers lg;
  lg.count_mers(reads, conf.kk);
}


void do_count(sf::packed_read_vector const& reads, configuration const& conf) {
  std::cout << "sequence\tcount\n";
  switch (conf.counting_method) {
    case counting_method_type::map: {
      count_stl_map cc;
      cc.count_mers(reads, conf.kk);
      output_all(cc);
      break;
    }

    case counting_method_type::gapmer_count: {
      gapmer_count cc(conf.kk);
      cc.count_mers(reads, conf.kk);
      output_all(cc);
      break;
    }

    case counting_method_type::partial_count: {
      partial_count cc{};
      cc.count_mers(reads, conf.kk);
      output_all(cc);
      break;
    };
  }
}

}  // namespace


int main(int argc, char** argv) {
  omp_set_num_threads(1);

  configuration conf;
  parse_arguments(argc, argv, conf);

  sf::packed_read_vector reads;
  read_input(conf.input_path, reads);

  switch (conf.action) {
    case action_type::list:
      do_list(reads, conf);
      break;
    case action_type::count:
      do_count(reads, conf);
      break;
  };

  return 0;
}
