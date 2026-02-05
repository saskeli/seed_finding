/*
 * Copyright (c) 2026 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <args.hxx>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "include/args.hpp"
#include "include/configuration.hpp"
#include "include/count_base.hpp"
#include "include/gapmer.hpp"
#include "include/packed_read.hpp"
#include "include/reader_adapter.hpp"
#include "include/version.hpp"

namespace {

typedef sf::gapmer<true, 15> gapmer_type;


struct configuration {
  std::string input_path;
  std::uint16_t kk{};
};


class counter final : public sf::count_base<gapmer_type> {
  virtual uint64_t offset(uint8_t gap_start, uint8_t gap_length) const override;

  virtual void increment_signal_count(gapmer_type gg) override;
  virtual inline void increment_signal_count_gapped(gapmer_type gg,
                                                    uint64_t off) override;

  virtual void increment_background_count(gapmer_type gg) override;
  virtual void increment_background_count_gapped(gapmer_type gg,
                                                 uint64_t off) override;
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
  args::ArgumentParser parser("List gapmers of given length in input FASTA/Q.");

  args::HelpFlag help_(parser, "help", "Display this help.", {'h', "help"});
  sf::args::version_flag version_(parser, "version",
                                  "Output the version number and exit.",
                                  {'V', "version"});
  args::CompletionFlag completion_(parser, {"complete"});

  args::Positional kk_(parser, "lookup_k", "Number of defined characters",
                       conf.kk, args::Options::Required);
  args::Positional input_path_(parser, "signal_path", "Signal FASTA file.",
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


uint64_t counter::offset(uint8_t gap_start, uint8_t gap_length) const {
  return 0;
}


void counter::increment_signal_count(gapmer_type gg) {
  std::cout << gg << '\n';
}


void counter::increment_signal_count_gapped(gapmer_type gg, uint64_t off) {
  increment_signal_count(gg);
}


void counter::increment_background_count(gapmer_type gg) {
  throw std::runtime_error("Should not be called");
}


void counter::increment_background_count_gapped(gapmer_type gg, uint64_t off) {
  throw std::runtime_error("Should not be called");
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
}  // namespace


int main(int argc, char** argv) {
  configuration conf;
  parse_arguments(argc, argv, conf);

  sf::packed_read_vector reads;
  read_input(conf.input_path, reads);

  counter cc;
  cc.count_mers(reads, conf.kk);

  return 0;
}
