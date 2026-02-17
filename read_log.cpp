/*
 * Copyright (c) 2026 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <args.hxx>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <libbio/assert.hh>
#include <print>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>


namespace {

template <typename t_value>
struct field_value {
  t_value value{};
  std::int32_t field{};
};

template <typename t_value>
struct field_value_reader {
  void operator()(std::string const& name, std::string const& value,
                  field_value<t_value>& dst) {
    std::string_view sv{value};
    auto const comma_pos{sv.find(',')};
    if (std::string_view::npos == comma_pos)
      throw args::ValidationError("Comma not found in passed field value");

    auto const first{sv.data()};
    auto const last{first + comma_pos};
    auto const res{std::from_chars(first, last, dst.field)};
    if (not(std::errc{} == res.ec && res.ptr == last))
      throw args::ValidationError("Unable to parse field");

    if (dst.field <= 0) {
      throw args::ValidationError("Field index must be positive");
    }

    --dst.field;
    dst.value = sv.substr(1 + comma_pos);
  }
};


void split_tabs(std::string const& input, std::vector<std::string_view>& dst) {
  std::string_view sv{input};
  std::size_t start{};
  dst.clear();
  while (start < sv.size()) {
    auto const pos{sv.find('\t', start)};
    if (std::string_view::npos == pos) {
      dst.emplace_back(sv.substr(start));
      break;
    }

    dst.emplace_back(sv.substr(start, pos - start));
    start = pos + 1;
  }
}
}  // namespace


int main(int argc, char** argv) {
  std::ios_base::sync_with_stdio(false);

  args::ArgumentParser parser(
      "Read a TSV-formatted log file from stdin and filter lines according to "
      "the given options.");
  args::ValueFlagList<std::size_t> output_field_indices_(
      parser, "index", "Output only the given fields", {'f', "field"});
  args::ValueFlagList<field_value<std::string>, std::vector,
                      field_value_reader<std::string>>
      field_equals_(parser, "equals", "Given field equals the passed value",
                    {"equals"});
  args::ValueFlagList<field_value<std::string>, std::vector,
                      field_value_reader<std::string>>
      field_contains_(parser, "contains",
                      "Given field contains the passed value", {"contains"});

  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help const&) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError const& e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  } catch (args::ValidationError const& e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  std::string buffer;
  std::vector<std::size_t> output_field_indices{};
  std::vector<std::string_view> fields;
  std::size_t lineno{};

  auto const output([&] -> void {
    if (output_field_indices.empty()) {
      std::cout << buffer << '\n';
      return;
    }

    bool is_first{true};
    for (auto const idx : output_field_indices) {
      if (not is_first) std::cout << '\t';
      is_first = false;
      std::cout << fields[idx];
    }
    std::cout << '\n';
  });

  std::size_t expected_field_count{};
  while (std::getline(std::cin, buffer)) {
    ++lineno;

    if (buffer.empty()) {
      std::cerr << "ERROR: Line " << lineno << " was empty.\n";
      return 1;
    }

    // Skip the header if there was one.
    if ('#' == buffer.front()) continue;

    split_tabs(buffer, fields);
    expected_field_count = fields.size();

    for (auto const idx : args::get(output_field_indices_)) {
      if (not idx) {
        std::cerr << "ERROR: Output field indices must be positive, got " << idx
                  << ".\n";
        return 1;
      }

      auto const idx_{idx - 1};
      if (expected_field_count <= idx_) {
        std::cerr << "ERROR: Got output field index " << idx
                  << " but expected field count is" << expected_field_count
                  << ".\n";
        return 1;
      }
      output_field_indices.push_back(idx_);
    }

    output();
    break;
  }

  // Check if a header was found.
  if (0 == lineno) return 0;

  auto const check_tested_field_indices{
      [expected_field_count](auto& args) -> void {
        for (auto const& fv : args::get(args)) {
          std::size_t field_idx(fv.field);
          if (expected_field_count <= field_idx) {
            std::cerr << "ERROR: Got output field index " << field_idx
                      << " but expected field count is" << expected_field_count
                      << ".\n";
            std::exit(1);
          }
        }
      }};

  check_tested_field_indices(field_equals_);
  check_tested_field_indices(field_contains_);

  auto const check_fields([&](auto& args, auto const& test) -> bool {
    for (auto const& fv : args::get(args)) {
      std::size_t field_idx(fv.field);
      libbio_assert_lt(field_idx, fields.size());
      if (not test(fv.value, fields[field_idx])) return false;
    }

    return true;
  });

  while (std::getline(std::cin, buffer)) {
    ++lineno;
    split_tabs(buffer, fields);

    if (fields.size() != expected_field_count) {
      std::print(std::cerr,
                 "WARNING: Unexpected field count on line {} (expected: {} "
                 "actual: {})\n",
                 lineno, expected_field_count, fields.size());
      continue;
    }

    if (not check_fields(field_equals_, std::equal_to{})) continue;
    if (not check_fields(field_contains_,
                         [](auto const& passed, auto const& actual) {
                           return actual.contains(passed);
                         }))
      continue;

    output();
  }

  return 0;
}
