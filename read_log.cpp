/*
 * Copyright (c) 2026 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <args.hxx>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <iostream>
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
      "the given "
      "options.");
  args::ValueFlagList<field_value<std::string>, std::vector,
                      field_value_reader<std::string>>
      field_equals_(parser, "equals", "Given field equals the passed value",
                    {"equals"});

  try {
    parser.ParseCLI(argc, argv);
  } catch (args::Help) {
    std::cout << parser;
    return 0;
  } catch (args::ParseError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  } catch (args::ValidationError e) {
    std::cerr << e.what() << std::endl;
    std::cerr << parser;
    return 1;
  }

  std::string buffer;
  std::vector<std::string_view> fields;
  std::size_t lineno{};

  auto const check_fields_eq([&] -> bool {
    for (auto const& eq : args::get(field_equals_)) {
      std::size_t field_idx(eq.field - 1); // Checked earlier.
      if (fields.size() <= field_idx) {
        std::print(std::cerr, "WARNING: Field count on line {} was {}.\n",
                   lineno, fields.size());
        return false;
      }

      if (eq.value != fields[field_idx]) return false;
    }

    return true;
  });

  while (std::getline(std::cin, buffer)) {
    ++lineno;
    split_tabs(buffer, fields);

    if (not check_fields_eq()) continue;

    std::cout << buffer << '\n';
  }

  return 0;
}
