#pragma once

#include <cstddef>
#include <cstdint>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/zip.hpp>
#include <span>
#include <utility>
#include <vector>

#include "../include/packed_character_iteration.hpp"
#include "nucleotide.hpp"
#include "test.hpp"

// FIXME: Some of the types used in the tests need to have typedefs for the
// associated types since the preprocessor has difficulties with commas as part
// of the type names. Some of the typedefs have very long names in order to
// encode the name of the test they have to do with. Separating the tests into
// multiple translation units could help with this.

namespace {

typedef std::vector<packed_nucleotide> packed_nucleotide_vector;


struct packed_character_iteration_test_input {
  packed_nucleotide_vector text;
  std::size_t start_pos{};
};


struct packed_character_pair_iteration_test_input {
  packed_nucleotide_vector text;
  std::size_t lhs{};
  std::size_t rhs{};
};


std::vector<std::uint64_t> packed_text(packed_nucleotide_vector const& src) {
  sf::packed_word_vector retval((src.size() + 31U) / 32U, 0);
  std::uint64_t pos{};
  for (auto const nt : src) {
    auto const word_idx{pos / 32U};
    auto const chr_idx{pos % 32U};
    std::uint64_t nt_(nt.value);
    nt_ <<= 62U - 2U * chr_idx;
    retval.at(word_idx) |= nt_; // Bounds check for extra safety.
    ++pos;
  }

  return retval;
}


template <typename t_cb>
void iterate_packed_charaters(
    packed_character_pair_iteration_test_input const& input, t_cb&& cb) {
  RC_ASSERT(0 == input.rhs || input.lhs < input.rhs);
  std::size_t const count{input.text.size() - input.rhs};
  for (std::size_t ii{}; ii < count; ++ii)
    cb(input.text[input.lhs + ii], input.text[input.rhs + ii]);
}


std::ostream& operator<<(std::ostream& os,
                         packed_character_iteration_test_input const& input) {
  os << "text: ";
  for (auto const cc : input.text) os << nucleotide(cc);
  os << " start_pos: " << input.start_pos;
  return os;
}


std::ostream& operator<<(
    std::ostream& os, packed_character_pair_iteration_test_input const& input) {
  os << "text: ";
  for (auto const cc : input.text) os << nucleotide(cc);
  os << " lhs: " << input.lhs << " rhs: " << input.rhs;
  return os;
}
}  // namespace


namespace rc {

template <>
struct Arbitrary<packed_character_iteration_test_input> {
  static Gen<packed_character_iteration_test_input> arbitrary() {
    return gen::mapcat(
        gen::arbitrary<std::vector<packed_nucleotide>>(),
        [](std::vector<packed_nucleotide>&& text) {
          std::size_t const length{
              text.size()}; // Determine the length before moving below.
          return gen::construct<packed_character_iteration_test_input>(
              gen::just(std::move(text)),
              length ? gen::inRange(std::size_t{}, length)
                     : gen::just(length) // Half-open range.
          );
        });
  }
};


template <>
struct Arbitrary<packed_character_pair_iteration_test_input> {
  static Gen<packed_character_pair_iteration_test_input> arbitrary() {
    return gen::mapcat(
        gen::arbitrary<std::vector<packed_nucleotide>>(),
        [](std::vector<packed_nucleotide>&& text) {
          std::size_t const length{
              text.size()}; // Determine the length before moving below.
          if (!length)
            return gen::just(packed_character_pair_iteration_test_input{});

          return gen::mapcat(
              1U < length ? gen::inRange(std::size_t{1U}, length)
                          : gen::just(std::size_t{1}),
              [text = std::move(text)](auto const rhs) { // Half-open range.
                return gen::construct<
                    packed_character_pair_iteration_test_input>(
                    gen::just(std::move(text)),
                    gen::inRange(std::size_t{}, rhs), // Half-open range.
                    gen::just(rhs));
              });
        });
  }
};
}  // namespace rc


namespace {

RC_GTEST_PROP(packed_character_iteration_arbitrary, IterateSingle,
              (packed_character_iteration_test_input const& input)) {
  SF_RC_TAG(input.text.size());

  std::vector<nucleotide> expected;
  std::vector<nucleotide> actual;

  {
    auto const span(std::span(input.text.begin(), input.text.end())
                        .subspan(input.start_pos));
    for (auto nt : span) expected.emplace_back(nt);
  }

  auto const packed_text_{packed_text(input.text)};
  sf::iterate_packed_characters(
      packed_text_, input.text.size(), input.start_pos,
      [&](auto const cc) { actual.emplace_back(packed_nucleotide(cc)); });

  RC_ASSERT(expected == actual);
}


RC_GTEST_PROP(packed_character_iteration_arbitrary, IteratePairs,
              (packed_character_pair_iteration_test_input const& input)) {
  SF_RC_TAG(input.text.size());

  std::vector<nucleotide> expected;
  std::vector<nucleotide> actual;

  iterate_packed_charaters(input, [&](auto const lhsc, auto const rhsc) {
    expected.emplace_back(lhsc);
    expected.emplace_back(rhsc);
  });

  auto const packed_text_{packed_text(input.text)};
  sf::iterate_packed_character_pairs(
      packed_text_, input.text.size(), input.lhs, input.rhs,
      [&](auto const lhsc, auto const rhsc) {
        actual.emplace_back(packed_nucleotide(lhsc));
        actual.emplace_back(packed_nucleotide(rhsc));
      });

  RC_LOG() << "Character pairs (expected length " << expected.size()
           << ", actual length " << actual.size() << "):\n";
  for (auto const& [idx, tup] :
       ranges::views::zip(expected, actual) | ranges::views::enumerate) {
    auto const& [ec, ac] = tup;
    RC_LOG() << idx << ":\t" << ec << ' ' << ac;
    if (ec != ac) RC_LOG() << " mismatch";
    RC_LOG() << '\n';
  }

  RC_ASSERT(expected == actual);
}
}  // namespace
