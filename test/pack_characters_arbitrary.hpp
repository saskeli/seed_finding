/*
 * Copyright (c) 2025-2026 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include <cstdint>
#include <range/v3/view/all.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/view/zip.hpp>
#include <string>
#include <vector>

#include "../include/pack_characters.hpp"
#include "gtest/gtest.h"
#include "nucleotide.hpp"
#include "test.hpp"


namespace {
template <sf::dna_alphabet t_alphabet>
struct pack_characters_arbitrary_input {
  typedef std::vector<std::vector<nucleotide_tpl<t_alphabet>>> text_type;

  constexpr static sf::dna_alphabet alphabet{t_alphabet};
  text_type text;
};
}  // namespace


namespace rc {
template <sf::dna_alphabet t_alphabet>
struct Arbitrary<pack_characters_arbitrary_input<t_alphabet>> {
  static Gen<pack_characters_arbitrary_input<t_alphabet>> arbitrary() {
    typedef pack_characters_arbitrary_input<t_alphabet> return_type;
    typedef typename return_type::text_type text_type;

    return gen::map(gen::arbitrary<text_type>(),
                    [](auto&& text) { return return_type{std::move(text)}; });
  }
};
}  // namespace rc


SF_RC_TEMPLATE_TEST(pack_characters_arbitrary, PackAndUnpack,
                    (TypeParam const& input),
                    pack_characters_arbitrary_input<sf::dna_alphabet::dna4>,
                    pack_characters_arbitrary_input<sf::dna_alphabet::dna16>) {
  {
    // Determine the sum of the lengths of the input texts.
    auto const total_length{[&] {
      std::uint64_t retval{};
      for (auto const& vec : input.text) retval += vec.size();
      return retval;
    }()};
    RC_TAG(input.text.size(), total_length);
  }

  std::vector<std::uint64_t> packed;
  std::string unpacked;

  std::string buffer;
  std::uint64_t pos{};

  // Pack the characters.
  for (auto const& vec : input.text) {
    buffer.clear();
    for (auto const nt : vec) buffer.push_back(char(nt));

    pos = sf::pack_characters_<input.alphabet>(buffer, packed, pos);
    EXPECT_NE(UINT64_MAX, pos);
  }

  // Unpack.
  sf::unpack_characters_<input.alphabet>(packed, pos, unpacked);

  // Since the input consists of multiple vectors of characters,
  // we flatten the enclosing range so that we can compare the
  // unpacked value (above) to the input character by character.
  auto rng{ranges::views::for_each(
      input.text, [](auto const& vec) { return ranges::views::all(vec); })};

  // Compare character by character.
  for (auto const& [lhs, rhs] : ranges::views::zip(rng, unpacked))
    EXPECT_EQ(lhs, rhs);
}
