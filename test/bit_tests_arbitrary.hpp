/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

#include "../include/bits.hpp"
#include "../include/string_buffer.hpp"
#include "../include/util.hpp"
#include "nucleotide.hpp"
#include "rapidcheck/Assertions.h"
#include "test.hpp"

namespace {

struct short_dna_string {
  std::string value;
};


std::string string_from_packed(uint64_t word, uint8_t const count) {
  std::string retval;
  retval.resize(count, 0);
  for (uint8_t ii{}; ii < count; ++ii) {
    auto const packed(word & 0x3);
    retval[count - ii - 1] = sf::v_to_nuc[packed];
    word >>= 2;
  }
  return retval;
}
}  // namespace


namespace rc {

template <>
struct Arbitrary<short_dna_string> {
  static Gen<short_dna_string> arbitrary() {
    return gen::mapcat(
        gen::inRange<std::size_t>(1, 9), [](std::size_t const length) {
          return gen::construct<short_dna_string>(gen::container<std::string>(
              length, gen::arbitrary<nucleotide>()));
        });
  }
};
}  // namespace rc


namespace sf {

SF_RC_TEMPLATE_TEST(bit_arbitrary_byteswap, ByteswapWorksAsExpected,
                    (TypeParam const value), uint16_t, uint32_t, uint64_t, int16_t, int32_t, int64_t) {
  constexpr static bool const stdlib_byteswap_available{
#if defined(__cpp_lib_byteswap)
      true
#else
      false
#endif
  };

  if constexpr (stdlib_byteswap_available) {
    auto const expected(std::byteswap(value));
    auto const res(sf::bits::detail::byteswap_generic(value));
    RC_ASSERT(res == expected);
  } else {
    GTEST_SKIP()
        << "Test skipped: Unable to test the byteswap implementation when "
           "the stdlib implementation is not available.";
  }
}


SF_RC_TEMPLATE_TEST(bit_arbitrary_pext, PEXTWorksAsExpected,
                    (TypeParam const value, TypeParam const mask), uint32_t,
                    uint64_t) {
  if constexpr (sf::bits::detail::pext_intrinsic_available_for_type_v<
                    TypeParam>) {
    auto const res(sf::bits::detail::pext_generic(value, mask));
    auto const expected(sf::bits::detail::pext_intrinsic(value, mask));
    RC_ASSERT(res == expected);
  } else {
    GTEST_SKIP() << "Test skipped: Unable to test the PEXT implementation when "
                    "the intrinsic is not available.";
  }
}


RC_GTEST_PROP(bit_arbitrary, ReadMultipleCharacters,
              (short_dna_string const &dna_str)) {
  typedef string_buffer<uint64_t> string_buffer_type;
  auto const length(dna_str.value.size());
  SF_RC_TAG(length);
  string_buffer_type sb{dna_str.value};
  auto word(bits::detail::read_multiple_dna_characters_generic(*sb.data()));
  word >>= 2 * (8 - length);
  auto res(string_from_packed(word, length));
  RC_ASSERT(dna_str.value == res);
}
}  // namespace sf
