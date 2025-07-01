/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <climits>
#include <concepts>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <span>

#if defined(__linux__)
#include <byteswap.h>
#endif

#if defined(__i386__) || defined(__x86_64__)
#include <immintrin.h>
#endif


namespace sf::bits::detail {

#if defined(__BMI2__)  // PEXT is part of BMI2.

template <typename t_type>
constexpr static inline bool const pext_intrinsic_available_for_type_v{
    std::is_same_v<t_type, std::uint32_t> ||
    std::is_same_v<t_type, std::uint64_t>};

inline uint32_t pext_intrinsic(uint32_t source, uint32_t mask) {
  return _pext_u32(source, mask);
}

inline uint64_t pext_intrinsic(uint64_t source, uint64_t mask) {
  return _pext_u64(source, mask);
}

#else

template <typename>
constexpr static inline bool const pext_intrinsic_available_for_type_v{false};

inline uint64_t pext_intrinsic(uint64_t, uint64_t) {
  return 0;
}  // FIXME: I don’t remember how to use if constexpr in such a way that the
   // discarded statement is not checked. The compiler currently checks both
   // branches in sf::bits::pext().

#endif


#if defined(__linux__)
inline uint16_t byteswap_linux(uint16_t source) { return bswap_16(source); }
inline uint32_t byteswap_linux(uint32_t source) { return bswap_32(source); }
inline uint64_t byteswap_linux(uint64_t source) { return bswap_64(source); }
#endif


template <typename t_type>
constexpr inline t_type byteswap_generic(t_type val) {
  auto val_(std::bit_cast<std::array<std::byte, sizeof(t_type)>>(val));
  std::reverse(val_.begin(), val_.end());
  return std::bit_cast<t_type>(val_);
}


template <std::unsigned_integral t_unsigned>
constexpr inline t_unsigned pext_generic(t_unsigned src, t_unsigned mask) {
  // A linear in the number of runs of set bits in the mask implementation of
  // PEXT.

  t_unsigned retval{};
  t_unsigned dst_idx{};
  t_unsigned clear_mask{};

  clear_mask = ~clear_mask;

  // As of clang-format 20, the tool breaks the formatting below unless
  // AlignTrailingComments's Kind is set to Always.
  // clang-format off
  while (mask) {
    t_unsigned const trailing_zeros(  // Find the next set bit in the mask.
        std::countr_zero(mask));      //
    src >>= trailing_zeros;           // Shift the source and the mask to
                                      // this position.
    mask >>= trailing_zeros;          //
    retval |= src << dst_idx;         // Extract.

    t_unsigned const trailing_ones(   // Find the next unset bit in the mask.
        std::countr_one(mask));       //
    src >>= trailing_ones;            // Skip to the end of the current run.
    mask >>= trailing_ones;           //
    dst_idx += trailing_ones;         // Set up the next target position.
    clear_mask <<= trailing_ones;     // Clear the bits outside the current run.
    retval &= ~clear_mask;
  }
  // clang-format on

  return retval;
}
}  // namespace sf::bits::detail


namespace sf::bits {

template <std::unsigned_integral t_unsigned>
constexpr inline t_unsigned byteswap(t_unsigned source) {
  if consteval {
#if defined(__cpp_lib_byteswap)
    return std::byteswap(source);
#else
    return detail::byteswap_generic(source);
#endif // defined(__cpp_lib_byteswap)
  } else {
#if defined(__cpp_lib_byteswap)
    return std::byteswap(source);
#else
#if defined(__linux__)
    if constexpr (2 <= sizeof(t_unsigned) && sizeof(t_unsigned) <= 8)
      return detail::byteswap_linux(source);
#endif // defined(__linux__)
    return detail::byteswap_generic(source);
#endif // defined(__cpp_lib_byteswap)
  }
}


template <std::unsigned_integral t_unsigned>
constexpr inline t_unsigned pext(t_unsigned source, t_unsigned mask) {
  if consteval {
    return detail::pext_generic(source, mask);
  } else {
    if constexpr (detail::pext_intrinsic_available_for_type_v<t_unsigned>)
      return detail::pext_intrinsic(source, mask);
    else
      return detail::pext_generic(source, mask);
  }
}


template <std::unsigned_integral t_value, std::size_t t_n,
          std::unsigned_integral t_count>
constexpr inline void shift_left(std::span<t_value, t_n> span,
                                 t_count const count) {
  if (!count) return;

  auto const value_bits{sizeof(t_value) * CHAR_BIT};
  assert(count < value_bits);  // The general case has not been implemented.

  auto const lower_mask([&]() constexpr {
    t_value mask{};
    mask = ~mask;
    mask >>= value_bits - count;
    return mask;
  }());
  auto const higher_mask(~lower_mask);

  t_value lower{};
  for (auto &word : span) {
    auto word_(std::rotl(word, count));
    word = word_ & higher_mask;
    word |= lower;
    lower = word_ & lower_mask;
  }
}

template <std::unsigned_integral t_value, std::size_t t_n,
          std::unsigned_integral t_count>
constexpr inline void shift_right(std::span<t_value, t_n> span,
                                  t_count const count) {
  if (!count) return;

  auto const value_bits{sizeof(t_value) * CHAR_BIT};
  assert(count < value_bits);  // The general case has not been implemented.

  auto const higher_mask([&]() constexpr {
    t_value mask{};
    mask = ~mask;
    mask <<= value_bits - count;
    return mask;
  }());
  auto const lower_mask(~higher_mask);

  t_value higher{};
  for (auto it(span.rbegin()); it != span.rend(); ++it) {
    auto &word(*it);
    auto word_(std::rotr(word, count));
    word = word_ & lower_mask;
    word |= higher;
    higher = word_ & higher_mask;
  }
}
}  // namespace sf::bits


namespace sf::bits::detail {

constexpr inline uint64_t const read_multiple_dna_characters_xor_mask{
    UINT64_C(0b0101010101010101)};


inline uint64_t read_multiple_dna_characters_generic(uint64_t word) {
  // Both Clang and GCC generate some 11 instructions for 64-bit ARM for this
  // function (excluding ret).

  // As of clang-format 20, the tool breaks the formatting below unless
  // AlignTrailingComments's Kind is set to Always.
  // clang-format off
  word = byteswap(word);                // Change the order to what
                                        // gapmer expects.
  word >>= 1;                           // Shift s.t. the relevant
                                        // bits are in the beginning
                                        // of each byte.
  word &= UINT64_C(0x303030303030303);  // Clear everything else.
  word |= word >> 6;                    // Group pairs of characters.
  word |= word >> 12;                   // Group quads of characters.
  word &= UINT64_C(0xFF000000FF);       // Clear everything else.
  word |= word >> 24;                   // Group the eight characters.
  word ^= (word >> 1) &                 // Fix the bit representation.
          read_multiple_dna_characters_xor_mask;
  word &= 0xFFFF;                       // Clear everything else.
  return word;
  // clang-format on
}


inline uint64_t read_multiple_dna_characters_pext(uint64_t const word) {
  constexpr uint64_t const pext_mask{UINT64_C(0x0606060606060606)};

  uint64_t retval{bits::byteswap(word)};
  retval = bits::pext(retval, pext_mask);
  retval ^= (retval >> 1) & read_multiple_dna_characters_xor_mask;
  return retval;
}
}  // namespace sf::bits::detail


namespace sf::bits {

/// Reads eight characters of our extended DNA alphabet ([ACGT.]) from
/// the given word and places them in reverse order and 2-bit encoded
/// to the return value.
inline uint64_t read_multiple_dna_characters(uint64_t const word) {
  if consteval {
    return detail::read_multiple_dna_characters_generic(word);
  } else {
    if constexpr (detail::pext_intrinsic_available_for_type_v<uint64_t>)
      return detail::read_multiple_dna_characters_pext(word);
    else
      return detail::read_multiple_dna_characters_generic(word);
  }
}
}  // namespace sf::bits
