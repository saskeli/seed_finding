/*
 * Copyright (c) 2025-2026 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <array>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>


namespace sf {
enum class dna_alphabet { dna4, dna16 };
typedef std::span<std::uint64_t const> packed_word_span;
typedef std::vector<std::uint64_t> packed_word_vector;


constexpr static std::array dna4_characters{'A', 'C', 'G', 'T'};

constexpr static std::array dna16_characters{'n', 'A', 'C', 'M', 'G', 'R',
                                             'S', 'V', 'T', 'W', 'Y', 'H',
                                             'K', 'D', 'B', 'N'};

static_assert(sizeof(dna16_characters) == 16);


template <dna_alphabet t_alphabet>
struct alphabet_size {};

template <>
struct alphabet_size<dna_alphabet::dna4> {
  constexpr static uint8_t value{dna4_characters.size()};
};

template <>
struct alphabet_size<dna_alphabet::dna16> {
  constexpr static uint8_t value{dna16_characters.size()};
};

template <dna_alphabet t_alphabet>
constexpr inline uint8_t alphabet_size_v{alphabet_size<t_alphabet>::value};


// Function templates.
template <dna_alphabet t_alphabet, bool t_is_lenient = false>
std::uint64_t pack_characters_(std::string_view sv, packed_word_vector& dst,
                               std::uint64_t dst_pos = 0);

template <dna_alphabet t_alphabet>
void unpack_characters_(packed_word_span const& src, std::uint64_t length,
                        std::string& dst);

template <dna_alphabet t_alphabet>
inline void unpack_characters_(packed_word_vector const& src,
                               std::uint64_t length, std::string& dst) {
  unpack_characters_<t_alphabet>(packed_word_span{src.data(), src.size()},
                                 length, dst);
}

// Forward declarations of instantiations.
extern template std::uint64_t pack_characters_<dna_alphabet::dna4>(
    std::string_view sv, packed_word_vector& dst, std::uint64_t dst_pos);

extern template std::uint64_t pack_characters_<dna_alphabet::dna4, true>(
    std::string_view sv, packed_word_vector& dst, std::uint64_t dst_pos);

extern template std::uint64_t pack_characters_<dna_alphabet::dna16>(
    std::string_view sv, packed_word_vector& dst, std::uint64_t dst_pos);

extern template void unpack_characters_<dna_alphabet::dna4>(
    packed_word_span const& src, std::uint64_t length, std::string& dst);

extern template void unpack_characters_<dna_alphabet::dna16>(
    packed_word_span const& src, std::uint64_t length, std::string& dst);


// DNA4
// The non-lenient version returns UINT64_MAX if an unexpected character is encountered.
[[nodiscard]] inline std::uint64_t pack_characters(std::string_view sv,
                                                   packed_word_vector& dst,
                                                   std::uint64_t dst_pos = 0) {
  return pack_characters_<dna_alphabet::dna4>(sv, dst, dst_pos);
}


inline std::uint64_t pack_characters_lenient(std::string_view sv,
                                             packed_word_vector& dst,
                                             std::uint64_t dst_pos = 0) {
  return pack_characters_<dna_alphabet::dna4, true>(sv, dst, dst_pos);
}


inline void unpack_characters(packed_word_vector const& src,
                              std::uint64_t length, std::string& dst) {
  return unpack_characters_<dna_alphabet::dna4>(src, length, dst);
}


// DNA16
inline std::uint64_t pack_characters_dna16(std::string_view sv,
                                           packed_word_vector& dst,
                                           std::uint64_t dst_pos = 0) {
  return pack_characters_<dna_alphabet::dna16>(sv, dst, dst_pos);
}

void print_read(std::span<std::uint64_t const> span, std::uint64_t len);
}  // namespace sf
