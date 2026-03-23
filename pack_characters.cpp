/*
 * Copyright (c) 2025-2026 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#include "pack_characters.hpp"

#include <bit>
#include <cstdint>
#include <iostream>
#include <libbio/assert.hh>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>


namespace sf {

void print_read(std::span<std::uint64_t const> span, std::uint64_t len) {
  std::cerr << "** read: ";
  for (auto word : span) {
    for (std::uint8_t ii{}; ii < 32U; ++ii) {
      word = std::rotl(word, 2);
      auto const cc{word & 0x3};
      switch (cc) {
        case 0:
          std::cerr << 'A';
          break;
        case 1:
          std::cerr << 'C';
          break;
        case 2:
          std::cerr << 'G';
          break;
        case 3:
          std::cerr << 'T';
          break;
        default:
          throw std::runtime_error("Unexpected character");
      }

      --len;
      if (0 == len) break;
    }
  }
  std::cerr << '\n';
}
}  // namespace sf


namespace {

template <sf::dna_alphabet t_alphabet, bool t_is_lenient>
struct pack_characters_trait {};


template <bool t_is_lenient>
struct pack_characters_trait<sf::dna_alphabet::dna4, t_is_lenient> {
  constexpr static std::uint8_t character_bits{0x2};
  constexpr static std::uint8_t character_max{0x3};

  static void push_unpacked(std::uint8_t const cc, std::string& dst) {
    libbio_assert_lt(cc, sf::dna4_characters.size());
    dst.push_back(sf::dna4_characters[cc]);
  }

  static bool push_packed(char const cc, std::uint64_t& dst_word,
                          std::uint64_t const dst_pos) {
    enum : uint64_t { A_ = 0x0, C_ = 0x1, G_ = 0x2, T_ = 0x3 };

    auto const chararcter_pos{dst_pos % 32U};
    std::uint64_t const shift_amt{62U - 2 * chararcter_pos};

    auto const append([&](uint64_t val) { dst_word |= val << shift_amt; });

    switch (cc) {
      case 'A':
      case 'a':
        return true;

      case 'C':
      case 'c':
        append(C_);
        return true;

      case 'G':
      case 'g':
        append(G_);
        return true;

      case 'T':
      case 't':
        append(T_);
        return true;

      [[unlikely]] default:
        if constexpr (t_is_lenient)
          return true;
        else
          return false;
    }
  }
};


template <>
struct pack_characters_trait<sf::dna_alphabet::dna16, false> {
  constexpr static std::uint8_t character_bits{0x4};
  constexpr static std::uint8_t character_max{0xf};

  static void push_unpacked(std::uint8_t const cc, std::string& dst) {
    libbio_assert_lt(cc, sf::dna16_characters.size());
    dst.push_back(sf::dna16_characters[cc]);
  }

  static bool push_packed(char const cc, std::uint64_t& dst_word,
                          std::uint64_t const dst_pos) {
    enum : uint64_t { A_ = 0x1, C_ = 0x2, G_ = 0x4, T_ = 0x8 };

    auto const chararcter_pos{dst_pos % 16U};
    std::uint64_t const shift_amt{60U - 4 * chararcter_pos};

    auto const append([&](uint64_t val) { dst_word |= val << shift_amt; });

    switch (cc) {
      case '.':
      case '-':
      case 'n': // Gap
        return true;

      case 'A':
        append(A_);
        return true;

      case 'C':
        append(C_);
        return true;

      case 'G':
        append(G_);
        return true;

      case 'T':
        append(T_);
        return true;

      case 'R':
        append(A_ | G_);
        return true;

      case 'Y':
        append(C_ | T_);
        return true;

      case 'S':
        append(G_ | C_);
        return true;

      case 'W':
        append(A_ | T_);
        return true;

      case 'K':
        append(G_ | T_);
        return true;

      case 'M':
        append(A_ | C_);
        return true;

      case 'B':
        append(C_ | G_ | T_);
        return true;

      case 'D':
        append(A_ | G_ | T_);
        return true;

      case 'H':
        append(A_ | C_ | T_);
        return true;

      case 'V':
        append(A_ | C_ | G_);
        return true;

      case 'N':
        append(A_ | C_ | G_ | T_);
        return true;

      [[unlikely]] default:
        return false;
    }
  }
};
}  // namespace


namespace sf {
template <dna_alphabet t_alphabet, bool t_is_lenient>
std::uint64_t pack_characters_(std::string_view sv, sf::packed_word_vector& dst,
                               std::uint64_t dst_pos) {
  pack_characters_trait<t_alphabet, t_is_lenient> trait;
  constexpr auto const character_count{UINT64_C(64) / trait.character_bits};

  auto it{sv.cbegin()};
  auto const end{sv.cend()};

  // If the current position is not divisible by character_count, dst
  // must be non-empty.
  if (dst_pos % character_count) {
    auto& dst_word{dst.back()};
    do {
      if (it == end) return dst_pos;

      if (!trait.push_packed(*it, dst_word, dst_pos)) [[unlikely]]
        return UINT64_MAX;

      ++dst_pos;
      ++it;
    } while (dst_pos % character_count);
  }

  if (it == end) return dst_pos;

  // 0 == dst_pos % character_count.
  while (true) {
    auto& dst_word{dst.emplace_back(0)};
    do {
      if (!trait.push_packed(*it, dst_word, dst_pos)) [[unlikely]]
        return UINT64_MAX;

      ++dst_pos;
      ++it;

      if (it == end) return dst_pos;
    } while (dst_pos % character_count);

    if (it == end) return dst_pos;
  }
}


template <dna_alphabet t_alphabet>
void unpack_characters_(packed_word_span const& src, std::uint64_t length,
                        std::string& dst) {
  if (src.empty()) return;

  pack_characters_trait<t_alphabet, false> trait;
  constexpr auto const character_count{UINT64_C(64) / trait.character_bits};

  std::uint64_t ii{};
  while (ii + character_count <= length) {
    auto word{src[ii / character_count]};
    for (std::uint8_t jj{}; jj < character_count; ++jj) {
      word = std::rotl(word, trait.character_bits);
      trait.push_unpacked(word & trait.character_max, dst);
    }

    ii += character_count;
  }

  auto word{src[ii / character_count]};
  while (ii < length) {
    word = std::rotl(word, trait.character_bits);
    trait.push_unpacked(word & trait.character_max, dst);
    ++ii;
  }
}


// Explicit instantiations.
template std::uint64_t pack_characters_<dna_alphabet::dna4>(
    std::string_view sv, packed_word_vector& dst, std::uint64_t dst_pos);

template std::uint64_t pack_characters_<dna_alphabet::dna4, true>(
    std::string_view sv, packed_word_vector& dst, std::uint64_t dst_pos);

template std::uint64_t pack_characters_<dna_alphabet::dna16>(
    std::string_view sv, packed_word_vector& dst, std::uint64_t dst_pos);

template void unpack_characters_<dna_alphabet::dna4>(
    packed_word_span const& src, std::uint64_t length, std::string& dst);

template void unpack_characters_<dna_alphabet::dna16>(
    packed_word_span const& src, std::uint64_t length, std::string& dst);

}  // namespace sf
