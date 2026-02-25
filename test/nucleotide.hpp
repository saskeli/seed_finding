/*
 * Copyright (c) 2025-2026 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <array>
#include <cstdint>
#include <ostream>
#include <stdexcept>

#include "../include/pack_characters.hpp"
#include "libbio/assert.hh"
#include "test.hpp"


namespace {

template <sf::dna_alphabet t_alphabet>
struct dna_value_to_character_ {};

template <>
struct dna_value_to_character_<sf::dna_alphabet::dna4> {
  char operator()(uint8_t vv) const {
    libbio_assert_lt(vv, 4);
    return sf::dna4_characters[vv];
  }
};

template <>
struct dna_value_to_character_<sf::dna_alphabet::dna16> {
  char operator()(uint8_t vv) const {
    libbio_assert_lt(vv, 16);
    return sf::dna16_characters[vv];
  }
};

template <sf::dna_alphabet t_alphabet>
char dna_value_to_character(uint8_t vv) {
  return dna_value_to_character_<t_alphabet>{}(vv);
}


template <sf::dna_alphabet t_dna_alphabet>
struct nucleotide_tpl {
  constexpr static uint8_t alphabet_size{sf::alphabet_size_v<t_dna_alphabet>};
  char value{};

  /* implicit */ operator char() const { return value; }
};


template <sf::dna_alphabet t_dna_alphabet>
struct packed_nucleotide_tpl {
  constexpr static uint8_t alphabet_size{sf::alphabet_size_v<t_dna_alphabet>};
  std::uint8_t value{};

  /* implicit */ operator char() const {
    return dna_value_to_character<t_dna_alphabet>(value);
  }
};


typedef nucleotide_tpl<sf::dna_alphabet::dna4> nucleotide;
typedef nucleotide_tpl<sf::dna_alphabet::dna16> dna16_nucleotide;
typedef packed_nucleotide_tpl<sf::dna_alphabet::dna4> packed_nucleotide;
typedef packed_nucleotide_tpl<sf::dna_alphabet::dna16> packed_dna16_nucleotide;


inline std::ostream& operator<<(std::ostream& os, nucleotide const nt) {
  os << char(nt);
  return os;
}


inline std::ostream& operator<<(std::ostream& os, packed_nucleotide const nt) {
  os << char(nt);
  return os;
}


char complement_nt(char const cc) {
  switch (cc) {
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
      return 'A';
    case '.':
      return '.';
    default:
      throw std::runtime_error("Unexpected character");
  }
}
}  // namespace


namespace rc {

template <sf::dna_alphabet t_alphabet>
struct Arbitrary<packed_nucleotide_tpl<t_alphabet>> {
  static Gen<packed_nucleotide_tpl<t_alphabet>> arbitrary() {
    typedef packed_nucleotide_tpl<t_alphabet> return_type;
    return gen::construct<return_type>(gen::inRange(
        std::uint8_t{}, return_type::alphabet_size)); // Half-open range.
  }
};


template <sf::dna_alphabet t_alphabet>
struct Arbitrary<nucleotide_tpl<t_alphabet>> {
  static Gen<nucleotide_tpl<t_alphabet>> arbitrary() {
    typedef packed_nucleotide_tpl<t_alphabet> packed_type;
    typedef nucleotide_tpl<t_alphabet> return_type;
    return gen::map(gen::arbitrary<packed_type>(),
                    [](auto const pn) { return return_type(pn); });
  }
};
}  // namespace rc
