/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <cstdint>
#include <ostream>
#include <stdexcept>

#include "test.hpp"


namespace {

struct nucleotide {
  char value{};

  /* implicit */ operator char() const { return value; }
};


struct packed_nucleotide {
  std::uint8_t value{};

  /* implicit */ inline operator char() const;
};


packed_nucleotide::operator char() const {
  switch (value) {
    case 0x0:
      return 'A';
    case 0x1:
      return 'C';
    case 0x2:
      return 'G';
    case 0x3:
      return 'T';
    default:
      throw std::runtime_error{"Unexpected packed value"};
  }
}


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

template <>
struct Arbitrary<packed_nucleotide> {
  static Gen<packed_nucleotide> arbitrary() {
    return gen::construct<packed_nucleotide>(
        gen::inRange(std::uint8_t{}, std::uint8_t{0x4})); // Half-open range.
  }
};


template <>
struct Arbitrary<nucleotide> {
  static Gen<nucleotide> arbitrary() {
    return gen::map(gen::arbitrary<packed_nucleotide>(),
                    [](auto const pn) { return nucleotide(pn); });
  }
};
}  // namespace rc
