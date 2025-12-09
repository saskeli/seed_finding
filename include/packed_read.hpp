#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "pack_characters.hpp"
#include "packed_character_iteration.hpp"

namespace sf {
struct packed_read {
  packed_word_vector packed_characters;
  std::uint64_t length{};

  void unpack(std::string& dst) const {
    unpack_characters(packed_characters, length, dst);
  }

  void print() const { print_read(packed_characters, length); }

  template <typename t_cb>
  void iterate_characters(std::uint64_t start_pos, t_cb&& cb) const {
    iterate_packed_characters(packed_characters, length, start_pos, cb);
  }

  template <typename t_cb>
  void iterate_character_pairs(std::uint64_t lhs_start, std::uint64_t rhs_start,
                               t_cb&& cb) const {
    iterate_packed_character_pairs(packed_characters, length, lhs_start,
                                   rhs_start, cb);
  }
};

typedef std::vector<packed_read> packed_read_vector;
}  // namespace sf
