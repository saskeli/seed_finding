#pragma once

#include <cstdint>

#include "../include/gapmer.hpp"
#include "../include/gapmer_count.hpp"

namespace {
template <bool t_middle_gap_only, uint8_t t_max_gap>
using gapmer_count_t =
    sf::gapmer_count<sf::gapmer<t_middle_gap_only, t_max_gap>>;
}

namespace sf {
TEST(GapmerCount, Offset1) {
  typedef gapmer_count_t<true, 6> G_C;
  ASSERT_EQ(G_C::offset(5, 0, 0), 0ull);
}

TEST(GapmerCount, Offset2) {
  typedef gapmer_count_t<true, 6> G_C;
  uint64_t c_val = 1;
  c_val <<= 10;
  ASSERT_EQ(G_C::offset(5, 2, 1), c_val);
}

TEST(GapmerCount, Offset3) {
  typedef gapmer_count_t<true, 6> G_C;
  uint64_t p_width = 1;
  p_width <<= 10;
  uint64_t c_val = p_width + p_width * 6;
  c_val += p_width;
  ASSERT_EQ(G_C::offset(5, 3, 2), c_val);
}

TEST(GapmerCount, ElemCount) {
  typedef gapmer_count_t<true, 6> G_C;
  uint64_t p_width = 1;
  p_width <<= 10;
  uint64_t c_val = p_width + p_width * 6 * 2;
  ASSERT_EQ(G_C::lookup_elems(5), c_val);
}

TEST(GapmerCount, Offset4) {
  typedef gapmer_count_t<false, 6> G_C;
  uint64_t c_val = 1;
  c_val <<= 10;
  ASSERT_EQ(G_C::offset(5, 1, 1), c_val);
}

TEST(GapmerCount, Offset5) {
  typedef gapmer_count_t<false, 6> G_C;
  uint64_t p_width = 1;
  p_width <<= 10;
  uint64_t c_val = p_width + p_width * 6 * 2;
  c_val += p_width;
  ASSERT_EQ(G_C::offset(5, 3, 2), c_val);
}

TEST(GapmerCount, ElemCount2) {
  typedef gapmer_count_t<false, 6> G_C;
  uint64_t p_width = 1;
  p_width <<= 10;
  uint64_t c_val = p_width + p_width * 6 * 4;
  ASSERT_EQ(G_C::lookup_elems(5), c_val);
}

}  // namespace sf
