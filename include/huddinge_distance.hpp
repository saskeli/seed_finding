/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <hwy/aligned_allocator.h>
#include <hwy/highway.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <libbio/assert.hh>
#include <libbio/hwy_apply.hh>
#include <libbio/utility.hh>
#include <ostream>
#include <span>
#include <tuple>
#include <type_traits>


namespace sf::detail {

struct huddinge_distance_return_value {
  uint16_t distance{UINT16_MAX};
  int16_t position{};

  auto to_tuple() const { return std::make_tuple(distance, position); }
  bool operator==(huddinge_distance_return_value other) const {
    return to_tuple() == other.to_tuple();
  }
};


inline std::ostream& operator<<(std::ostream& os,
                                huddinge_distance_return_value rv) {
  os << "distance: " << rv.distance << " position: " << rv.position;
  return os;
}


template <typename t_tag>
class huddinge_distance_calculator {
 public:
  typedef libbio::hwy_apply<t_tag> hwy_apply_type;
  typedef typename hwy_apply_type::value_type value_type;
  typedef typename hwy_apply_type::vector_type vector_type;
  // The value type needs to be the same as the value type used for the
  // SIMD operations to make parallel storing possible.
  typedef std::array<value_type, 32> score_array;
  typedef std::array<value_type, hwy_apply_type::lanes> fast_path_score_array;

  template <bool t_uses_fast_path>
  using score_array_t =
      std::conditional_t<t_uses_fast_path, fast_path_score_array, score_array>;

  template <std::uint8_t t_length>
  using span_t = std::span<value_type const, t_length>;

 private:
  // Automatic formatting makes the following arrays less readable, so we turn
  // it off.
  // clang-format off
  alignas(HWY_ALIGNMENT) constexpr static value_type s_shift_left_amounts[]{
    0U,  2U,  4U,  6U,
    8U,  10U, 12U, 14U,
    16U, 18U, 20U, 22U,
    24U, 26U, 28U, 30U,
    32U, 34U, 36U, 38U,
    40U, 42U, 44U, 46U,
    48U, 50U, 52U, 54U,
    56U, 58U, 60U, 62U
  };

  // NEON does not seem to have rotates, at least implemented in Highway,
  // so we have the right shift amounts here.
  // Since counts not in [0, sizeof(T) * 8) yield implementation-defined results,
  // we shift in two phases to zero the first value.
  alignas(HWY_ALIGNMENT) constexpr static value_type s_shift_right_amounts[]{
    62U, 60U, 58U, 56U,
    54U, 52U, 50U, 48U,
    46U, 44U, 42U, 40U,
    38U, 36U, 34U, 32U,
    30U, 28U, 26U, 24U,
    22U, 20U, 18U, 16U,
    14U, 12U, 10U, 8U,
    6U,  4U,  2U,  0U
  };
  // clang-format on

  constexpr static value_type s_odd_mask{UINT64_C(0xAAAA'AAAA'AAAA'AAAA)};

  static_assert(32 == libbio::array_size(s_shift_left_amounts));
  static_assert(32 == libbio::array_size(s_shift_right_amounts));
  static_assert(32 == score_array{}.size());

  HWY_LANES_CONSTEXPR static hwy_apply_type s_apply{};

  huddinge_distance_return_value m_best_alignment{};
  value_type m_max_defined_characters{};

 private:
  template <typename t_cb>
  inline vector_type combine_diff_bits(t_cb& cb, vector_type const diff) const;

  template <typename t_cb>
  inline void update_scores_from_diff(t_cb& cb, vector_type const combined_diff,
                                      vector_type const combined_masks,
                                      score_array& scores) const;

  template <bool t_forward_direction>
  inline void update_best_alignment(score_array& scores, int16_t base_position);

  template <bool t_forward_direction, typename t_cb>
  inline void update_best_alignment_fast_path(
      t_cb& cb, vector_type const combined_diff,
      vector_type const combined_masks,
      vector_type const max_defined_characters, fast_path_score_array& scores);

  template <bool t_forward_direction, bool t_uses_fast_path>
  inline void shift_word_and_compare_(score_array_t<t_uses_fast_path>& scores,
                                      uint64_t lhs, uint64_t rhs,
                                      uint64_t lhs_mask, uint64_t rhs_mask,
                                      uint16_t lhsl);

 public:
  huddinge_distance_return_value best_alignment() const {
    return m_best_alignment;
  }

  template <std::size_t t_n>
  inline uint16_t count_defined_characters(span_t<t_n> mask) const;

  template <std::size_t t_lhsn, std::size_t t_rhsn>
  inline void determine_max_defined_characters(span_t<t_lhsn> lhs_mask,
                                               span_t<t_rhsn> rhs_mask);

  template <bool t_forward_direction, std::size_t t_lhsn, std::size_t t_rhsn>
  inline void shift_word_pair_and_compare(
      span_t<t_lhsn> lhs, span_t<t_rhsn> rhs, span_t<t_lhsn> lhs_mask,
      span_t<t_rhsn> rhs_mask, std::size_t const lhs_start_idx, uint16_t lhsl);

  template <bool t_forward_direction, std::size_t t_lhsn>
  inline void shift_word_and_compare(uint64_t lhs, uint64_t rhs,
                                     uint64_t lhs_mask, uint64_t rhs_mask,
                                     uint16_t lhsl);
};


template <typename t_tag>
template <std::size_t t_n>
uint16_t huddinge_distance_calculator<t_tag>::count_defined_characters(
    span_t<t_n> mask) const {
  // Likely no need to use vector instructions since t_n ≤ 2.
  uint16_t retval{};
  for (auto const word : mask) retval += std::popcount(word);
  return retval;
}


template <typename t_tag>
template <std::size_t t_lhsn, std::size_t t_rhsn>
void huddinge_distance_calculator<t_tag>::determine_max_defined_characters(
    span_t<t_lhsn> lhs_mask, span_t<t_rhsn> rhs_mask) {
  m_max_defined_characters = std::max(count_defined_characters(lhs_mask),
                                      count_defined_characters(rhs_mask));
}


template <typename t_tag>
template <typename t_cb>
auto huddinge_distance_calculator<t_tag>::combine_diff_bits(
    t_cb& cb, vector_type const diff) const -> vector_type {
  namespace hn = hwy::HWY_NAMESPACE;

  // Mask the odd bits, shift right and OR with the even bits.
  // Consequently the even bits in the return value will be of interest.
  auto const odd_mask{cb.set(s_odd_mask)};
  auto const odd_diff{diff & odd_mask};
  auto const odd_diff_{hn::ShiftRight<1>(odd_diff)};
  return diff | odd_diff_;
}


template <typename t_tag>
template <typename t_cb>
void huddinge_distance_calculator<t_tag>::update_scores_from_diff(
    t_cb& cb, vector_type const combined_diff, vector_type const combined_masks,
    score_array& scores) const {
  namespace hn = hwy::HWY_NAMESPACE;

  // Flip the bits in combined_diff to get ones for even positions of the
  // matching characters, mask and count. We assumed that only the even bits are
  // set in the mask.
  auto const matching_positions{hn::AndNot(combined_diff, combined_masks)};
  auto const matching_counts{hn::PopulationCount(matching_positions)};

  // Add to the scores.
  auto const old_scores{cb.load(scores.data())};
  auto const new_scores{old_scores + matching_counts};
  cb.store(new_scores, scores.data());
}


template <typename t_tag>
template <bool t_forward_direction, typename t_cb>
void huddinge_distance_calculator<t_tag>::update_best_alignment_fast_path(
    t_cb& cb, vector_type const combined_diff, vector_type const combined_masks,
    vector_type const max_defined_characters,
    fast_path_score_array& distances) {
  // Do the same as in update_scores_from_diff() but since we need not calculate
  // the score in multiple iterations, we can determine the distance
  // immediately.

  namespace hn = hwy::HWY_NAMESPACE;

  auto const matching_positions{hn::AndNot(combined_diff, combined_masks)};
  auto const matching_counts{hn::PopulationCount(matching_positions)};
  auto const distances_{max_defined_characters - matching_counts};
  cb.store_(distances_, distances.data());

  // We check the limit from cb.count() b.c. we may be currently handling less
  // values than there are SIMD lanes.
  for (std::size_t jj{}; jj < cb.count(); ++jj) {
    if (distances[jj] < m_best_alignment.distance) {
      m_best_alignment.distance = distances[jj];
      m_best_alignment.position = cb.ii + jj;

      if constexpr (t_forward_direction)
        m_best_alignment.position = -m_best_alignment.position;
    }
  }
}


template <typename t_tag>
template <bool t_forward_direction>
void huddinge_distance_calculator<t_tag>::update_best_alignment(
    score_array& scores, int16_t base_position) {
  // Calculate the distances from the scores and find the best score.
  std::conditional_t<t_forward_direction, std::plus<>, std::minus<>> op{};

  auto const max_defined_characters_{s_apply.set(m_max_defined_characters)};
  alignas(HWY_ALIGNMENT) std::array<uint64_t, s_apply.lanes> distances{};
  s_apply(32U, std::false_type{}, [&](auto& cb) {
    auto const scores_{cb.load(scores.data())};
    auto const distances_{max_defined_characters_ - scores_};
    cb.store_(distances_, distances.data());

    for (std::size_t jj{}; jj < distances.size(); ++jj) {
      if (distances[jj] < m_best_alignment.distance) {
        m_best_alignment.distance = distances[jj];
        m_best_alignment.position = op(base_position, cb.ii + jj);
      }
    }
  });
}


template <typename t_tag>
template <bool t_forward_direction, std::size_t t_lhsn, std::size_t t_rhsn>
void huddinge_distance_calculator<t_tag>::shift_word_pair_and_compare(
    span_t<t_lhsn> lhs, span_t<t_rhsn> rhs, span_t<t_lhsn> lhs_mask,
    span_t<t_rhsn> rhs_mask, std::size_t const lhs_start_idx, uint16_t lhsl) {
  namespace hn = hwy::HWY_NAMESPACE;

  alignas(HWY_ALIGNMENT) score_array scores{};

  // Compare the characters starting from lhs_start_idx in such a way that it is
  // followed by a (possibly non-full) word. Since the current word always has
  // 32 characters (and the next one possibly less), we can do 32 shifts without
  // checking. We assume that HWY_HAVE_CONSTEXPR_LANES is defined to one.
  static_assert(0 == 32 % hwy_apply_type::lanes);
  constexpr auto needs_to_consider_remaining_values{std::false_type{}};

  std::size_t lhs_word_idx{lhs_start_idx};
  std::size_t rhs_word_idx{};
  while (lhs_word_idx < t_lhsn - 1 && rhs_word_idx < t_rhsn) {
    auto const lhs_word{s_apply.set(lhs[lhs_word_idx])};
    auto const lhs_next_word{s_apply.set(lhs[lhs_word_idx + 1])};
    auto const rhs_word{s_apply.set(rhs[rhs_word_idx])};
    auto const lhs_word_mask{s_apply.set(lhs_mask[lhs_word_idx])};
    auto const lhs_next_word_mask{s_apply.set(lhs_mask[lhs_word_idx + 1])};
    auto const rhs_word_mask{s_apply.set(rhs_mask[rhs_word_idx])};

    s_apply(32U, needs_to_consider_remaining_values, [&](auto& cb) {
      // Shift to the left since we have one word remaining.
      auto const lhs_shift_amounts{cb.load(s_shift_left_amounts)};
      auto const lhs_shifted{lhs_word << lhs_shift_amounts};

      // Same with the next word. We have one due to lhs_word_idx_ < t_lhsn - 1.
      auto const lhs_next_shift_amounts{cb.load(s_shift_right_amounts)};
      auto const lhs_next_shifted{hn::ShiftRight<2>(lhs_next_word) >>
                                  lhs_next_shift_amounts};

      // Combine.
      auto const lhs_combined{lhs_shifted | lhs_next_shifted};

      // Compare by using XOR and get a single bit for each character indicating
      // whether there is a mismatch.
      auto const diff{lhs_combined ^ rhs_word};
      auto const combined_diff{combine_diff_bits(cb, diff)};

      // Shift the masks and combine.
      auto const lhs_word_mask_shifted{lhs_word_mask << lhs_shift_amounts};
      auto const lhs_next_word_mask_shifted{
          hn::ShiftRight<2>(lhs_next_word_mask) >> lhs_next_shift_amounts};
      auto const combined_masks{
          (lhs_word_mask_shifted | lhs_next_word_mask_shifted) & rhs_word_mask};

      update_scores_from_diff(cb, combined_diff, combined_masks, scores);
    });

    ++lhs_word_idx;
    ++rhs_word_idx;
  }

  if (rhs_word_idx != t_rhsn)
    shift_word_and_compare_<t_forward_direction, false>(
        scores, lhs[lhs_word_idx], rhs[rhs_word_idx], lhs_mask[lhs_word_idx],
        rhs_mask[rhs_word_idx], lhsl);

  update_best_alignment<t_forward_direction>(
      scores, (t_forward_direction ? -1 : 1) * 32 * int16_t(lhs_start_idx));
}


template <typename t_tag>
template <bool t_forward_direction, bool t_uses_fast_path>
void huddinge_distance_calculator<t_tag>::shift_word_and_compare_(
    score_array_t<t_uses_fast_path>& scores, uint64_t lhs, uint64_t rhs,
    uint64_t lhs_mask, uint64_t rhs_mask, uint16_t lhsl) {
  auto const lhs_{s_apply.set(lhs)};
  auto const rhs_{s_apply.set(rhs)};
  auto const lhs_mask_{s_apply.set(lhs_mask)};
  auto const rhs_mask_{s_apply.set(rhs_mask)};
  auto const max_defined_characters([&] {
    struct empty {};

    if constexpr (t_uses_fast_path)
      return s_apply.set(m_max_defined_characters);
    else
      return empty{};
  }());

  s_apply(uint16_t((lhsl % 32) ?: 32), [&](auto& cb) {
    // Shift to the left since we have one word remaining.
    auto const lhs_shift_amounts{cb.load(s_shift_left_amounts)};
    auto const lhs_shifted{lhs_ << lhs_shift_amounts};

    // Compare by using XOR and get a single bit for each character indicating
    // whether there is a mismatch.
    auto const diff{lhs_shifted ^ rhs_};
    auto const combined_diff{combine_diff_bits(cb, diff)};

    // Shift the masks and combine.
    auto const lhs_word_mask_shifted{lhs_mask_ << lhs_shift_amounts};
    auto const combined_masks{lhs_word_mask_shifted & rhs_mask_};

    if constexpr (t_uses_fast_path)
      update_best_alignment_fast_path<t_forward_direction>(
          cb, combined_diff, combined_masks, max_defined_characters, scores);
    else
      update_scores_from_diff(cb, combined_diff, combined_masks, scores);
  });
}


template <typename t_tag>
template <bool t_forward_direction, std::size_t t_lhsn>
void huddinge_distance_calculator<t_tag>::shift_word_and_compare(
    uint64_t lhs, uint64_t rhs, uint64_t lhs_mask, uint64_t rhs_mask,
    uint16_t lhsl) {
  constexpr bool uses_fast_path{1 == t_lhsn};
  alignas(HWY_ALIGNMENT) score_array_t<uses_fast_path> scores{};

  shift_word_and_compare_<t_forward_direction, uses_fast_path>(
      scores, lhs, rhs, lhs_mask, rhs_mask, lhsl);

  if constexpr (not uses_fast_path)
    update_best_alignment<t_forward_direction>(
        scores, (t_forward_direction ? -1 : 1) * 32 * int16_t(t_lhsn - 1));
}
}  // namespace sf::detail


namespace sf {

typedef detail::huddinge_distance_return_value huddinge_distance_return_value;

// Calculate the Huddinge distance from the spans that contain packed
// characters. Gaps are indicated in the mask spans in such a way that only the
// even bit of each non-gap character is set. We assume that the spans are
// fairly short (at most two words in practice) so we need not have a threshold
// value for t_lhsn and t_rhsn to produce generic code.
template <std::size_t t_lhsn, std::size_t t_rhsn>
HWY_ATTR huddinge_distance_return_value huddinge_distance(
    std::span<uint64_t const, t_lhsn> lhs,
    std::span<uint64_t const, t_rhsn> rhs,
    std::span<uint64_t const, t_lhsn> lhs_mask,
    std::span<uint64_t const, t_rhsn> rhs_mask, uint16_t lhsl, uint16_t rhsl) {
  static_assert(0 < t_lhsn);
  static_assert(0 < t_rhsn);

  // Sanity check; we expect the span sizes not to be larger than needed.
  libbio_assert_eq(t_lhsn, (lhsl + 31) / 32);
  libbio_assert_eq(t_rhsn, (rhsl + 31) / 32);

  namespace hn = hwy::HWY_NAMESPACE;
  typedef hn::ScalableTag<std::uint64_t> tag_type;
  typedef detail::huddinge_distance_calculator<tag_type>
      distance_calculator_type;
  distance_calculator_type dc;

  dc.determine_max_defined_characters(lhs_mask, rhs_mask);

  // Case where lhs is shifted to the left, i.e. has negative or zero position
  // w.r.t. rhs.
  for (std::size_t lhs_word_idx{}; lhs_word_idx < t_lhsn - 1; ++lhs_word_idx)
    dc.shift_word_pair_and_compare<false>(lhs, rhs, lhs_mask, rhs_mask,
                                          lhs_word_idx, lhsl);

  dc.shift_word_and_compare<false, t_lhsn>(
      lhs.back(), rhs.front(), lhs_mask.back(), rhs_mask.front(), lhsl);

  // Case where lhs is shifted to the right, i.e. has positive or zero position
  // w.r.t. rhs. (To avoid complexity, we handle position zero twice.)
  for (std::size_t rhs_word_idx{}; rhs_word_idx < t_rhsn - 1; ++rhs_word_idx)
    dc.shift_word_pair_and_compare<true>(rhs, lhs, rhs_mask, lhs_mask,
                                         rhs_word_idx, rhsl);

  dc.shift_word_and_compare<true, t_rhsn>(
      rhs.back(), lhs.front(), rhs_mask.back(), lhs_mask.front(), rhsl);

  return dc.best_alignment();
}
}  // namespace sf
