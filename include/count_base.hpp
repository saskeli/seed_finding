#pragma once

#include <cstdint>
#include <libbio/algorithm.hh>

#include "packed_read.hpp"


namespace sf {

template <typename t_gapmer>
struct count_base {
  typedef t_gapmer gapmer_type;
  constexpr static bool middle_gap_only{gapmer_type::middle_gap_only};

  // gapmer uses std::uint16_t for max_gap but we try to use a smaller type
  // here.
  static_assert(gapmer_type::max_gap <= UINT8_MAX);
  constexpr static uint8_t max_gap{gapmer_type::max_gap};

  using increment_fn_type = void (count_base::*)(gapmer_type);
  using increment_gapped_fn_type = void (count_base::*)(gapmer_type, uint64_t);

  virtual ~count_base() {}

  virtual void increment_signal_count(gapmer_type gg) = 0;
  virtual void increment_background_count(gapmer_type gg) = 0;

  virtual uint64_t offset(uint8_t gap_start, uint8_t gap_length) const = 0;

  virtual void increment_signal_count_gapped(gapmer_type gg, uint64_t off) = 0;
  virtual void increment_background_count_gapped(gapmer_type gg, uint64_t off) = 0;

  // Member function pointers for the template version of count_mers.
  struct increment_signal_counts {
    constexpr static inline increment_fn_type increment_fn{
        &count_base::increment_signal_count};
    constexpr static inline increment_gapped_fn_type increment_gapped_fn{
        &count_base::increment_signal_count_gapped};
  };

  struct increment_background_counts {
    constexpr static inline increment_fn_type increment_fn{
        &count_base::increment_background_count};
    constexpr static inline increment_gapped_fn_type increment_gapped_fn{
        &count_base::increment_background_count_gapped};
  };

  // Increment using the pointer-to-member functions in t_increment_count.
  // (Since the function pointers are constexpr (see above), at least LLVM 21
  // should completely optimise this away in subclasses marked final.)
  template <typename t_increment_count>
  void count_mers(packed_read_vector const& reads, uint8_t kk);

  void count_mers(packed_read_vector const& signal_reads,
                  packed_read_vector const& background_reads, uint8_t kk);
};


template <typename t_gapmer>
template <typename t_increment_count>
void count_base<t_gapmer>::count_mers(packed_read_vector const& reads,
                                      uint8_t kk) {
  // Helpers for calling the member functions via the pointers.
  // Not repeating the parameter lists in case the types need to be changed.
  auto const increment_count([this](auto&&... args) {
    (this->*t_increment_count::increment_fn)(
        std::forward<decltype(args)>(args)...);
  });

  auto const increment_gapped_count([this](auto&&... args) {
    (this->*t_increment_count::increment_gapped_fn)(
        std::forward<decltype(args)>(args)...);
  });

  // OpenMP’s parallel for should work with the range b.c. it is essentially a
  // memory range and the iterator is a pointer.
#pragma omp parallel for
  for (auto const& read : reads) {
    if (read.length < kk) continue;

    gapmer_type gg(read.packed_characters.front() >> (64U - 2U * kk), kk);
    increment_count(gg);
    read.iterate_characters(kk, [&](std::uint8_t const cc) {
      gg = gg.next_(cc);
      increment_count(gg);
    });

    uint8_t gap_s = middle_gap_only ? kk / 2 : 1;
    auto gap_lim =
        libbio::min_ct(read.length, middle_gap_only ? (kk + 3) / 2 : kk);

    for (; gap_s < gap_lim; ++gap_s) {
      for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
        auto const off{offset(gap_s, gap_l)};
        gg = gapmer_type(read.packed_characters, kk, gap_s, gap_l);
        increment_gapped_count(gg, off);
        read.iterate_character_pairs(
            gap_s, kk + gap_l,
            [&](std::uint8_t const lhsc, std::uint8_t const rhsc) {
              gg = gg.next_(lhsc, rhsc);
              increment_gapped_count(gg, off);
            });
      }
    }
  }  // for (auto const &read : reads)
}


template <typename t_gapmer>
void count_base<t_gapmer>::count_mers(
    packed_read_vector const& signal_reads,
    packed_read_vector const& background_reads, uint8_t kk) {
  count_mers<increment_signal_counts>(signal_reads, kk);
  count_mers<increment_background_counts>(background_reads, kk);
}
}  // namespace sf
