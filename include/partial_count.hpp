#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <libbio/algorithm.hh>

#include "count_base.hpp"
#include "lossy_hash_map.hpp"

namespace sf {

// t_base used only for having a type name for the base class.
template <class t_gapmer, typename t_base = count_base<t_gapmer>>
class partial_count final : public t_base {
 public:
  typedef t_base base_type;
  typedef base_type::gapmer_type gapmer_type;
  typedef base_type::counting_context counting_context;
  using base_type::count_mers;
  using base_type::max_gap;
  using base_type::middle_gap_only;

  typedef uint64_t value_type;

  struct count_pair {
    value_type signal_count{};
    value_type background_count{};
  };

  struct smooth_count_pair {
    double signal_count{};
    double background_count{};
  };

 private:
  typedef lossy_hash_map<gapmer_type, count_pair, typename gapmer_type::hash>
      count_map;
  count_map counts_;

 public:
  partial_count() : counts_(typename count_map::allocate_tag{}) {}

  void init(gapmer_type gg) { counts_.init(gg); }
  void clear() { counts_.clear(); }

  void increment_signal_count(gapmer_type gg,
                              counting_context const&) override {
    if (auto vp{counts_.find(gg)}; vp) {
#pragma omp atomic relaxed
      ++vp->signal_count;
    }
  }

  void increment_background_count(gapmer_type gg,
                                  counting_context const&) override {
    if (auto vp{counts_.find(gg)}; vp) {
#pragma omp atomic relaxed
      ++vp->background_count;
    }
  }

  uint64_t offset(uint8_t, uint8_t) const override { return 0; }

  void increment_signal_count_gapped(gapmer_type gg, uint64_t,
                                     counting_context const& ctx) override {
    increment_signal_count(gg, ctx);
  }

  void increment_background_count_gapped(gapmer_type gg, uint64_t,
                                         counting_context const& ctx) override {
    increment_background_count(gg, ctx);
  }

  count_pair count(gapmer_type gg) const {
    if (auto vp{counts_.find(gg)}; vp) {
      return {vp->signal_count, vp->background_count};
    }
    return {0, 0};
  }

  smooth_count_pair smooth_count(gapmer_type gg) const {
    auto const counts{count(gg)};
    double sig_val(counts.signal_count);
    double bg_val(counts.background_count);
    std::array<uint64_t, 10> addables{};
    std::array<uint64_t, 10> bddables{};
    uint16_t smallest = 0;

    gg.hamming_neighbours([&](gapmer_type oo) {
      auto const o_counts{count(oo)};
      if (addables[smallest] < o_counts.signal_count) {
        addables[smallest] = o_counts.signal_count;
        bddables[smallest] = o_counts.background_count;
        smallest = 0;
        for (uint16_t i = 1; i < 10; ++i) {
          smallest = addables[i] < addables[smallest] ? i : smallest;
        }
      }
    });

    for (uint16_t i = 0; i < 10; ++i) {
      sig_val += addables[i];
      bg_val += bddables[i];
    }
    return {sig_val / 11, bg_val / 11};
  }

  double fill_rate() const { return counts_.fill_rate(); }
};

}  // namespace sf
