#pragma once

#include <stdlib.h>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <libbio/algorithm.hh>
#include <sdsl/bit_vectors.hpp>
#include <vector>

#include "count_base.hpp"
#include "packed_read.hpp"

namespace sf {

// t_base used only for having a type name for the base class.
template <typename t_gapmer, typename t_base = count_base <t_gapmer>>
class gapmer_count final : public t_base {
 public:
  typedef t_base base_type;
  typedef double value_type;
  const static constexpr uint64_t ONE = 1;

  typedef base_type::gapmer_type gapmer_type;
  using base_type::max_gap;
  using base_type::middle_gap_only;
  using base_type::count_mers;

  struct count_pair {
    value_type signal_count{};
    value_type background_count{};
  };

  static constexpr uint64_t lookup_elems(uint8_t k) {
    uint64_t ret = ONE << (k * 2);
    uint64_t add = ret * max_gap;
    uint64_t mul = k - 1;
    if constexpr (middle_gap_only) {
      mul = 1 + k % 2;
    }
    add *= mul;
    ret += add;
    return ret;
  }

  static constexpr uint64_t lookup_bytes(uint8_t k) {
    uint64_t ret = lookup_elems(k);
    ret += (ret + 63) / 64 + ret;
    ret += ONE << (k * 2 + 1);
    return ret * sizeof(uint64_t);
  }

  static constexpr uint64_t offset(uint8_t k, uint8_t gap_s, uint8_t gap_l) {
    uint64_t min_gap = middle_gap_only ? k / 2 : 1;
    uint64_t step = ONE << (k * 2);
    uint64_t ret = step + (gap_s - min_gap) * step * max_gap;
    ret += step * (gap_l - 1);
    ret *= gap_s > 0;
    return ret;
  }

 private:
  typedef std::vector<value_type> value_vector;

  value_vector sig_counts_;
  value_vector bg_counts_;
  sdsl::bit_vector discarded_;
  uint8_t k_{};

 private:
  void increment_signal_count(gapmer_type gg) override {
    uint64_t cv{gg.value()};
#pragma omp atomic relaxed
    ++sig_counts_[cv];
  }

  void increment_background_count(gapmer_type gg) override {
    uint64_t cv{gg.value()};
#pragma omp atomic relaxed
    ++bg_counts_[cv];
  }

  void increment_signal_count_gapped(gapmer_type gg, uint64_t off) override {
    uint64_t cv{off + gg.value()};
#pragma omp atomic relaxed
    ++sig_counts_[cv];
  }

  void increment_background_count_gapped(gapmer_type gg,
                                         uint64_t off) override {
    uint64_t cv{off + gg.value()};
#pragma omp atomic relaxed
    ++bg_counts_[cv];
  }

  private:
   gapmer_count(packed_read_vector const& sig_reads,
                packed_read_vector const& bg_reads, uint64_t size, uint8_t k)
       : sig_counts_(size), bg_counts_(size), discarded_(size), k_(k) {
     count_mers(sig_reads, bg_reads, k_);
   }

  public:
   gapmer_count(packed_read_vector const& sig_reads,
                packed_read_vector const& bg_reads, uint8_t k)
       : gapmer_count(sig_reads, bg_reads, lookup_elems(k), k) {}

   gapmer_count() = default;

 private:
  void smooth(value_type* arr, value_type* brr, value_type* sig_scratch, value_type* bg_scratch,
              uint64_t v_lim) {
#pragma omp parallel for
    for (uint64_t v = 0; v < v_lim; ++v) {
      sig_scratch[v] = arr[v];
      bg_scratch[v] = brr[v];
      std::array<uint64_t, 10> addables{};
      std::array<uint64_t, 10> bddables{};
      uint16_t smallest = 0;
      for (uint64_t x_v = 1; x_v < 4; ++x_v) {
        for (uint64_t shl = 0; shl < 2 * k_; shl += 2) {
          uint64_t o_v = v ^ (x_v << shl);
          uint64_t o_c = arr[o_v];
          if (addables[smallest] < o_c) {
            addables[smallest] = o_c;
            bddables[smallest] = brr[o_v];
            smallest = 0;
            for (uint16_t i = 1; i < 10; ++i) {
              smallest = addables[i] < addables[smallest] ? i : smallest;
            }
          }
        }
      }
      for (uint16_t i = 0; i < 10; ++i) {
        sig_scratch[v] += addables[i];
        bg_scratch[v] += bddables[i];
      }
      sig_scratch[v] /= 11;
      bg_scratch[v] /= 11;
    }
  }

 public:
  void smooth() {
    uint64_t v_lim = ONE << (k_ * 2);
    value_type* sig_scratch = (value_type*)calloc(v_lim, sizeof(value_type));
    value_type* bg_scratch = (value_type*)calloc(v_lim, sizeof(value_type));
    smooth(sig_counts_.data(), bg_counts_.data(), sig_scratch, bg_scratch, v_lim);
    std::memcpy(sig_counts_.data(), sig_scratch, v_lim * sizeof(value_type));
    std::memcpy(bg_counts_.data(), bg_scratch, v_lim * sizeof(value_type));
    uint8_t gap_s = middle_gap_only ? k_ / 2 : 1;
    uint8_t gap_lim = middle_gap_only ? (k_ + 3) / 2 : k_;
    for (; gap_s < gap_lim; ++gap_s) {
      for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
        uint64_t off = offset(gap_s, gap_l);
        value_type* l_count = sig_counts_.data() + off;
        value_type* r_count = bg_counts_.data() + off;
        std::memset(sig_scratch, 0, v_lim * sizeof(value_type));
        std::memset(bg_scratch, 0, v_lim * sizeof(value_type));
        smooth(l_count, r_count, sig_scratch, bg_scratch, v_lim);
        std::memcpy(l_count, sig_scratch, v_lim * sizeof(value_type));
        std::memcpy(r_count, bg_scratch, v_lim * sizeof(value_type));
      }
    }
    free(sig_scratch);
    free(bg_scratch);
  }

  count_pair count(gapmer_type gg) const {
    uint64_t off = offset(gg.gap_start(), gg.gap_length());
    return count(gg, off);
  }

  count_pair count(gapmer_type gg, uint64_t off) const {
    uint64_t cv{off + gg.value()};
    return {sig_counts_[cv], bg_counts_[cv]};
  }

  uint64_t offset(uint8_t gap_s, uint8_t gap_l) const override {
    return offset(k_, gap_s, gap_l);
  }

  bool is_discarded(gapmer_type gg, uint64_t off) const {
    return is_discarded_(gg.value(), off);
  }

  void mark_discarded(gapmer_type gg, uint64_t off) {
    mark_discarded_(gg.value(), off);
  }

  // gapmer_type is implicitly convertible to uint64_t, so
  // we refrain from overloading is_discarded() and mark_discarded().
  bool is_discarded_(uint64_t vv, uint64_t off) const {
    return discarded_[vv + off];
  }

  void mark_discarded_(uint64_t vv, uint64_t off) {
    discarded_[vv + off] = true;
  }
};
}  // namespace sf
