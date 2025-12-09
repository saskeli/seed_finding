#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <libbio/algorithm.hh>
#include <utility>
#include "packed_read.hpp"

namespace sf {

template <class gapmer_t>
class partial_count {
 private:
  static const constexpr uint64_t MOD = 59999999;
  static const constexpr uint64_t STEP = 40960001;
  struct elem {
    gapmer_t mer;
    uint64_t sig_count;
    uint64_t bg_count;
  };

  elem* data_;
  uint64_t size_;

  uint64_t find(gapmer_t g) const {
    uint64_t g_val(g);
    uint64_t idx = g_val % MOD;
    uint64_t val(data_[idx].mer);
    while (val > 0 && val != g_val) [[unlikely]] {
      idx += STEP;
      idx -= idx >= MOD ? MOD : 0;
      val = data_[idx].mer;
    }
    return idx;
  }

 public:
  partial_count(): size_(0) {
    data_ = (elem*)calloc(MOD, sizeof(elem));
  }

  ~partial_count() {
    if (data_ != nullptr) {
      free(data_);
    }
  }

  void init(gapmer_t g) {
    uint64_t idx = find(g);
    if (data_[idx].mer != g) {
      data_[idx] = {g, 0, 0};
      ++size_;
    }
  }

  void clear() {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
    std::memset(data_, 0, MOD * sizeof(elem));
#pragma GCC diagnostic pop
    size_ = 0;
  }

  void inc_sig(gapmer_t g) {
    uint64_t idx = find(g);
    if (data_[idx].mer == g) {
#pragma omp atomic
      data_[idx].sig_count = data_[idx].sig_count + 1;
    }
  }

  void inc_bg(gapmer_t g) {
    uint64_t idx = find(g);
    if (data_[idx].mer == g) {
#pragma omp atomic
      data_[idx].bg_count = data_[idx].bg_count + 1;
    }
  }

  template <bool middle_gap_only, uint16_t max_gap>
  void count_mers(packed_read_vector const &sig_reads, packed_read_vector const &bg_reads,
                  uint16_t k) {
    auto const do_count([&](packed_read_vector const &reads, bool should_check_k, auto&& inc) {
      // OpenMP’s parallel for should work with the range b.c. it is essentially a
      // memory range and the iterator is a pointer.
#pragma omp parallel for
      for (auto const &read : reads) {
        if (should_check_k && read.length < k) continue;

        gapmer_t g(read.packed_characters.front() >> (64U - 2U * k), k);
        inc(g);
        read.iterate_characters(k, [&](std::uint8_t const cc) {
          g = g.next_(cc);
          inc(g);
        });

        uint8_t gap_s = middle_gap_only ? k / 2 : 1;
        auto gap_lim = libbio::min_ct(read.length,
                                      middle_gap_only ? (k + 3) / 2 : k);

        for (; gap_s < gap_lim; ++gap_s) {
          for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
            g = gapmer_t(read.packed_characters.front(), k, gap_s, gap_l);
            inc(g);
            read.iterate_character_pairs(
                gap_s, k + gap_l,
                [&](std::uint8_t const lhsc, std::uint8_t const rhsc) {
                  g = g.next_(lhsc, rhsc);
                  inc(g);
                });
          }
        }
      } // for (auto const &read : reads)
    });

    do_count(sig_reads, true, [this](auto const gg) { inc_sig(gg); });
    do_count(bg_reads, false, [this](auto const gg) { inc_bg(gg); });
  }

  std::pair<uint64_t, uint64_t> count(gapmer_t g) const {
    uint64_t idx = find(g);
    if (data_[idx].mer == g) {
      return {data_[idx].sig_count, data_[idx].bg_count};
    }
    return {0, 0};
  }

  std::pair<double, double> smooth_count(gapmer_t g) const {
    uint64_t idx = find(g);
    double sig_val = data_[idx].sig_count;
    double bg_val = data_[idx].bg_count;
    std::array<uint64_t, 10> addables{};
    std::array<uint64_t, 10> bddables{};
    uint16_t smallest = 0;
    auto callback = [&](gapmer_t o) {
      uint64_t o_idx = find(o);
      uint64_t o_sig = data_[o_idx].sig_count;
      if (addables[smallest] < o_sig) {
        addables[smallest] = o_sig;
        bddables[smallest] = data_[o_idx].bg_count;
        smallest = 0;
        for (uint16_t i = 1; i < 10; ++i) {
          smallest = addables[i] < addables[smallest] ? i : smallest;
        }
      }
    };
    g.hamming_neighbours(callback);
    for (uint16_t i = 0; i < 10; ++i) {
      sig_val += addables[i];
      bg_val += bddables[i];
    }
    return {sig_val / 11, bg_val / 11};
  }

  double fill_rate() const { return double(size_) / MOD; }
};

}  // namespace sf
