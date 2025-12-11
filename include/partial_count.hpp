#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <libbio/algorithm.hh>
#include <utility>
#include "count_base.hpp"

namespace sf {

// t_base used only for having a type name for the base class.
template <class t_gapmer, typename t_base = count_base <t_gapmer>>
class partial_count final : public t_base {
 public:
  typedef t_base base_type;
  typedef base_type::gapmer_type gapmer_type;
  using base_type::max_gap;
  using base_type::middle_gap_only;
  using base_type::count_mers;

 private:
  static const constexpr uint64_t MOD = 59999999;
  static const constexpr uint64_t STEP = 40960001;
  struct elem {
    t_gapmer mer;
    uint64_t sig_count;
    uint64_t bg_count;
  };

  elem* data_;
  uint64_t size_;

  uint64_t find(t_gapmer g) const {
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

  void init(gapmer_type g) {
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

  void increment_signal_count(gapmer_type gg) override {
    uint64_t idx = find(gg);
    if (data_[idx].mer == gg) {
#pragma omp atomic relaxed
      ++data_[idx].sig_count;
    }
  }

  void increment_background_count(gapmer_type gg) override {
    uint64_t idx = find(gg);
    if (data_[idx].mer == gg) {
#pragma omp atomic relaxed
      ++data_[idx].bg_count;
    }
  }

  uint64_t offset(uint8_t, uint8_t) const override { return 0; }

  void increment_signal_count_gapped(gapmer_type gg, uint64_t) override {
    increment_signal_count(gg);
  }

  void increment_background_count_gapped(gapmer_type gg, uint64_t) override {
    increment_background_count(gg);
  }

  std::pair<uint64_t, uint64_t> count(t_gapmer g) const {
    uint64_t idx = find(g);
    if (data_[idx].mer == g) {
      return {data_[idx].sig_count, data_[idx].bg_count};
    }
    return {0, 0};
  }

  std::pair<double, double> smooth_count(t_gapmer g) const {
    uint64_t idx = find(g);
    double sig_val = data_[idx].sig_count;
    double bg_val = data_[idx].bg_count;
    std::array<uint64_t, 10> addables{};
    std::array<uint64_t, 10> bddables{};
    uint16_t smallest = 0;
    auto callback = [&](t_gapmer o) {
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
