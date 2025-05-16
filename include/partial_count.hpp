#pragma once

#include <stdlib.h>

#include <array>
#include <cstdint>
#include <utility>

namespace sf {

template <class gapmer_t>
class partial_count {
 private:
  static const constexpr uint64_t MOD = 59999999;
  static const constexpr uint64_t STEP = 127;
  struct elem {
    gapmer_t mer;
    uint64_t sig_count;
    uint64_t bg_count;
  };

  elem* data_;
  uint64_t size_;

  uint64_t find(gapmer g) {
    uint64_t g_val(g);
    uint64_t idx = g_val % MOD;
    uint64_t val(data_[idx].mer);
    while (val > 0 && val != g_val) {
      idx += STEP;
      idx -= idx >= MOD ? MOD : 0;
      val = data_[idx].mer;
    }
    return idx;
  }

 public:
  partial_count() : size(0) { data_ = (elem*)calloc(MOD, sizeof(elem)); }

  ~partial_count() {
    if (data_ != nullptr) {
      free(data_);
    }
  }

  void init(gapmer_t g) {
    uint64_t idx = find(g);
    if (data_[idx].mer != g) {
      data_[idx] = {g, 0};
      ++size_;
    }
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

  std::pair<uint64_t, uint64_t> count(gapmer_t g) {
    uint64_t idx = find(g);
    if (data_[idx].mer == g) {
      return data[idx].count;
    }
    return 0;
  }

  std::pair<double, double> smooth_count(gapmer_t g) {
    uint64_t idx = find(g);
    double sig_val = data_[idx].sig_count;
    double bg_val = data_[idx].bg_count;
    std::array<uint64_t, 10> addables{};
    std::array<uint64_t, 10> bddables{};
    uint16_t smallest = 0;
    auto callback = [&](gapmer o) {
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
};

}  // namespace sf
