#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <SeqIO/SeqIO.hh>
#include <string>
#include <string_view>
#include <utility>

#include "string_buffer.hpp"

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

  uint64_t find(gapmer_t g) const {
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
  partial_count() : size_(0) { data_ = (elem*)calloc(MOD, sizeof(elem)); }

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
  void count_mers(const std::string& sig_path, const std::string& bg_path,
                  uint16_t k) {
    seq_io::Reader_x sr(sig_path);
    sr.enable_reverse_complements();
#pragma omp parallel
    while (true) {
      uint64_t len{};
      string_buffer<uint64_t> buffer;
#pragma omp critical
      {
        len = sr.get_next_read_to_buffer();
        buffer = std::string_view{sr.read_buf, len};
      }
      if (len == 0) {
        break;
      }
      if (len < k) {
        continue;
      }
      std::string_view const ss{buffer};
      gapmer_t g(buffer.data(), k);
      inc_sig(g);
      for (uint32_t next = k; next < len; ++next) {
        g = g.next(ss[next]);
        inc_sig(g);
      }
      uint8_t gap_s = middle_gap_only ? k / 2 : 1;
      uint8_t gap_lim = middle_gap_only ? (k + 3) / 2 : k;
      for (; gap_s < gap_lim; ++gap_s) {
        for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
          if (k + gap_l >= len) {
            break;
          }
          g = {buffer.data(), uint8_t(k), gap_s, gap_l};
          inc_sig(g);
          for (uint32_t mid = gap_s, next = k + gap_l; next < len;
               ++next, ++mid) {
            g = g.next(ss[mid], ss[next]);
            inc_sig(g);
          }
        }
      }
    }
    seq_io::Reader_x br(bg_path);
    br.enable_reverse_complements();
#pragma omp parallel
    while (true) {
      uint64_t len{};
      string_buffer<uint64_t> buffer;
#pragma omp critical
      {
        len = br.get_next_read_to_buffer();
        buffer = std::string_view{br.read_buf, len};
      }
      if (len == 0) {
        break;
      }
      std::string_view const ss{buffer};
      gapmer_t g(buffer.data(), k);
      inc_bg(g);
      for (uint32_t next = k; next < len; ++next) {
        g = g.next(ss[next]);
        inc_bg(g);
      }
      uint8_t gap_s = middle_gap_only ? k / 2 : 1;
      uint8_t gap_lim = middle_gap_only ? (k + 3) / 2 : k;
      for (; gap_s < gap_lim; ++gap_s) {
        for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
          if (k + gap_l >= len) {
            break;
          }
          g = {buffer.data(), uint8_t(k), gap_s, gap_l};
          inc_bg(g);
          for (uint32_t mid = gap_s, next = k + gap_l; next < len;
               ++next, ++mid) {
            g = g.next(ss[mid], ss[next]);
            inc_bg(g);
          }
        }
      }
    }
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

  double fill_rate() const {
    return double(size_) / MOD;
  }
};

}  // namespace sf
