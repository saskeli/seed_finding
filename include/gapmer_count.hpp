#pragma once

#include <stdlib.h>

#include <SeqIO/SeqIO.hh>
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <sdsl/bit_vectors.hpp>
#include <string>
#include <string_view>
#include <utility>

#include "gapmer.hpp"
#include "string_buffer.hpp"

namespace sf {

template <bool middle_gap_only, uint8_t max_gap>
class gapmer_count {
 public:
  const static constexpr uint64_t ONE = 1;

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

  double* sig_counts;
  double* bg_counts;
  sdsl::bit_vector discarded;
  uint8_t k_;

 private:
  static void count_mers(std::string& fasta_path, double* counts, uint8_t k) {
    seq_io::Reader_x r(fasta_path);
    r.enable_reverse_complements();
#pragma omp parallel
    while (true) {
      uint64_t len{};
      string_buffer<uint64_t> buffer;
#pragma omp critical
      {
        len = r.get_next_read_to_buffer();
        buffer = std::string_view{r.read_buf, len};
      }
      if (len == 0) {
        break;
      }
      if (len < k) {
        continue;
      }
      std::string_view const ss{buffer};
      gapmer<middle_gap_only, max_gap> g(buffer.data(), k);
      uint64_t cv = g.value();
#pragma omp atomic
      counts[cv] = counts[cv] + 1;
      for (uint32_t next = k; next < len; ++next) {
        g = g.next(ss[next]);
        cv = g.value();
#pragma omp atomic
        counts[cv] = counts[cv] + 1;
      }
      uint8_t gap_s = middle_gap_only ? k / 2 : 1;
      uint8_t gap_lim = middle_gap_only ? (k + 3) / 2 : k;
      for (; gap_s < gap_lim; ++gap_s) {
        for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
          if (k + gap_l >= len) {
            break;
          }
          uint64_t off = offset(k, gap_s, gap_l);
          g = {buffer.data(), k, gap_s, gap_l};
          cv = off + g.value();
#pragma omp atomic
          counts[cv] = counts[cv] + 1;
          for (uint32_t mid = gap_s, next = k + gap_l; next < len;
               ++next, ++mid) {
            g = g.next(ss[mid], ss[next]);
            cv = off + g.value();
#pragma omp atomic
            counts[cv] = counts[cv] + 1;
          }
        }
      }
    }
  }

 public:
  gapmer_count(std::string& sig_fasta_path, std::string& bg_fasta_path,
               uint8_t k)
      : sig_counts((double*)calloc(lookup_elems(k), sizeof(double))),
        bg_counts((double*)calloc(lookup_elems(k), sizeof(double))),
        discarded(lookup_elems(k)),
        k_(k) {
    count_mers(sig_fasta_path, sig_counts, k_);
    count_mers(bg_fasta_path, bg_counts, k_);
  }

  gapmer_count()
      : sig_counts(nullptr), bg_counts(nullptr), discarded(0), k_(0) {}

  gapmer_count(gapmer_count&& other) : discarded(0) {
    k_ = std::exchange(other.k_, 0);
    sig_counts = std::exchange(other.sig_counts, nullptr);
    bg_counts = std::exchange(other.bg_counts, nullptr);
    discarded = std::exchange(other.discarded, discarded);
  }

  gapmer_count(const gapmer_count&) = delete;

  gapmer_count& operator=(gapmer_count&& other) {
    std::swap(other.k_, k_);
    std::swap(other.sig_counts, sig_counts);
    std::swap(other.bg_counts, bg_counts);
    std::swap(other.discarded, discarded);
    return *this;
  }

  gapmer_count& operator=(const gapmer_count&) = delete;

 private:
  void smooth(double* arr, double* brr, double* sig_scratch, double* bg_scratch,
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
    double* sig_scratch = (double*)calloc(v_lim, sizeof(double));
    double* bg_scratch = (double*)calloc(v_lim, sizeof(double));
    smooth(sig_counts, bg_counts, sig_scratch, bg_scratch, v_lim);
    std::memcpy(sig_counts, sig_scratch, v_lim * sizeof(double));
    std::memcpy(bg_counts, bg_scratch, v_lim * sizeof(double));
    uint8_t gap_s = middle_gap_only ? k_ / 2 : 1;
    uint8_t gap_lim = middle_gap_only ? (k_ + 3) / 2 : k_;
    for (; gap_s < gap_lim; ++gap_s) {
      for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
        uint64_t off = offset(gap_s, gap_l);
        double* l_count = sig_counts + off;
        double* r_count = bg_counts + off;
        std::memset(sig_scratch, 0, v_lim * sizeof(double));
        std::memset(bg_scratch, 0, v_lim * sizeof(double));
        smooth(l_count, r_count, sig_scratch, bg_scratch, v_lim);
        std::memcpy(l_count, sig_scratch, v_lim * sizeof(double));
        std::memcpy(r_count, bg_scratch, v_lim * sizeof(double));
      }
    }
    free(sig_scratch);
    free(bg_scratch);
  }

  template <class gapmer>
  std::pair<double, double> count(gapmer g) const {
    uint64_t off = offset(g.gap_start(), g.gap_length());
    off += g.value();
    return {sig_counts[off], bg_counts[off]};
  }

  uint64_t offset(uint8_t gap_s, uint8_t gap_l) const {
    return offset(k_, gap_s, gap_l);
  }

  ~gapmer_count() {
    if (sig_counts != nullptr) {
      free(sig_counts);
    }
    if (bg_counts != nullptr) {
      free(bg_counts);
    }
  }
};
}  // namespace sf
