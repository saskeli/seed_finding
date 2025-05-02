#pragma once

#include <stdlib.h>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <utility>

#include "SeqIO/SeqIO.hh"
#include "gapmer.hpp"
#include "sdsl/bit_vectors.hpp"

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
    ret += (ret + 63) / 64;
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

  uint64_t* counts;
  sdsl::bit_vector discarded;
  uint8_t k_;

  gapmer_count(const char* fasta_path, uint8_t k)
      : counts((uint64_t*)calloc(lookup_elems(k), sizeof(uint64_t))),
        discarded(lookup_elems(k)),
        k_(k) {
    seq_io::Reader r(fasta_path);
    r.enable_reverse_complements();
    while (true) {
      uint64_t len = r.get_next_read_to_buffer();
      if (len == 0) {
        break;
      }
      gapmer<middle_gap_only, max_gap> g(r.read_buf, k_);
      uint64_t cv = g.value();
      counts[cv] = counts[cv] + 1;
      for (uint32_t next = k_; next < len; ++next) {
        g = g.next(r.read_buf[next]);
        cv = g.value();
        counts[cv] = counts[cv] + 1;
      }
      uint8_t gap_s = middle_gap_only ? k_ / 2 : 1;
      uint8_t gap_lim = middle_gap_only ? (k_ + 3) / 2 : k_;
      for (; gap_s < gap_lim; ++gap_s) {
        for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
          if (k_ + gap_l >= len) {
            break;
          }
          uint64_t off = offset(gap_s, gap_l);
          g = {r.read_buf, k, gap_s, gap_l};
          cv = off + g.value();
          counts[cv] = counts[cv] + 1;
          for (uint32_t mid = gap_s, next = k_ + gap_l; next < len;
               ++next, ++mid) {
            g = g.next(r.read_buf[mid], r.read_buf[next]);
            cv = off + g.value();
            counts[cv] = counts[cv] + 1;
          }
        }
      }
    }
  }

  gapmer_count() : counts(nullptr), discarded(0), k_(0) {}

  gapmer_count(gapmer_count&& other) : discarded(0) {
    k_ = std::exchange(other.k_, 0);
    counts = std::exchange(other.counts, nullptr);
    discarded = std::exchange(other.discarded, discarded);
  }

  gapmer_count(const gapmer_count&) = delete;

  gapmer_count& operator=(gapmer_count&& other) {
    std::swap(other.k_, k_);
    std::swap(other.counts, counts);
    std::swap(other.discarded, discarded);
    return *this;
  }

  gapmer_count& operator=(const gapmer_count&) = delete;

  template <class gapmer>
  uint64_t count(gapmer g) const {
    uint64_t off = offset(g.gap_start(), g.gap_length());
    return counts[off + g.value()];
  }

  uint64_t offset(uint8_t gap_s, uint8_t gap_l) const {
    return offset(k_, gap_s, gap_l);
  }

  ~gapmer_count() {
    if (counts != nullptr) {
      free(counts);
    }
  }
};
}  // namespace sf
