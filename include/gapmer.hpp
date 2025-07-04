#pragma once

#include <unistd.h>

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>

#include "bits.hpp"
#include "util.hpp"

#ifdef DEBUG
#include <bitset>
#include <cassert>
#include <iostream>
#endif

namespace sf {
template <bool middle_gap_only = false, uint16_t t_max_gap = 10>
class gapmer {
 public:
  const static constexpr uint64_t max_k = 24;
  const static constexpr uint16_t max_gap = t_max_gap;

 private:
  const static constexpr uint64_t ONE = 1;
  const static constexpr uint64_t value_mask = (ONE << (max_k * 2)) - 1;
  const static constexpr uint64_t meta_mask = 0b11111;
  const static constexpr uint8_t meta_length = 5;

  static_assert(meta_length);
  static_assert(max_gap <= (1 << meta_length) - 1);

  // data_ formatted as follows. We denote M = max_k. The length does not
  // include the gap length. The value is stored such that the rightmost
  // character is at position zero, the next one at position 2 and so on. Note
  // that the padding must be zeroed to make operator== etc. work.
  //
  //  0                    2M           2M + 5       2M + 10      2M + 15  63
  // +--------------------+------------+------------+------------+-----------+
  // | value              | length     | gap_start  | gap_length | padding   |
  // +--------------------+------------+------------+------------+-----------+
  uint64_t data_{};

  template <bool t_has_gap>
  static inline uint64_t read_word_aligned_characters(uint64_t const* d_ptr,
                                                      uint8_t kk,
                                                      uint8_t gap_start = 0,
                                                      uint8_t gap_length = 0);

  static inline uint64_t from_packed_characters(
      uint64_t data,
      uint8_t kk);  //< Returns the packed representation including the length.
  static inline uint64_t from_packed_characters(
      uint64_t data, uint8_t kk, uint8_t gap_start,
      uint8_t gap_length);  //< Returns the packed representation including the
                            // lengths and the gap position.

  explicit gapmer(uint64_t data)
      : data_(data) {}  //< Construct from the given data.
  gapmer(uint64_t prefix, uint64_t suffix, uint8_t p_len, uint8_t s_len,
         uint8_t gap_s, uint8_t gap_l);  //< Construct from the given prefix,
                                         // suffix, gap position and lengths.

  uint64_t prefix() const;
  uint64_t suffix() const;
  uint64_t prefix(uint64_t i) const;
  uint64_t suffix(uint64_t i) const;
  gapmer hamming(uint64_t v, auto& callback) const;

 public:  // FIXME: all_gap_neighbours, middle_gap_neighbours should be private.
          // For now they are public so that unit tests can access them.
  template <bool no_smaller, bool no_same, bool no_larger>
  void middle_gap_neighbours(auto&& callback) const;

  template <bool no_smaller, bool no_same, bool no_larger>
  void all_gap_neighbours(auto&& callback) const;

  /// Writes the 2-bit compressed string representation to the given buffer.
  /// When the span is dynamically sized, the caller is responsible for
  /// reserving enough memory.
  template <std::size_t t_n>
  void write_2bit_coded_to_buffer(std::span<uint64_t, t_n> span) const
    requires((2U * (max_k + max_gap) + 7U) / 8U <= t_n);

 public:
  constexpr gapmer() = default;  //< Construct an empty value.
  gapmer(uint64_t v, uint8_t k)
      : data_(from_packed_characters(v, k)) {
  }  //< Construct from the given packed data and lengths.
  gapmer(uint64_t v, uint8_t k, uint8_t gap_start, uint8_t gap_length)
      : data_(from_packed_characters(v, k, gap_start, gap_length)) {
  }  //< Construct from the given packed data, gap position and lengths.

  /// Construct from the given character data ([ACGT.]*) aligned to 8 bytes.
  gapmer(uint64_t const* d_ptr, uint8_t kk)
      : gapmer(read_word_aligned_characters<false>(d_ptr, kk), kk) {}

  /// Construct from the given character data ([ACGT.]*) aligned to 8 bytes.
  gapmer(uint64_t const* d_ptr, uint8_t kk, uint8_t gap_start,
         uint8_t gap_length)
      : gapmer((0 == gap_length ? read_word_aligned_characters<false>(d_ptr, kk)
                                : read_word_aligned_characters<true>(
                                      d_ptr, kk, gap_start, gap_length)),
               kk, gap_start, gap_length) {}

  uint64_t data() const { return data_; }       //< Get the packed data.
  operator uint64_t() const { return data(); }  //< Get the packed data.
  bool operator==(const gapmer& rhs) const {
    return data_ == rhs.data_;
  }  //< Compare packed bytes.
  bool operator!=(const gapmer& rhs) const {
    return data_ != rhs.data_;
  }  //< Compare packed bytes.

  template <bool compare_rc = false, bool debug = false>
  bool is_neighbour(const gapmer& other) const;

  /// Check if this aligns to another gapmer.
  /// Returns true iff.
  ///  1. data and lengths match exactly or
  ///  2. other’s length is strictly greater and
  ///     2.1 other’s gap length is non-zero, this’s gap length is zero and this
  ///     is a substring of other’s prefix or suffix or 2.2 gap lengths match
  ///     and this’s prefix is a suffix of other’s prefix and this’s suffix is a
  ///     prefix of other’s suffix or 2.3 both gap lengths are zero and this is
  ///     a substring of other.
  template <bool compare_rc = true, bool debug = false>
  bool aligns_to(gapmer other) const;

  size_t compute_offset(const gapmer& other, int& out) const;

  uint16_t length() const {
    return (data_ >> (max_k * 2)) & meta_mask;
  }  //< Get the count of the defined bases.
  uint16_t gap_start() const;  //< Get the starting position of the gap.
  uint16_t gap_length() const {
    return data_ >> (max_k * 2 + 10);
  }  //< Get the gap length.
  uint8_t nuc(uint8_t i) const;  //< Get the i-th 2-bit encoded nucleotide (gap
                                 // not taken into account).
  uint8_t get_c(uint64_t i)
      const;  //< Get the unpacked character or gap at the i-th position.
  uint64_t value() const {
    return data_ & value_mask;
  }  //< Get the encoded value.
  gapmer next(
      char c) const;  //< For a non-gapped gapmer, returns a copy of it with c
                      // appended and the leftmost character removed. If this is
                      // empty, returns an empty gapmer.
  gapmer next(char c1, char c2)
      const;  //< For a gapped gapmer, returns a copy of it with c1 and c2
              // appended respectively to the prefix and to the suffix, with the
              // leftmost characters removed.
  void hamming_neighbours(auto& callback) const;

  template <bool no_smaller = false, bool no_same = false,
            bool no_larger = false>
  void huddinge_neighbours(auto&& callback) const;

  uint16_t huddinge_distance(
      gapmer const other, int& out) const;  //< Calculate the Huddinge distance.
  std::string to_string() const;  //< Returns the sequence as std::string.
  bool is_canonical()
      const;  //< Returns true iff. this is lexicographically equal or smaller
              // than its reverse complement. // FIXME: is the statement true?
  gapmer reverse_complement() const;  //< Get this’s reverse complement.
  bool is_valid() const;
  std::bitset<64> bits() const { return data_; }
};

template <bool middle_gap_only, uint16_t t_max_gap>
gapmer<middle_gap_only, t_max_gap>::gapmer(uint64_t prefix, uint64_t suffix,
                                           uint8_t p_len, uint8_t s_len,
                                           uint8_t gap_s, uint8_t gap_l) {
  data_ = prefix << (2 * s_len);
  data_ |= suffix;
  data_ |= uint64_t(p_len + s_len) << (max_k * 2);
  data_ |= uint64_t(gap_s) << (max_k * 2 + 5);
  data_ |= uint64_t(gap_l) << (max_k * 2 + 10);
}

template <bool middle_gap_only, uint16_t t_max_gap>
uint64_t gapmer<middle_gap_only, t_max_gap>::from_packed_characters(
    uint64_t data, uint8_t kk) {
  assert(kk <= max_k);
  uint64_t meta{kk};
  meta <<= 2 * max_k;
  return meta | data;
}

template <bool middle_gap_only, uint16_t t_max_gap>
uint64_t gapmer<middle_gap_only, t_max_gap>::from_packed_characters(
    uint64_t data, uint8_t kk, uint8_t gap_start, uint8_t gap_length) {
  assert(kk <= max_k);
  assert((0 == gap_start && 0 == gap_length) ||
         (gap_start < kk && gap_length <= max_gap));

  uint64_t meta{gap_length};
  meta <<= meta_length;
  meta |= gap_start;
  meta <<= meta_length;
  meta |= kk;
  meta <<= 2 * max_k;
  return meta | data;
}

template <bool middle_gap_only, uint16_t t_max_gap>
template <bool t_has_gap>
uint64_t gapmer<middle_gap_only, t_max_gap>::read_word_aligned_characters(
    uint64_t const* d_ptr, uint8_t const kk, uint8_t const gap_start,
    uint8_t const gap_length) {
  assert(kk <= max_k);
  assert(gap_start < kk);
  assert(gap_length <= max_gap);

  uint64_t retval{};
  uint8_t ii{};
  uint8_t const limit((kk + gap_length + 7) /
                      8);  // Since d_ptr is uint64_t const *, we are bound to
                           // have a multiple of 8 bytes.
  if constexpr (t_has_gap) {
    assert(gap_start);
    assert(gap_length);
    uint64_t iv{};
    // The loop condition is an unrolling hint for the compiler.
    // (Checked the output from Compiler Explorer. GCC unrolls the
    // loop and the one later when using -O3. Clang unrolls the loops
    // even with -O2.)
    while (ii < (max_k + max_gap + 7) / 8) {
      iv = bits::read_multiple_dna_characters(d_ptr[ii]);
      if (gap_start < 8 * (ii + 1)) break;
      retval |= iv;
      retval <<= 16;
      ++ii;
    }

    {
      // Remove the gap characters if needed.
      // Higher bits.
      {
        uint64_t const count(gap_start - 8 * ii);
        uint64_t const mask{0xffff & (UINT64_C(0xffff) << (2 * (8 - count)))};
        retval |= iv & mask;
        retval <<= 2 * count;
      }

      // Lower bits.
      if (8 * (ii + 1) <= gap_start + gap_length) {
        ii = (gap_start + gap_length) / 8;
        iv = bits::read_multiple_dna_characters(d_ptr[ii]);
      }

      {
        uint64_t const count(8 * (ii + 1) - gap_start - gap_length);
        uint64_t mask{UINT64_C(0xffff) >> (2 * (8 - count))};
        iv &= mask;

        ++ii;
        if (ii == limit) {
          iv <<= 2 * (8 - count);
          retval |= iv;
          retval >>= 2 * (8 - (kk - gap_start));
          return retval;
        }

        iv <<= 16 - 2 * count;
        retval |= iv;
        retval <<= 2 * count;
      }
    }
  } else {
    assert(0 == gap_start);
    assert(0 == gap_length);
    if (ii == limit) return retval;
  }

  // ii < limit (checked before).
  // The loop condition is an unrolling hint for the compiler.
  while (ii < (max_k + max_gap + 7) / 8) {
    auto iv(bits::read_multiple_dna_characters(d_ptr[ii]));
    retval |= iv;
    ++ii;
    if (ii == limit) break;
    retval <<= 16;
  }

  {
    auto const remaining_count((kk + gap_length) % 8);
    if (remaining_count) retval >>= 2 * (8 - remaining_count);
  }

  // Metadata not encoded.
  return retval;
}

template <bool middle_gap_only, uint16_t t_max_gap>
uint64_t gapmer<middle_gap_only, t_max_gap>::prefix() const {
  uint64_t v = length() - gap_start();
  return value() >> (v * 2);
}

template <bool middle_gap_only, uint16_t t_max_gap>
uint64_t gapmer<middle_gap_only, t_max_gap>::suffix() const {
  uint64_t v = 1;
  v <<= (length() - gap_start()) * 2;
  --v;
  return value() & v;
}

template <bool middle_gap_only, uint16_t t_max_gap>
uint64_t gapmer<middle_gap_only, t_max_gap>::prefix(uint64_t i) const {
  uint64_t v = length() - i;
  return value() >> (v * 2);
}

template <bool middle_gap_only, uint16_t t_max_gap>
uint64_t gapmer<middle_gap_only, t_max_gap>::suffix(uint64_t i) const {
  uint64_t v = 1;
  v <<= i * 2;
  --v;
  return value() & v;
}

template <bool middle_gap_only, uint16_t t_max_gap>
auto gapmer<middle_gap_only, t_max_gap>::hamming(uint64_t v,
                                                 auto& callback) const
    -> gapmer {
  callback(*this);
  callback(gapmer(data_ ^ v));
  uint64_t vv = v << 1;
  callback(gapmer(data_ ^ vv));
  vv |= v;
  callback(gapmer(data_ ^ vv));
  return *this;
}

template <bool middle_gap_only, uint16_t t_max_gap>
template <bool no_smaller, bool no_same, bool no_larger>
void gapmer<middle_gap_only, t_max_gap>::middle_gap_neighbours(
    auto&& callback) const {
  const uint64_t val = value();
  const uint8_t len = length();
  assert(no_larger ||
         len < max_k);  // It is the user’s responsibility to check the length
                        // before listing the neighbours.
  assert(no_smaller || 5 <= len);
  const uint8_t gap_s = gap_start();
  const uint8_t prefix_len = gap_s;
  const uint8_t suffix_len = len - gap_s;
  const uint8_t gap_l = gap_length();

  // Version of callback that checks that the new mer is not the same as the
  // old one
  auto val_cb = [&](gapmer o) {
    if (*this != o) {
      callback(o);
    }
  };

  uint8_t new_len = len + 1;
  uint64_t h_val = ONE << (len * 2);

  if constexpr (no_larger == false) {
    // For defined bases + 1
    // Add one nuc to each gap (implicit gap at start and end)

    // Start:
    // Middle gapped mers can only be added to at the start if
    // length is even or unven length with gap closer to start.

    // ACGTGC -> nACGTGC
    // AC.TGC -> nAC.TGC
    // ACG.TGC -> nACG.TGC
    if (gap_l == 0 || prefix_len == len / 2) {
      uint8_t n_gap_s = gap_s + (gap_l ? 1 : 0);
      gapmer(val, new_len, n_gap_s, gap_l).hamming(h_val, callback);
    }
    // Middle:
    // Middle gaps can be shortened from both ends if length is even
    // and from only one end if length is not even.
    // Obvioustly splitting gaps can't be done here.
    if (gap_l) {
      h_val = ONE << (suffix_len * 2);
      uint64_t n_val = prefix();
      n_val <<= suffix_len * 2 + 2;
      n_val |= suffix();
      if (gap_l == 1) {
        // ACG.TGT -> ACGnTGT
        gapmer(n_val, new_len).hamming(h_val, callback);
      } else if (len % 2 == 0) {
        // ACG..TGT -> ACGn.TGT, ACG.nTGT
        gapmer(n_val, new_len, gap_s + 1, gap_l - 1).hamming(h_val, callback);
        gapmer(n_val, new_len, gap_s, gap_l - 1).hamming(h_val, callback);
      } else if (suffix_len > prefix_len) {
        // AC..TGT -> ACn.TGT
        gapmer(n_val, new_len, gap_s + 1, gap_l - 1).hamming(h_val, callback);
      } else {
        // ACG..TG -> ACG.nTG
        gapmer(n_val, new_len, gap_s, gap_l - 1).hamming(h_val, callback);
      }
    }

    // End:
    // Middle gapped mers can only be added to at the end if
    // length is even or unven length with gap closer to end.
    if (gap_l == 0 || suffix_len == len / 2) {
      // ACGTAT -> ACGTATn
      // ACG.TAT -> ACG.TATn
      // ACG.TA -> ACG.TAn
      gapmer(val << 2, new_len, gap_s, gap_l).hamming(ONE, callback);
    }
  }

  // For defined bases + 0
  new_len = len;
  if constexpr (no_same == false) {
    // Hamming
    for (uint64_t n = 1; n < 4; ++n) {
      uint64_t xor_val = n;
      for (size_t i = 0; i < len; ++i) {
        callback(gapmer(data_ ^ xor_val));
        xor_val <<= 2;
      }
    }

    // Add one to any gap and gap any existing position

    // Start:
    // if a gap exists, add one to start and extend the gap
    // Else create a new gap either in the center or in every
    // position
    h_val = ONE << (len * 2 - 2);
    if (gap_l) {
      uint64_t pref = prefix();
      if (gap_l < max_gap) {
        // AAA.GGG -> nAA..GGG
        // AAA.GG -> nAA..GG
        // AA.GGG -> nA.GGG
        gapmer(pref >> 2, suffix(), prefix_len, suffix_len, gap_s, gap_l + 1)
            .hamming(h_val, callback);
        if (suffix_len > prefix_len) {
          // AA.GGT -> nAA..GT
          gapmer(pref, suffix(suffix_len - 1), prefix_len + 1, suffix_len - 1,
                 gap_s + 1, gap_l + 1)
              .hamming(h_val, callback);
          // AA.GGT -> nAA.GG
        }
      }
      if (suffix_len > prefix_len) {
        gapmer(pref, suffix() >> 2, prefix_len + 1, suffix_len - 1, gap_s + 1,
               gap_l)
            .hamming(h_val, callback);
      }
    } else {
      // ACTGTA -> nACTGT ***can be same***
      gapmer(val >> 2, new_len).hamming(h_val, val_cb);
      // ACTGTA -> nAC.GTA
      // ACTGT -> nA.TGT
      gapmer(prefix(len / 2 - 1), suffix((len + 1) / 2), len / 2, (len + 1) / 2,
             len / 2, 1)
          .hamming(h_val, callback);
      if (len % 2 != 0) {
        // ACTGT -> nAC.GT
        gapmer(prefix(len / 2), suffix(len / 2), len / 2 + 1, len / 2,
               len / 2 + 1, 1)
            .hamming(h_val, callback);
      }
    }

    // Middle:
    // Gapping start and end is allowed
    if (gap_l == 1) {
      // ACG.TGT -> ACGnTG
      h_val = ONE << (suffix_len * 2 - 2);
      gapmer(prefix(), suffix() >> 2, prefix_len, suffix_len, 0, 0)
          .hamming(h_val, callback);
      // ACG.TGT -> CGnTGT
      h_val <<= 2;
      gapmer(suffix(len - 1) >> (suffix_len * 2), suffix(), prefix_len - 1,
             suffix_len + 1, 0, 0)
          .hamming(h_val, callback);
      if (prefix_len < suffix_len) {
        // AC.GTG -> ACn.TG
        h_val >>= 2;
        gapmer(prefix(), suffix(prefix_len), prefix_len, suffix_len, gap_s + 1,
               gap_l)
            .hamming(h_val, callback);
      } else if (suffix_len < prefix_len) {
        // ACG.TG -> AC.nTG
        gapmer(prefix(suffix_len), suffix(), suffix_len, prefix_len, gap_s - 1,
               gap_l)
            .hamming(h_val, callback);
      }
    } else if (gap_l > 1) {
      // ACG...TGT -> ACG..nTG
      // AC...GTG -> AC..nGT
      // ACG...TG -> ACG..nT
      h_val = ONE << (suffix_len * 2 - 2);
      uint64_t pref = prefix();
      uint64_t suf = suffix() >> 2;
      gapmer(pref, suf, prefix_len, suffix_len, gap_s, gap_l - 1)
          .hamming(h_val, callback);
      if (prefix_len < suffix_len) {
        // AC...GTG -> ACn..GT
        gapmer(pref, suf, prefix_len, suffix_len, gap_s + 1, gap_l - 1)
            .hamming(h_val, callback);
        // AC...GTG -> ACn...TG
        gapmer(val, len, gap_s + 1, gap_l).hamming(h_val, callback);
      }
      // ACG...TGT -> CGn..TGT
      // AC...GTG -> Cn..GTG
      // ACG...TG -> CGn..TG
      h_val <<= 2;
      pref <<= 2;
      pref &= (ONE << (prefix_len * 2)) - ONE;
      suf = suffix();
      gapmer(pref, suf, prefix_len, suffix_len, gap_s, gap_l - 1)
          .hamming(h_val, callback);
      if (suffix_len < prefix_len) {
        // ACG...TG -> AC...nTG
        gapmer(val, len, gap_s - 1, gap_l).hamming(h_val, callback);
        // ACG...TG -> CG..nTG
        gapmer(pref >> 2, suf, prefix_len - 1, suffix_len + 1, gap_s - 1,
               gap_l - 1)
            .hamming(h_val, callback);
      }
    }

    // End:
    // if a gap exists, add one to end and extend the gap
    // Else create a new gap either in the center or in every
    // position
    if (gap_l) {
      uint8_t l_g_l = gap_l + 1;
      if (l_g_l <= max_gap) {
        // ACG..TG -> ACG...Gn
        // ACG..TGC -> ACG...GCn
        // AC..TGC -> AC...GCn
        uint64_t suf = suffix(suffix_len - 1) << 2;
        gapmer(prefix(), suf, prefix_len, suffix_len, gap_s, l_g_l)
            .hamming(ONE, callback);
        if (prefix_len > suffix_len) {
          // ACG..TG -> AC...TGn
          uint64_t p = prefix(prefix_len - 1);
          gapmer(p, suffix() << 2, prefix_len - 1, suffix_len + 1, gap_s - 1,
                 l_g_l)
              .hamming(ONE, callback);
        }
      }
      if (prefix_len > suffix_len) {
        // ACG..TG -> CG...TGn
        uint64_t n_val = val << 2;
        n_val &= (ONE << (len * 2)) - ONE;
        gapmer(n_val, len, gap_s - 1, gap_l).hamming(ONE, callback);
      }
    } else {
      // ACGTGC -> CGTGCn
      uint64_t n_val = suffix(len - 1);
      n_val <<= 2;
      gapmer(n_val, len).hamming(ONE, val_cb);
      // ACGTGC -> ACG.GCn
      // ACGTG -> AC.TGn
      gapmer(prefix(len / 2), suffix((len - 1) / 2) << 2, len / 2,
             len - len / 2, len / 2, 1)
          .hamming(ONE, callback);
      if (len % 2 > 0) {
        // ACGTG -> ACG.Gn
        gapmer(prefix(len / 2 + 1), suffix(len / 2 - 1) << 2, len / 2 + 1,
               len / 2, len / 2 + 1, 1)
            .hamming(ONE, callback);
      }
    }
  }

  if constexpr (no_smaller == false) {
    // For defined bases - 1
    // Gap any existing position
    // Start:
    if (gap_s == 0) {
      // ACGTAC -> CGTAC
      callback(gapmer(suffix(len - 1), len - 1));
    } else {
      // ACG..TC -> CG..TC
      // ACG..TCG -> CG..TCG
      if (len / 2 == suffix_len) {
        callback(gapmer(suffix(len - 1), len - 1, gap_s - 1, gap_l));
      }
    }

    // Middle
    // Extend gap at edges or add new gap
    if (gap_l == 0) {
      uint8_t half = len / 2;
      uint8_t rem = half - 1;
      if (len % 2 == 0) {
        // ACGTAC -> ACG.AC, AC.TAC
        callback(gapmer(prefix(half), suffix(rem), half, rem, half, 1));
        callback(gapmer(prefix(rem), suffix(half), rem, half, rem, 1));
      } else {
        // ACGTA -> AC.TA
        callback(gapmer(prefix(half), suffix(half), half, half, half, 1));
      }
    } else {
      if (gap_l < max_gap) {
        uint8_t half = len / 2;
        uint8_t s_l = len - half - 1;
        if (prefix_len == len / 2) {
          // ACG..TGC -> ACG...GC
          // AC..TGC -> AC...GC
          callback(
              gapmer(prefix(half), suffix(s_l), half, s_l, gap_s, gap_l + 1));
        }
        if (suffix_len == len / 2) {
          // ACG..TGC -> AC...TGC
          // ACG..TG -> AC...TG
          callback(gapmer(prefix(s_l), suffix(half), s_l, half, gap_s - 1,
                          gap_l + 1));
        }
      }
    }

    // End
    if (gap_s == 0) {
      // CATTATT -> CATTAT
      callback(gapmer(prefix(len - 1), len - 1));
    } else {
      // CAT..ATT -> CAT..AT
      // CA..TAT -> CA..TA
      if (prefix_len == len / 2) {
        callback(gapmer(prefix(len - 1), len - 1, gap_s, gap_l));
      }
    }
  }
}

template <bool middle_gap_only, uint16_t t_max_gap>
template <bool no_smaller, bool no_same, bool no_larger>
void gapmer<middle_gap_only, t_max_gap>::all_gap_neighbours(
    auto&& callback) const {
  const uint64_t val = value();
  const uint8_t len = length();
  assert(no_larger ||
         len < max_k);  // It is the user’s responsibility to check the length
                        // before listing the neighbours.
  assert(no_smaller || 5 <= len);
  const uint8_t gap_s = gap_start();
  const uint8_t prefix_len = gap_s;
  const uint8_t suffix_len = len - gap_s;
  const uint8_t gap_l = gap_length();

  // Version of callback that checks that the new mer is not the same as the
  // old one
  auto val_cb = [&](gapmer o) {
    if (*this != o) {
      callback(o);
    }
  };

  uint8_t new_len = len + 1;
  uint64_t h_val = ONE << (len * 2);

  if constexpr (no_larger == false) {
    // For defined bases + 1
    // Add one nuc to each gap (implicit gap at start and end)

    // Start:
    // Middle gapped mers can only be added to at the start if
    // length is even or unven length with gap closer to start.
    if (gap_l == 0) {
      // ACGTT -> nACGTT, n.*ACGTT
      gapmer(val, new_len).hamming(h_val, callback);
      for (uint8_t gl = 1; gl <= max_gap; ++gl) {
        gapmer(val, new_len, 1, gl).hamming(h_val, callback);
      }
    } else {
      // ACGT.G -> nACGT.G
      gapmer(val, new_len, gap_s + 1, gap_l).hamming(h_val, callback);
    }

    // Middle:
    // Middle gaps can be shortened from both ends if length is even
    // and from only one end if length is not even.
    // Obvioustly splitting gaps can't be done here.
    if (gap_l) {
      h_val = ONE << (suffix_len * 2);
      uint64_t n_val = prefix();
      n_val <<= suffix_len * 2 + 2;
      n_val |= suffix();
      if (gap_l == 1) {
        // A.CGTGT -> AnCGTGT
        gapmer(n_val, new_len).hamming(h_val, callback);
      } else {
        // A..CGTGT -> An.CGTGT, A.nCGTGT
        gapmer(n_val, new_len, gap_s + 1, gap_l - 1).hamming(h_val, callback);
        gapmer(n_val, new_len, gap_s, gap_l - 1).hamming(h_val, callback);
      }
    }

    // End:
    // Middle gapped mers can only be added to at the end if
    // length is even or unven length with gap closer to end.
    uint64_t n_val = val << 2;
    if (gap_l == 0) {
      // ACGTAT -> ACGTATn, ACGTAT.*n
      gapmer(n_val, new_len).hamming(ONE, callback);
      for (uint8_t gl = 1; gl <= max_gap; ++gl) {
        gapmer(n_val, new_len, len, gl).hamming(ONE, callback);
      }
    } else {
      // ACG...TAT -> ACG..TATn
      // ACG...TA -> ACG...TAn
      // AC...TAT -> AC...TATn
      gapmer(n_val, new_len, gap_s, gap_l).hamming(ONE, callback);
    }
  }

  // For defined bases + 0
  new_len = len;
  if constexpr (no_same == false) {
    // Hamming
    for (uint64_t n = 1; n < 4; ++n) {
      uint64_t xor_val = n;
      for (size_t i = 0; i < len; ++i) {
        callback(gapmer(data_ ^ xor_val));
        xor_val <<= 2;
      }
    }

    // Add one to any gap and gap any existing position

    // Start:
    // if a gap exists, add one to start and extend the gap
    // Else create a new gap either in the center or in every
    // position
    h_val = ONE << (len * 2 - 2);
    if (gap_l) {
      if (gap_s == len - 1) {
        // ACGT.G -> nACGT
        uint64_t n_val = prefix();
        gapmer(n_val, new_len).hamming(h_val, callback);
        // ACGT.G -> n.*ACGT
        for (uint8_t gl = 1; gl <= max_gap; ++gl) {
          gapmer(n_val, new_len, 1, gl).hamming(h_val, callback);
        }
      } else {
        // ACG.TGT -> nACG.TG
        gapmer(val >> 2, new_len, gap_s + 1, gap_l).hamming(h_val, callback);
        if (gap_s == 1) {
          // A.CGTG -> n..*CGTG
          uint64_t suf = suffix();
          for (uint8_t gl = gap_l + 1; gl <= max_gap; ++gl) {
            gapmer(suf, len, 1, gl).hamming(h_val, callback);
          }
        }
      }
      if (gap_l < max_gap) {
        if (gap_s == 1) {
          // A.CGTG -> n..CGTG
          gapmer(prefix() >> 2, suffix(), prefix_len, suffix_len, 1, gap_l + 1)
              .hamming(h_val, callback);
          // A.CGTG -> nA..GTG
          gapmer(prefix(), suffix(suffix_len - 1), prefix_len + 1,
                 suffix_len - 1, 2, gap_l + 1)
              .hamming(h_val, callback);
        } else if (gap_s == len - 1) {
          // ACGT.G -> nACG..G
          gapmer(prefix(prefix_len - 1), suffix(), prefix_len, 1, gap_s,
                 gap_l + 1)
              .hamming(h_val, callback);
        } else {
          // ACG.TGT -> nACG..GT
          gapmer(prefix(), suffix(suffix_len - 1), prefix_len + 1,
                 suffix_len - 1, gap_s + 1, gap_l + 1)
              .hamming(h_val, callback);
          // ACG.TGT -> nAC..TGT
          gapmer(prefix(prefix_len - 1), suffix(suffix_len), prefix_len,
                 suffix_len, gap_s, gap_l + 1)
              .hamming(h_val, callback);
        }
      }
    } else {
      // ACGTGT -> nACGTG
      uint64_t n_val = prefix(len - 1);
      gapmer(n_val, new_len).hamming(h_val, val_cb);
      // ACGTGT -> n.*CGTGT
      for (uint8_t gl = 1; gl <= max_gap; ++gl) {
        gapmer(n_val, new_len, 1, gl).hamming(h_val, callback);
      }
      n_val = suffix(len - 1);
      for (uint8_t gl = 1; gl <= max_gap; ++gl) {
        gapmer(n_val, new_len, 1, gl).hamming(h_val, callback);
      }
      // ACGTGT -> nA.GTGT, nAC.TGT, nACG.GT, nACGT.T
      for (uint8_t gs = 2; gs < len; ++gs) {
        // TODO: This is slower than it needs to be.
        gapmer(prefix(gs - 1), suffix(len - gs), gs, len - gs, gs, 1)
            .hamming(h_val, callback);
      }
    }

    // Middle:
    // Gapping start and end is allowed
    if (gap_l == 1) {
      if (gap_s == 1) {
        // A.GTGC -> nGTGC
        h_val = ONE << (len * 2 - 2);
        gapmer(val, len).hamming(h_val, callback);
        // A.GTGC -> An.TGC
        h_val >>= 2;
        gapmer(val, len, 2, 1).hamming(h_val, callback);
        // A.GTGC -> AnG.GC, AnGT.C
        uint64_t pref = prefix() << 4;
        uint8_t n_loc = 1;
        pref |= nuc(n_loc++);
        uint64_t suf = suffix(suffix_len - 2);
        uint8_t s_len = suffix_len - 2;
        uint64_t s_mask = (ONE << (s_len * 2)) - 1;
        for (; s_len > 0; --s_len) {
          gapmer(pref, suf, len - s_len, s_len, len - s_len, 1)
              .hamming(h_val, callback);
          pref <<= 2;
          pref |= nuc(n_loc++);
          s_mask >>= 2;
          suf &= s_mask;
        }
        // A.GTGC -> AnGTG
        gapmer(prefix() << 2, suffix() >> 2, 2, len - 2, 0, 0)
            .hamming(h_val, callback);
      } else if (gap_s == len - 1) {
        // ACGTG.C -> ACGTGn
        gapmer(val, len).hamming(ONE, callback);
        // ACGTG.C -> ACGT.nC
        h_val = ONE << 2;
        gapmer(val, len, gap_s - 1, 1).hamming(h_val, callback);
        // ACGTG.C -> A.GTGnC, AC.TGnC, ACG.GnC
        uint64_t pref = prefix(len - 2);
        uint64_t suf = nuc(len - 1) | (nuc(len - 2) << 4);
        uint8_t sl = 3;
        for (uint8_t pl = len - 3; pl > 0; --pl) {
          gapmer(pref >> 2, suf, pl, sl, pl, 1).hamming(h_val, callback);
          suf |= (pref & 0b11) << (sl * 2);
          pref >>= 2;
          ++sl;
        }
        gapmer(suf, len).hamming(h_val, callback);
        // ACGTG.C -> CGTGnC
      } else {
        // AGC.GATG -> GCnGATG
        h_val = ONE << (suffix_len * 2);
        uint64_t pref = prefix() << 2;
        pref &= (ONE << (prefix_len * 2)) - ONE;
        gapmer(pref, suffix(), prefix_len, suffix_len, 0, 0)
            .hamming(h_val, callback);
        // AGC.GATG -> AG.nGATG
        gapmer(val, len, gap_s - 1, 1).hamming(h_val, callback);
        // AGC.GATG -> AGCn.ATG
        h_val >>= 2;
        gapmer(val, len, gap_s + 1, 1).hamming(h_val, callback);
        // AGC.GATG -> AGCnGAT
        gapmer(prefix(), suffix() >> 2, prefix_len, suffix_len, 0, 0)
            .hamming(h_val, callback);
        // AGC.GATG -> A.CnGATG, AGCnG.TG, AGCnGA.G
        h_val <<= 2;
        uint64_t suf = suffix() | (prefix() << (suffix_len * 2 + 2));
        uint64_t suf_m = (ONE << (len * 2)) - 1;
        suf &= suf_m;
        pref = nuc(0);
        for (uint8_t gs = 1; gs < len; ++gs) {
          uint64_t mc = suf >> (2 * (len - gs));
          suf_m >>= 2;
          suf &= suf_m;
          if (gs < gap_s - 1 || gs > gap_s + 1) {
            gapmer(pref, suf, gs, len - gs, gs, 1).hamming(h_val, callback);
          } else if (gs == gap_s) {
            h_val >>= 2;
          }
          pref <<= 2;
          pref |= mc;
        }
      }
    } else if (gap_l > 1) {
      if (gap_s == 1) {
        // A...CGTGT -> n..CGTGT, n.CGTGT, nCGTGT
        h_val = ONE << (suffix_len * 2);
        for (uint8_t gl = gap_l - 1; gl > 0; --gl) {
          gapmer(val, len, 1, gl).hamming(h_val, callback);
        }
        gapmer(val, len, 0, 0).hamming(h_val, callback);
        // A...CGTGT -> An...GTGT
        h_val >>= 2;
        gapmer(val, len, 2, gap_l).hamming(h_val, callback);
        // A...CGTGT -> An..CGTG, A..nCGTG
        uint64_t suf = suffix() >> 2;
        suf |= uint64_t(nuc(0)) << (suffix_len * 2);
        gapmer(suf, len, 1, gap_l - 1).hamming(h_val, callback);
        gapmer(suf, len, 2, gap_l - 1).hamming(h_val, callback);
      } else if (gap_s == len - 1) {
        // ACGTG...T -> ACGTGn, ACGTG.n, ACGTG..n
        gapmer(val, len, 0, 0).hamming(ONE, callback);
        for (uint8_t gl = 1; gl < gap_l; ++gl) {
          gapmer(val, len, len - 1, gl).hamming(ONE, callback);
        }
        h_val = ONE << 2;
        // ACGTG...T -> CGTGn..T, CGTG..nT, ACGT...nT
        gapmer(prefix(prefix_len - 1), suffix(), prefix_len - 1, 2, gap_s - 1,
               gap_l)
            .hamming(h_val, callback);
        // ACGTG...T -> CGTGn..T, CGTG..nT
        uint64_t pref = prefix() << 4;
        pref |= nuc(len - 1);
        pref &= (ONE << (len * 2)) - ONE;
        gapmer(pref, len, gap_s, gap_l - 1).hamming(h_val, callback);
        gapmer(pref, len, gap_s - 1, gap_l - 1).hamming(h_val, callback);
      } else {
        // ACG...TGT -> ACGn...GT
        h_val = ONE << (suffix_len * 2 - 2);
        gapmer(val, len, gap_s + 1, gap_l).hamming(h_val, callback);
        // ACG...TGT -> AC...nTGT
        h_val <<= 2;
        gapmer(val, len, gap_s - 1, gap_l).hamming(h_val, callback);
        // ACG...TGT -> CGn..TGT, CG..nTGT
        uint64_t n_val = prefix() << (suffix_len * 2 + 2);
        n_val |= suffix();
        n_val &= (ONE << (len * 2)) - ONE;
        gapmer(n_val, len, gap_s, gap_l - 1).hamming(h_val, callback);
        gapmer(n_val, len, gap_s - 1, gap_l - 1).hamming(h_val, callback);
        // ACG...TGT -> ACGn..TG, ACG..nTG
        n_val = prefix() << (suffix_len * 2);
        n_val |= suffix() >> 2;
        h_val >>= 2;
        gapmer(n_val, len, gap_s + 1, gap_l - 1).hamming(h_val, callback);
        gapmer(n_val, len, gap_s, gap_l - 1).hamming(h_val, callback);
      }
    }

    // End:
    // if a gap exists, add one to end and extend the gap
    // Else create a new gap either in the center or in every
    // position
    if (gap_l) {
      if (gap_s == 1) {
        // A..CGTA -> CGTAn, A...GTAn, CGTA.*n
        uint64_t s = suffix() << 2;
        uint64_t p = nuc(0);
        uint64_t s_s = s & ((ONE << (suffix_len * 2)) - 1);
        gapmer(s, len).hamming(ONE, callback);
        for (uint8_t i = 1; i <= max_gap; ++i) {
          gapmer(s, len, uint8_t(len - 1), i).hamming(ONE, callback);
        }
        if (gap_l < max_gap) {
          gapmer(p, s_s, 1, suffix_len, 1, uint8_t(gap_l + 1))
              .hamming(ONE, callback);
        }
      } else if (gap_s == len - 1) {
        // ACGT.A -> ACG..An, ACGT..n, ACGT..*n
        uint64_t p = prefix() << 2;
        for (uint8_t i = gap_l + 1; i <= max_gap; ++i) {
          gapmer(p, len, gap_s, i).hamming(ONE, callback);
        }
        if (gap_l < max_gap) {
          gapmer(prefix() >> 2, suffix() << 2, prefix_len - 1, 2, gap_s - 1,
                 gap_l + 1)
              .hamming(ONE, callback);
        }
        // ACGT.A -> CGT.An
        gapmer(suffix(len - 1) << 2, new_len, gap_s - 1, gap_l)
            .hamming(ONE, callback);
      } else {
        // ACG..TGC -> AC...TGCn, ACG...GCn
        uint64_t p = prefix();
        uint64_t s = suffix() << 2;
        uint64_t s_s = s & ((ONE << (suffix_len * 2)) - 1);
        uint8_t l_g_l = gap_l + 1;
        if (l_g_l <= max_gap) {
          gapmer(p >> 2, s, prefix_len - 1, suffix_len + 1, gap_s - 1, l_g_l)
              .hamming(ONE, callback);
          gapmer(p, s_s, prefix_len, suffix_len, gap_s, l_g_l)
              .hamming(ONE, callback);
        }
        // ACG..TGC -> CG..TGCn
        gapmer(suffix(len - 1) << 2, len, gap_s - 1, gap_l)
            .hamming(ONE, callback);
      }
    } else {
      uint64_t n_val = suffix(len - 1) << 2;
      uint64_t suf = n_val;
      uint64_t pref = val & ~uint64_t(0b11);
      // ACGTGC -> CGTGCn
      gapmer(suf, len).hamming(ONE, val_cb);
      // ACGTGC -> CGTGC.*n ACGTG.*n
      for (uint8_t i = 1; i <= max_gap; ++i) {
        gapmer(suf, len, len - 1, i).hamming(ONE, callback);
        gapmer(pref, len, len - 1, i).hamming(ONE, callback);
      }
      // ACGTGC -> A.GTGCn, AC.TGCn, ACG.GCn, ACGT.Cn
      uint64_t p = nuc(0);
      uint8_t s_l = len - 2;
      uint64_t s_l_mask = ONE << (len * 2 - 2);
      --s_l_mask;
      for (uint8_t i = 1; i < len - 1; ++i) {
        uint64_t g_c = suf >> (2 * s_l + 2);
        suf &= s_l_mask;
        gapmer(p, suf, i, uint8_t(len - i), i, 1).hamming(ONE, callback);
        p <<= 2;
        p |= g_c;
        s_l_mask >>= 2;
        --s_l;
      }
    }
  }

  if constexpr (no_smaller == false) {
    // For defined bases - 1
    // Gap any existing position
    // Start:
    if (gap_s == 1) {
      // A.CGTA -> CGTA
      callback(gapmer(suffix(), suffix_len));
    } else if (gap_s == 0) {
      // ACGTAC -> CGTAC
      callback(gapmer(suffix(len - 1), len - 1));
    } else {
      // AC..GTC -> C..GTC
      callback(gapmer(suffix(len - 1), len - 1, gap_s - 1, gap_l));
    }

    // Middle
    // Extend gap at edges or add new gap
    if (gap_l == 0) {
      // ACGTAC -> A.GTAC, AC.TAC, ACG.AC, ACGT.C
      uint8_t rem = len - 2;
      for (uint8_t i = 1; i < len - 1; ++i) {
        callback(gapmer(prefix(i), suffix(rem), i, rem, i, 1));
        --rem;
      }
    } else {
      if (gap_l < max_gap) {
        if (gap_s > 1) {
          // AA..TGC -> A...TGC
          callback(gapmer(prefix(prefix_len - 1), suffix(), prefix_len - 1,
                          suffix_len, gap_s - 1, gap_l + 1));
        }
        if (gap_s < len - 1) {
          // AAT..GC -> AAT...C
          callback(gapmer(prefix(), suffix(suffix_len - 1), prefix_len,
                          suffix_len - 1, gap_s, gap_l + 1));
        }
      }
    }
    // End
    if (gap_s == len - 1) {
      // CATTA.A -> CATTA
      callback(gapmer(prefix(), prefix_len));
    } else if (gap_s == 0) {
      // CATTATT -> CATTAT
      callback(gapmer(prefix(len - 1), len - 1));
    } else {
      // CAT..ATT -> CAT..AT
      // CA..TAT -> CA..TA
      // CAT..AT -> CAT..A
      callback(gapmer(prefix(len - 1), len - 1, gap_s, gap_l));
    }
  }
}

template <bool middle_gap_only, uint16_t t_max_gap>
template <bool compare_rc, bool debug>
bool gapmer<middle_gap_only, t_max_gap>::is_neighbour(
    const gapmer& other) const {
  if (other.data_ == data_) {
    return false;
  }
  gapmer a = other;
  gapmer b = *this;
  if (length() > other.length()) {
    a = *this;
    b = other;
  }
  uint64_t a_len = a.length();
  uint64_t b_len = b.length();
  if constexpr (debug) {
    std::cerr << "Comparing " << a.to_string() << " to " << b.to_string()
              << "\n"
              << "   len a: " << a_len << ", len b: " << b_len << std::endl;
  }
  if (a_len - b_len > 1) {
    return false;
  }
  uint16_t a_gl = a.gap_length();
  uint16_t a_gs = a.gap_start();
  uint16_t b_gl = b.gap_length();
  uint16_t b_gs = b.gap_start();
  if constexpr (debug) {
    std::cerr << "   gs a: " << a_gs << ", gl a: " << a_gl << "\n"
              << "   gs b: " << b_gs << ", gl b: " << b_gl << std::endl;
  }
  // a_len - 1 bases of b must align with bases of a.
  // The first base of b matches either the
  // first or second base of a...
  uint64_t matches = 0;
  for (uint64_t i = 0; i < b_len; ++i) {
    uint64_t ii = i + (i >= b_gs) * b_gl;
    auto ac = a.get_c(ii);
    auto bc = b.get_c(ii);
    if constexpr (debug) {
      std::cerr << i << ", " << ii << ": " << ac << ", " << bc << std::endl;
    }
    matches += ac == bc;
  }
  if constexpr (debug) {
    std::cerr << matches << " matches\n"
              << a.to_string() << "\n"
              << b.to_string() << std::endl;
  }
  if (matches >= a_len - 1) {
    return true;
  }
  matches = 0;
  uint64_t a_offset = 1 + (a_gs == 1) * a_gl;
  for (uint64_t i = 0; i < b_len; ++i) {
    uint64_t ii = i + (i >= b_gs) * b_gl;
    auto ac = a.get_c(ii + a_offset);
    auto bc = b.get_c(ii);
    if constexpr (debug) {
      std::cerr << i << ", " << ii << ": " << ac << ", " << bc << std::endl;
    }
    matches += ac == bc;
  }
  if constexpr (debug) {
    std::cerr << matches << " matches vs. " << a_len << "\n"
              << a.to_string() << "\n";
    for (uint64_t i = 0; i < a_offset; ++i) {
      std::cerr << " ";
    }
    std::cerr << b.to_string() << std::endl;
  }
  if (matches >= a_len - 1) {
    return true;
  }
  if (b_len == a_len) {
    // or len - 1 bases align some other way.
    --a_offset;
    while (a_offset) {
      matches = 0;
      for (uint64_t i = 0; i < b_len; ++i) {
        uint64_t ii = i + (i >= b_gs) * b_gl;
        auto ac = a.get_c(ii + a_offset);
        auto bc = b.get_c(ii);
        if constexpr (debug) {
          std::cerr << i << ", " << ii << ": " << ac << ", " << bc << std::endl;
        }
        matches += ac == bc;
      }
      if constexpr (debug) {
        std::cerr << matches << " matches\n" << a.to_string() << "\n";
        for (uint64_t i = 0; i < a_offset; ++i) {
          std::cerr << " ";
        }
        std::cerr << b.to_string() << std::endl;
      }
      if (matches >= a_len - 1) {
        return true;
      }
      --a_offset;
    }
    uint64_t b_offset = 1 + (b_gs == 1) * b_gl;
    while (b_offset) {
      matches = 0;
      for (uint64_t i = 0; i < a_len; ++i) {
        uint64_t ii = i + (i >= a_gs) * a_gl;
        auto ac = a.get_c(ii);
        auto bc = b.get_c(ii + b_offset);
        if constexpr (debug) {
          std::cerr << i << ", " << ii << ": " << ac << ", " << bc << std::endl;
        }
        matches += ac == bc;
      }
      if constexpr (debug) {
        std::cerr << matches << " matches\n";
        for (uint64_t i = 0; i < b_offset; ++i) {
          std::cerr << " ";
        }
        std::cerr << a.to_string() << "\n" << b.to_string() << std::endl;
      }
      if (matches >= a_len - 1) {
        return true;
      }
      --b_offset;
    }
  }
  if constexpr (compare_rc) {
    return is_neighbour<false, debug>(other.reverse_complement());
  }
  if constexpr (debug) {
    std::cerr << "Default false" << std::endl;
  }
  return false;
}

template <bool middle_gap_only, uint16_t t_max_gap>
template <bool compare_rc, bool debug>
bool gapmer<middle_gap_only, t_max_gap>::aligns_to(gapmer other) const {
  if constexpr (debug) {
    std::cerr << "Aligns_to comparison " << to_string() << " to\n"
              << other.to_string() << std::endl;
  }
  uint16_t len = length();
  uint16_t o_len = other.length();
  if (len == o_len) {
    if (other.data_ == data_) {
      return true;
    }
    if constexpr (compare_rc) {
      return aligns_to<false, debug>(other.reverse_complement());
    }
    return false;
  }

  if (len > o_len) {
    return false;
  }
  uint16_t gl = gap_length();
  uint16_t gs = gap_start();
  uint16_t o_gl = other.gap_length();
  uint16_t o_gs = other.gap_start();
  if (gl) {
    // Gap in *this
    uint16_t v_len = len + gl;
    if (o_gl) {
      // gap also on other
      uint16_t ov_len = o_len + o_gl;
      if (v_len > ov_len) {
        return false;
      }
      for (uint16_t offset = 0; offset <= ov_len - v_len; ++offset) {
        bool ok = true;
        if constexpr (debug) {
          for (uint16_t sc = 0; sc < offset; ++sc) {
            std::cerr << " ";
          }
          std::cerr << to_string() << std::endl;
        }
        for (uint16_t i = 0; i < len; ++i) {
          uint16_t o_i = i + offset;
          o_i += i >= gs ? gl : 0;
          if (o_i >= o_gs) {
            if (o_i < o_gs + o_gl) {
              ok = false;
              break;
            }
            o_i -= o_gl;
          }
          if (nuc(i) != other.nuc(o_i)) {
            ok = false;
            break;
          }
        }
        if (ok) {
          return true;
        }
      }
    } else {
      // no gap in other
      if (v_len > o_len) {
        return false;
      }
      for (uint16_t offset = 0; offset <= o_len - v_len; ++offset) {
        bool ok = true;
        if constexpr (debug) {
          for (uint16_t sc = 0; sc < offset; ++sc) {
            std::cerr << " ";
          }
          std::cerr << to_string() << std::endl;
        }
        for (uint16_t i = 0; i < len; ++i) {
          uint16_t o_i = i + offset;
          o_i += i >= gs ? gl : 0;
          if (o_i >= o_gs) {
            if (nuc(i) != other.nuc(o_i)) {
              ok = false;
              break;
            }
          }
          if (ok) {
            return true;
          }
        }
      }
      if constexpr (compare_rc) {
        return aligns_to<false, debug>(other.reverse_complement());
      }
      return false;
    }
  }
  // else no gaps so brute force align is easy.
  // len bases of this must align with bases of other.
  // The first base of this matches one of the first o_len - len bases of
  // other.
  for (uint16_t offset = 0; offset <= o_len - len; ++offset) {
    bool matches = true;
    if constexpr (debug) {
      for (uint16_t space = 0; space < offset; ++space) {
        std::cerr << " ";
      }
      std::cerr << to_string();
    }
    for (uint16_t i = 0; i < len; ++i) {
      if (nuc(i) != other.nuc(i + offset)) {
        matches = false;
        break;
      }
    }
    if constexpr (debug) {
      std::cerr << "\t" << matches << std::endl;
    }
    if (matches) {
      return true;
    }
  }
  if constexpr (compare_rc) {
    return reverse_complement().template aligns_to<false, debug>(other);
  }
  return false;
}

template <bool middle_gap_only, uint16_t t_max_gap>
uint16_t gapmer<middle_gap_only, t_max_gap>::gap_start() const {
  uint64_t ret = data_ >> (max_k * 2 + 5);
  return ret & meta_mask;
}

template <bool middle_gap_only, uint16_t t_max_gap>
uint8_t gapmer<middle_gap_only, t_max_gap>::nuc(uint8_t i) const {
  uint64_t v = data_ >> (length() - i - 1) * 2;
  return v & 0b11;
}

template <bool middle_gap_only, uint16_t t_max_gap>
uint8_t gapmer<middle_gap_only, t_max_gap>::get_c(uint64_t i) const {
  uint8_t gs = gap_start();
  uint8_t gl = gap_length();
  if (i >= gs && i < gs + gl) {
    return '.';
  }
  i -= (i >= gs) * gl;
  return i >= length() ? '.' : v_to_nuc[nuc(i)];
}

template <bool middle_gap_only, uint16_t t_max_gap>
auto gapmer<middle_gap_only, t_max_gap>::next(char c) const -> gapmer {
#ifdef DEBUG
  assert(gap_length() == 0);
#endif
  uint64_t v = data_ << 2;
  v |= nuc_to_v[c];
  v &= (ONE << (length() * 2)) - 1;
  return gapmer{(data_ & ~value_mask) | v};
}

template <bool middle_gap_only, uint16_t t_max_gap>
auto gapmer<middle_gap_only, t_max_gap>::next(char c1, char c2) const
    -> gapmer {
#ifdef DEBUG
  assert(gap_length() > 0);
#endif
  uint64_t suf_len = (length() - gap_start()) * 2;
  uint64_t gap_mask = uint64_t(0b11) << suf_len;
  uint64_t v = (data_ << 2) & ~gap_mask;
  v |= uint64_t(nuc_to_v[c1]) << suf_len;
  v |= nuc_to_v[c2];
  v &= (ONE << (length() * 2)) - 1;
  return gapmer{(data_ & ~value_mask) | v};
}

template <bool middle_gap_only, uint16_t t_max_gap>
void gapmer<middle_gap_only, t_max_gap>::hamming_neighbours(
    auto& callback) const {
  uint16_t len = length();
  for (uint64_t n = 1; n < 4; ++n) {
    uint64_t xor_val = n;
    for (size_t i = 0; i < len; ++i) {
      callback(gapmer(data_ ^ xor_val));
      xor_val <<= 2;
    }
  }
}

template <bool middle_gap_only, uint16_t t_max_gap>
template <bool no_smaller, bool no_same, bool no_larger>
void gapmer<middle_gap_only, t_max_gap>::huddinge_neighbours(
    auto&& callback) const {
  if constexpr (middle_gap_only) {
    middle_gap_neighbours<no_smaller, no_same, no_larger>(callback);
  } else {
    all_gap_neighbours<no_smaller, no_same, no_larger>(callback);
  }
}

template <bool middle_gap_only, uint16_t t_max_gap>
std::string gapmer<middle_gap_only, t_max_gap>::to_string() const {
  std::string ret;
  uint16_t i;
  for (i = 0; i < gap_start(); ++i) {
    uint16_t v = nuc(i);
    auto c = v_to_nuc[v];
    ret.push_back(c);
  }
  for (uint16_t j = 0; j < gap_length(); ++j) {
    ret.push_back('.');
  }
  for (; i < length(); ++i) {
    uint16_t v = nuc(i);
    auto c = v_to_nuc[v];
    ret.push_back(c);
  }
  return ret;
}

template <bool middle_gap_only, uint16_t t_max_gap>
bool gapmer<middle_gap_only, t_max_gap>::is_canonical() const {
  uint8_t a = 0;
  uint8_t b = length() - 1;
  while (a <= b) {
    auto a_n = nuc(a);
    auto b_n = (~nuc(b)) & 0b11;
    if (a_n != b_n) {
      return a_n < b_n;
    }
    ++a;
    --b;
  }
  return true;
}

template <bool middle_gap_only, uint16_t t_max_gap>
auto gapmer<middle_gap_only, t_max_gap>::reverse_complement() const -> gapmer {
  uint64_t v = value();
  uint64_t n_v = 0;
  uint8_t l = length();
  for (uint16_t i = 4; i <= max_k; i += 4) {
    n_v <<= 8;
    n_v |= rc_byte[uint8_t(v)];
    v >>= 8;
    if (i >= l) {
      break;
    }
  }
  uint8_t lm = l % 4;
  if (lm) {
    n_v >>= 2 * (4 - lm);
  }
  uint8_t g_s = gap_start();
  if (g_s == 0) {
    return {n_v, l};
  }
  g_s = l - g_s;
  uint8_t g_l = gap_length();
  return {n_v, l, g_s, g_l};
}

template <bool middle_gap_only, uint16_t t_max_gap>
bool gapmer<middle_gap_only, t_max_gap>::is_valid() const {
  uint8_t len = length();
  uint64_t val = suffix(len);
  uint8_t g_s = gap_start();
  uint8_t g_l = gap_length();
  if (g_s == 0 && g_l > 0) {
    std::cerr << "Gapmer validation error gap length != 0: " << int(g_l)
              << std::endl;
    return false;
  }
  if (g_l == 0 && g_s > 0) {
    std::cerr << "Gapmer validation error gap start != 0: " << int(g_s)
              << std::endl;
    return false;
  }
  if (len && g_s >= len) {  // g_s == 0 && g_l > 0 checked earlier.
    std::cerr << "Gapmer validation error gap start >= " << int(len) << ": "
              << int(g_s) << std::endl;
    return false;
  }
  if (g_l > max_gap) {
    std::cerr << "Gapmer validation error gap length >= " << int(max_gap)
              << ": " << int(g_l) << std::endl;
    return false;
  }
  gapmer o(val, len, g_s, g_l);
  bool valid = data_ == o.data_;
  if (not valid) {
    std::cerr << "Gapmer validation error:\n"
              << to_string() << ", " << bits() << "\n"
              << o.to_string() << ", " << o.bits() << "\nShould be equal"
              << std::endl;
  }
  return valid;
}

template <bool middle_gap_only, uint16_t t_max_gap>
template <std::size_t t_n>
void gapmer<middle_gap_only, t_max_gap>::write_2bit_coded_to_buffer(
    std::span<uint64_t, t_n> span) const
  requires((2U * (max_k + max_gap) + 7U) / 8U <= t_n)
{
  auto const ll(length());
  auto const gs(gap_start());
  auto const gl(gap_length());
  span.front() = prefix();
  bits::shift_left(span, 2U * (ll - gs + gl));  // Prefix length + gap length
  span.front() |= suffix();
}

template <bool middle_gap_only, uint16_t t_max_gap>
uint16_t gapmer<middle_gap_only, t_max_gap>::huddinge_distance(
    gapmer const other, int& out) const {
  out = 0;
  // Given this's and other's (gapped) lengths l_1 and l_2, we iterate over the
  // l_1 + l_2 - 1 alignments, check for matching characters with XOR and count
  // them with PEXT and POPCOUNT. This could possibly be made faster by using
  // large registers instead of arrays and applying SIMD.
  auto const llen{length()};
  auto const rlen{other.length()};
  if (!(llen && rlen)) return std::max(llen, rlen);

  // Set up the defined character masks.
  auto const lplen{gap_start()};
  auto const rplen{other.gap_start()};
  auto const lglen{gap_length()};
  auto const rglen{other.gap_length()};
  auto const lslen{llen - lplen};
  auto const rslen{rlen - rplen};
  auto const lmask{~(~(~(uint64_t(-1) << lplen) << lglen) << lslen)};
  auto const rmask{~(~(~(uint64_t(-1) << rplen) << rglen) << rslen)};

  // Set up the buffers.
  constexpr auto const max_aln_length(2 * (max_k + max_gap) - 1U);
  uint64_t buf1[(2 * max_aln_length + 7) / 8]{};
  uint64_t buf2[(2 * max_aln_length + 7) / 8]{};
  uint64_t maskbuf2[(max_aln_length + 7) / 8]{};
  std::span b1s{buf1};
  std::span b2s{buf2};
  std::span m2s{maskbuf2};

  // The longer one goes to buf1.
  auto const m1([&]() {
    if (llen + lglen <= rlen + rglen) {
      write_2bit_coded_to_buffer(b2s);
      other.write_2bit_coded_to_buffer(b1s);
      m2s.front() = lmask;
      bits::shift_left(b2s, 2 * (rlen + rglen - 1U));
      bits::shift_left(m2s, rlen + rglen - 1U);
      return rmask;
    } else {
      write_2bit_coded_to_buffer(b1s);
      other.write_2bit_coded_to_buffer(b2s);
      m2s.front() = rmask;
      bits::shift_left(b2s, 2 * (llen + lglen - 1U));
      bits::shift_left(m2s, llen + lglen - 1U);
      return lmask;
    }
  }());

  constexpr auto const pm1{UINT64_C(0x5555'5555'5555'5555)};  // PEXT mask 1
  constexpr auto const pm2{UINT64_C(0xAAAA'AAAA'AAAA'AAAA)};  // PEXT mask 2
  std::size_t const count{llen + lglen + rlen + rglen - 1U};
  uint16_t max_score{};
  for (std::size_t i{}; i < count; ++i) {
    // Compare.
    auto const cm{b1s.front() ^ b2s.front()};
    auto const p1{bits::pext(cm, pm1)};
    auto const p2{bits::pext(cm, pm2)};
    auto const cm_{p1 | p2};  // Get the non-matching character mask.

    // Get the relevant (aligned) characters.
    auto const mm(m1 & m2s.front());
    auto const cm__((~cm_) & mm);

    // Update the score if needed.
    uint16_t const score(std::popcount(cm__));
    if (score > max_score) {
      max_score = score;
      out = i;
    }

    // Shift the buffers.
    bits::shift_right(b2s, 2U);
    bits::shift_right(m2s, 1U);
  }
  if (llen + lglen <= rlen + rglen) {
    out -= llen + lglen - 1;
    out *= -1;
  } else {
    out -= rlen + rglen - 1;
  }

  return std::max(llen, rlen) - max_score;
}

}  // namespace sf
