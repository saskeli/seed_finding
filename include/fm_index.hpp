#include <algorithm>
#include <array>
#include <cstdint>
#include <SeqIO/SeqIO.hh>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/hyb_vector.hpp>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
#include <cstdint>

#ifdef DEBUG
#include <unordered_map>
#endif

#include "gapmer.hpp"

namespace sf {

template <uint16_t SA_gap_size = 64, uint16_t block_size = 2048>
class fm_index {
 private:
  static const constexpr uint64_t max_run_length = 0b111111;

  std::vector<uint64_t> sample_locations_;
  std::vector<std::array<uint64_t, 4>> partial_sums_;
  uint64_t* bwt_;
  std::array<uint64_t, 5> C_ = {0, 0, 0, 0, 0};
  sdsl::hyb_vector<> samples_;
  sdsl::hyb_vector<>::rank_1_type samples_rs_;
  sdsl::hyb_vector<> seq_starts_;
  sdsl::hyb_vector<>::rank_1_type seq_starts_rs_;

  uint8_t at(uint64_t i) const {
    uint64_t w = bwt_[i / 32];
    w >>= (i % 32) * 2;
    return w & 0b11;
  }

  uint64_t rank(uint64_t i, uint8_t v) const {
    uint64_t block = i / block_size;
    uint64_t count = partial_sums_[block][v];
    for (uint64_t idx = block * block_size; idx < i; ++idx) {
      count += at(idx) == v;
    }
    return count;
  }

  uint64_t lf(const uint64_t& i) const {
    uint8_t c = at(i);
    return C_[c] + rank(i, c);
  }

  std::pair<uint64_t, uint64_t> find(const std::string& s) const {
    uint64_t i = s.size() - 1;
    uint8_t c = nuc_to_v[s[i--]];
    uint64_t a = C_[c];
    uint64_t b = C_[c + 1];
    for (; i < s.size(); --i) {
      c = nuc_to_v[s[i]];
      a = C_[c] + rank(a, c);
      b = C_[c] + rank(b, c);
    }
    return {a, b};
  }

  template <class gapmer>
  std::pair<uint64_t, uint64_t> find(const gapmer& mer, uint16_t start,
                                     uint16_t end) const {
    uint64_t i = end - 1;
    uint8_t c = mer.nuc(i--);
    uint64_t a = C_[c];
    uint64_t b = C_[c + 1];
    for (; i < end; --i) {
      if (i < start) {
        break;
      }
      c = mer.nuc(i);
      a = C_[c] + rank(a, c);
      b = C_[c] + rank(b, c);
    }
    return {a, b};
  }

  fm_index(fm_index&) = delete;
  fm_index(fm_index&&) = delete;
  fm_index& operator=(fm_index&) = delete;
  fm_index& operator=(fm_index&&) = delete;

 public:
  ~fm_index() { free(bwt_); }

  fm_index(const char* fasta_path) : sample_locations_(), partial_sums_() {
    std::vector<uint64_t> starts;
    seq_io::Reader r(fasta_path);
    r.enable_reverse_complements();
    string_buffer <uint64_t> seq;
    while (true) {
      uint64_t len = r.get_next_read_to_buffer();
      if (len == 0) {
        break;
      }
      starts.push_back(seq.size);
      seq.append(std::string_view{r.read_buf, len});
    }
    C_[4] = seq.size;
    std::cerr << "Creating index from " << fasta_path << " with " << C_[4]
              << " total characters" << std::endl;
    sdsl::bit_vector samples(C_[4]);
    sdsl::bit_vector bwt_starts(C_[4]);
    sdsl::bit_vector seq_starts(C_[4]);
    for (auto v : starts) {
      seq_starts[v] = 1;
    }
    seq_starts_ = sdsl::hyb_vector<>(seq_starts);
    seq_starts_rs_ = sdsl::hyb_vector<>::rank_1_type(&seq_starts_);

    bwt_ = (uint64_t*)calloc((2 * C_[4] + 63) / 64, sizeof(uint64_t));
    uint64_t bwt_position = 0;
    seq.append(seq.c_str(), 16);
    std::array<uint64_t, 4> partial{};
    std::array<std::vector<uint64_t>, 256> buckets{};
    for (uint64_t c_i = 0; c_i < 4; ++c_i) {
      if (c_i > 0) {
        C_[c_i] = C_[c_i - 1];
      }
      char nuc = v_to_nuc[c_i];

      for (uint64_t i = 0; i < 256; ++i) {
        buckets[i].clear();
      }
      for (uint64_t i = 0; i < C_[4]; ++i) {
        if (seq[i] == nuc) {
          ++C_[c_i];
          sf::gapmer km(seq.c_str() + i, 5);
          uint64_t v = km.value() % 256;
          buckets[v].push_back(i);
        }
      }
#pragma omp parallel for
      for (uint16_t i = 0; i < 256; ++i) {
        std::sort(buckets[i].begin(), buckets[i].end(),
                  [&](uint64_t lhs, uint64_t rhs) {
                    uint64_t a_comp = lhs + 5;
                    uint64_t b_comp = rhs + 5;
                    while (true) {
                      a_comp -= a_comp >= C_[4] ? C_[4] : 0;
                      b_comp -= b_comp >= C_[4] ? C_[4] : 0;
                      if (seq[a_comp] == seq[b_comp]) {
                        ++a_comp;
                        ++b_comp;
                      } else {
                        return seq[a_comp] < seq[b_comp];
                      }
                    }
                  });
      }
      for (uint16_t i = 0; i < 256; ++i) {
        for (auto v : buckets[i]) {
          if (bwt_position % block_size == 0) [[unlikely]] {
            partial_sums_.push_back(partial);
          }
          uint64_t bwt_i = v ? v - 1 : C_[4] - 1;
          if (bwt_i % SA_gap_size == 0) [[unlikely]] {
            sample_locations_.push_back(bwt_i);
            samples[bwt_position] = 1;
          }
          uint64_t c_v = nuc_to_v[seq[bwt_i]];
          ++partial[c_v];
          c_v <<= ((bwt_position % 32) * 2);
          bwt_[bwt_position / 32] |= c_v;
          ++bwt_position;
        }
      }
    }
    if (bwt_position % block_size == 0) [[unlikely]] {
      partial_sums_.push_back(partial);
    }

#ifdef DEBUG
    assert(C_[3] == C_[4]);
#endif
    for (uint16_t i = 4; i > 0; --i) {
      C_[i] = C_[i - 1];
    }
    C_[0] = 0;

    samples_ = sdsl::hyb_vector<>(samples);
    samples_rs_ = sdsl::hyb_vector<>::rank_1_type(&samples_);
  }

  template <class gapmer>
  uint64_t count(const gapmer& mer, uint16_t gap_len) const {
    // TODO: Would be better to compute a locations first.
    // can easily check if the match is impossible due to overlapping
    // a read end if we start with the prefix.
    auto a_I = find(mer, 0, mer.gap_start());
    auto b_I = find(mer, mer.gap_start(), mer.length());
    uint32_t gap_e = gap_len + mer.gap_start();

    // get all the suffix start locations in vector
    std::vector<uint64_t> b_loc;
    for (uint64_t s_idx = b_I.first; s_idx < b_I.second; ++s_idx) {
      uint64_t idx = s_idx;
      uint64_t steps = 0;
      for (; steps < SA_gap_size; ++steps) {
        if (samples_[idx]) [[unlikely]] {
          break;
        }
        idx = lf(idx);
      }
      uint64_t loc = samples_rs_.rank(idx);
      loc = sample_locations_[loc];
      loc += steps + 1;
      loc -= loc >= C_[4] ? C_[4] : 0;
      if (loc + mer.length() - mer.gap_start() > C_[4]) [[unlikely]] {
        continue;
      }
      if (seq_starts_rs_.rank(loc + 1) ==
          seq_starts_rs_.rank(loc + (mer.length() - mer.gap_start()))) {
        b_loc.push_back(loc);
      }
    }
    std::sort(b_loc.begin(), b_loc.end());

    uint64_t acc = 0;
    // match prefix locations to suffis locations.
    for (uint64_t s_idx = a_I.first; s_idx < a_I.second; ++s_idx) {
      uint64_t idx = s_idx;
      uint64_t steps = 0;
      for (; steps < SA_gap_size; ++steps) {
        if (samples_[idx]) [[unlikely]] {
          break;
        }
        idx = lf(idx);
      }
      uint64_t loc = samples_rs_.rank(idx);
      loc = sample_locations_[loc];
      loc += steps + 1;
      loc -= loc >= C_[4] ? C_[4] : 0;
      if (loc + mer.gap_start() > C_[4]) [[unlikely]] {
        continue;
      }
      if (seq_starts_rs_.rank(loc + 1) ==
          seq_starts_rs_.rank(loc + mer.gap_start())) {
        if (std::binary_search(b_loc.begin(), b_loc.end(), loc + gap_e)) {
          if (seq_starts_rs_.rank(loc + 1) ==
              seq_starts_rs_.rank(loc + gap_e + 1)) {
            ++acc;
          }
        }
      }
    }
    return acc;
  }

  template <class gapmer>
  uint64_t count(const gapmer& mer) const {
    auto gap_len = mer.gap_length();
    if (gap_len) {
      return count(mer, gap_len);
    }
    auto mer_len = mer.length();
    auto I = find(mer, 0, mer_len);
    uint64_t acc = 0;
    for (uint64_t s_idx = I.first; s_idx < I.second; ++s_idx) {
      uint64_t idx = s_idx;
      uint64_t steps = 0;
      for (; steps < SA_gap_size; ++steps) {
        if (samples_[idx]) [[unlikely]] {
          break;
        }
        idx = lf(idx);
      }
      uint64_t loc = samples_rs_.rank(idx);
      loc = sample_locations_[loc];
      loc += steps + 1;
      loc -= loc >= C_[4] ? C_[4] : 0;
      if (loc + mer.length() > C_[4]) [[unlikely]] {
        continue;
      }
      acc += seq_starts_rs_.rank(loc + 1) ==
             seq_starts_rs_.rank(loc + mer.length());
    }
    return acc;
  }

  std::vector<uint64_t> locate(const std::string& s) const {
    auto gap_s = s.find(".");
    std::vector<uint64_t> ret;
    if (gap_s != std::string::npos) {
      uint64_t gap_e = gap_s + 1;
      for (; gap_e < s.size(); ++gap_e) {
        if (s[gap_e] != '.') {
          break;
        }
      }
      if (gap_s == 0) [[unlikely]] {
        if (gap_e == s.size()) {
          return ret;
        }
        return locate(s.substr(gap_e));
      } else if (gap_e == s.size()) {
        return locate(s.substr(0, gap_s));
      }
      auto a_locs = locate(s.substr(0, gap_s));
      auto b_locs = locate(s.substr(gap_e));
      std::sort(b_locs.begin(), b_locs.end());
      for (auto a_loc : a_locs) {
        if (std::binary_search(b_locs.begin(), b_locs.end(), a_loc + gap_e)) {
          if (seq_starts_rs_.rank(a_loc + 1) ==
              seq_starts_rs_.rank(a_loc + gap_e + 1)) {
            ret.push_back(a_loc);
          }
        }
      }
      return ret;
    }
    auto I = find(s);

    for (uint64_t s_idx = I.first; s_idx < I.second; ++s_idx) {
      uint64_t idx = s_idx;
      uint64_t steps = 0;
      for (; steps < SA_gap_size; ++steps) {
        if (samples_[idx]) [[unlikely]] {
          break;
        }
        idx = lf(idx);
      }
      uint64_t loc = samples_rs_.rank(idx);
      loc = sample_locations_[loc];
      loc += steps + 1;
      loc -= loc >= C_[4] ? C_[4] : 0;
      if (loc + s.size() > C_[4]) [[unlikely]] {
        continue;
      }
      if (seq_starts_rs_.rank(loc + 1) == seq_starts_rs_.rank(loc + s.size())) {
        ret.push_back(loc);
      }
    }
    return ret;
  }

  /// Return the length of the original text.
  uint64_t size() const { return C_[4]; }
};
}  // namespace sf
