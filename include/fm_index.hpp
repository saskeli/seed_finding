#include <algorithm>
#include <array>
#include <vector>

#include "SeqIO/SeqIO.hh"
#include "gapmer.hpp"

#include "sdsl/bit_vectors.hpp"
#include "sdsl/hyb_vector.hpp"

namespace sf {

template <uint16_t SA_gap_size = 64, uint64_t block_size = 2048> 
class fm_index {
 private:
  std::vector<uint64_t> sample_locations_;
  std::array<uint64_t, 5> C_ = {0, 0, 0, 0, 0};
  sdsl::hyb_vector<> samples_;
  sdsl::hyb_vector<> bwt_starts_;

 public:
  fm_index(const char* fasta_path): sample_locations_() {
    std::vector<uint64_t> starts;
    seq_io::Reader r(fasta_path);
    std::string seq;
    while (true) {
      uint64_t len = r.get_next_read_to_buffer();
      if (len == 0) {
        break;
      }
      starts.push_back(seq.size());
      seq.append(r.read_buf, len);
    }
    C_[4] = seq.size();
    sdsl::bit_vector samples(C_[4]);
    sdsl::bit_vector bwt_starts(C_[4]);
    sdsl::bit_vector seq_starts(C_[4]);
    for (auto v : starts) {
        seq_starts[v] = 1;
    }

    std::string bwt = "";
    seq.append(seq.c_str(), 4);
    for (uint64_t c_i = 0; c_i < 4; ++c_i) {
      if (c_i > 0) {
        C_[c_i] = C_[c_i - 1];
      }
      char nuc = v_to_nuc[c_i];
      std::array<std::vector<uint64_t>, 256> buckets;
      for (uint64_t i = 0; i < C_[4]; ++i) {
        if (seq[i] == nuc) {
          ++C_[c_i];
          sf::gapmer km(seq.c_str() + i, 5);
          buckets[km.value() % 256].push_back(i);
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
          uint64_t bwt_i = v ? v - 1 : C_[4] - 1;
          if (seq_starts[bwt_i]) [[unlikely]] {
            bwt_starts[bwt.size()] = 1;
          }
          if (bwt_i % SA_gap_size == 0) [[unlikely]] {
            sample_locations_.push_back(bwt_i);
            samples[bwt.size()] = 1;
          }
          bwt.push_back(seq[bwt_i]);
        }
      }
    }

    samples_ = sdsl::hyb_vector<>(samples);
    bwt_starts_ = sdsl::hyb_vector<>(bwt_starts);
    //TODO: Make blockRLBWT
    
  }
};
}  // namespace sf
