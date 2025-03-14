#include <vector>
#include <array>

#include "SeqIO/SeqIO.hh"

#include "gapmer.hpp"

namespace sf {

class fm_index {
 public:
  fm_index(const char* fasta_path) {
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
    uint64_t total_len = seq.size();
    seq.append(seq.c_str(), 4);
    std::array<std::vector<uint64_t>, 256> buckets;
    for (uint64_t i = 0; i < total_len; ++i) {
        if (seq[i] == 'A') {
            sf::gapmer km(seq.c_str() + i, 5);
            buckets[km.value()].push_back(i);
        }
    }
    
  }
};
} // namespace sf
