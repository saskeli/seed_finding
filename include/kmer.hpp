#include <immintrin.h>

#include <bitset>
#include <cstdint>
#include <string>
#include <type_traits>

namespace sf {
template <uint16_t max_k = 15>
class kmer {
 private:
  typedef std::conditional<15 < max_k, uint64_t, uint32_t>::type uint_t;

  uint_t data_;
  uint16_t size_;

 public:
  kmer(std::string& s, uint16_t k) : data_(0), size_(k) {
    uint64_t iv;
    const uint64_t xor_mask = 0b0101010101010101;
    const uint64_t pext_mask = 0x0606060606060606;
    for (uint16_t i = 0; i < max_k / 8; ++i) {
      iv = __bswap_64(reinterpret_cast<uint64_t*>(s.data())[i]);
      iv = _pext_u64(iv, pext_mask);
      iv ^= (iv >> 1) & xor_mask;
      data_ |= iv;
      if (k <= (1 + i) * 8) {
        break;
      }
      data_ <<= 16;
    }
    uint16_t shr = k % 8;
    data_ >>= shr ? (16 - 2 * shr) : 0;
  }

  uint_t value() const { return data_; }

  void print() const {
    std::cout << size_ << ", " << std::bitset<15 < max_k ? 64 : 32>(data_)
              << std::endl;
  }
};
}  // namespace sf
