#include <immintrin.h>

#include <array>
#include <bitset>
#include <cstdint>
#include <string>
#include <type_traits>
#ifdef DEBUG
#include <cassert>
#endif

namespace sf {
class gapmer {
 private:
  const static constexpr uint64_t max_k = 24;
  const static constexpr uint64_t value_mask = (uint64_t(1) << (max_k * 2)) - 1;
  const static constexpr uint64_t xor_mask = 0b0101010101010101;
  const static constexpr uint64_t pext_mask = 0x0606060606060606;
  const static constexpr uint64_t meta_mask = 0x11111;
  const static constexpr std::array<char, 4> nucs = {'A', 'C', 'G', 'T'};
  class huddingeIterator;

  uint64_t data_;

 public:
  gapmer(std::string& s, uint8_t k) : data_(0) {
#ifdef DEBUG
    assert(k > 0);
    assert(k <= 24);
#endif
    uint64_t iv;
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
    data_ |= uint64_t(k) << (max_k * 2);
  }

  gapmer(std::string& s, uint8_t k, uint8_t gap_start, uint8_t gap_length) : data_(0) {
#ifdef DEBUG
    assert(k > 0);
    assert(k <= 24);
    assert(gap_start < 24);
    assert(gap_length > 0);
    assert(gap_length < 32);
#endif
    uint64_t iv;
    uint64_t* d_ptr = reinterpret_cast<uint64_t*>(s.data());
    for (uint16_t i = 0; i < max_k / 8; ++i) {
      iv = __bswap_64(d_ptr[i]);
      iv = _pext_u64(iv, pext_mask);
      iv ^= (iv >> 1) & xor_mask;
      data_ |= iv;
      if (gap_start <= (1 + i) * 8) {
        break;
      }
      data_ <<= 16;
    }
    uint16_t shr = gap_start % 8;
    data_ >>= shr ? (16 - 2 * shr) : 0;

    uint8_t k_suf = k - gap_start;
    d_ptr = reinterpret_cast<uint64_t*>(s.data() + gap_start + gap_length);
    for (uint16_t i = 0; i < max_k / 8; ++i) {
      data_ <<= 16;
      iv = __bswap_64(d_ptr[i]);
      iv = _pext_u64(iv, pext_mask);
      iv ^= (iv >> 1) & xor_mask;
      data_ |= iv;
      if (k_suf <= (1 + i) * 8) {
        break;
      }
    }
    shr = k_suf % 8;
    data_ >>= shr ? (16 - 2 * shr) : 0;
    data_ |= uint64_t(k) << (max_k * 2);
    data_ |= uint64_t(gap_start) << (max_k * 2 + 5);
    data_ |= uint64_t(gap_length) << (max_k * 2 + 10);
  }

  uint8_t length() const { return (data_ >> (max_k * 2)) & meta_mask; }

  uint8_t gap_start() const { return (data_ >> (max_k * 2 + 5)) & meta_mask; }

  uint8_t gap_length() const { return (data_ >> (max_k * 2 + 10)) & meta_mask; }

  uint8_t nuc(uint8_t i) const {
    uint64_t v = data_ >> (length() - i - 1) * 2;
    return v & 0x11;
  }

  uint64_t value() const { return data_ & value_mask; }

  std::string to_string() const {
    std::string ret;
    uint16_t i;
    for (i = 0; i < gap_start(); ++i) {
      ret.push_back(nucs[nuc(i)]);
    }
    for (uint16_t j = 0; j < gap_length(); ++j) {
      ret.push_back('n');
    }
    for (; i < length(); ++i) {
      ret.push_back(nucs[nuc(i)]);
    }
    return ret;
  }
};

class gapmer::huddingeIterator {};
}  // namespace sf
