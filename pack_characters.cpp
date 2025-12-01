#include <bit>
#include <cstdint>
#include <iostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>


namespace sf {

void print_read(std::span<std::uint64_t const> span, std::uint64_t len) {
  std::cerr << "** read: ";
  for (auto word : span) {
    for (std::uint8_t ii{}; ii < 32U; ++ii) {
      word = std::rotl(word, 2);
      auto const cc{word & 0x3};
      switch (cc) {
        case 0:
          std::cerr << 'A';
          break;
        case 1:
          std::cerr << 'C';
          break;
        case 2:
          std::cerr << 'G';
          break;
        case 3:
          std::cerr << 'T';
          break;
        default:
          throw std::runtime_error("Unexpected character");
      }

      --len;
      if (0 == len) break;
    }
  }
  std::cerr << '\n';
}


std::uint64_t pack_characters(std::string_view sv,
                              std::vector<std::uint64_t>& dst,
                              std::uint64_t dst_pos) {
  auto const push_packed([](char const cc, std::uint64_t& dst_word,
                            std::uint64_t const dst_pos_) -> bool {
    auto const chararcter_pos{dst_pos_ % 32U};
    std::uint64_t const shift_amt{62U - 2 * chararcter_pos};
    switch (cc) {
      case 'A':
      case 'a':
        break;

      case 'C':
      case 'c':
        dst_word |= UINT64_C(0x1) << shift_amt;
        break;

      case 'G':
      case 'g':
        dst_word |= UINT64_C(0x2) << shift_amt;
        break;

      case 'T':
      case 't':
        dst_word |= UINT64_C(0x3) << shift_amt;
        break;

      [[unlikely]] default:
        return false;
    }

    return true;
  });

  auto it{sv.cbegin()};
  auto const end{sv.cend()};

                // If the current position is not divisible by 32, dst
                // must be non-empty.
  if (dst_pos % 32U) {
    auto& dst_word{dst.back()};
    do {
      if (it == end) return dst_pos;

      if (!push_packed(*it, dst_word, dst_pos)) [[unlikely]]
        return UINT64_MAX;

      ++dst_pos;
      ++it;
    } while (dst_pos % 32U);
  }

  if (it == end) return dst_pos;

                // 0 == dst_pos % 32U.
  while (true) {
    auto& dst_word{dst.emplace_back(0)};
    do {
      if (!push_packed(*it, dst_word, dst_pos)) [[unlikely]]
        return UINT64_MAX;

      ++dst_pos;
      ++it;

      if (it == end) return dst_pos;
    } while (dst_pos % 32U);

    if (it == end) return dst_pos;
  }
}


void unpack_characters(std::vector<std::uint64_t> const& src,
                       std::uint64_t length, std::string& dst) {
  if (src.empty()) return;

  auto const append_character([&dst](std::uint64_t const word) {
    switch (word & UINT64_C(0x3)) {
      case 0x0:
        dst.push_back('A');
        break;
      case 0x1:
        dst.push_back('C');
        break;
      case 0x2:
        dst.push_back('G');
        break;
      case 0x3:
        dst.push_back('T');
        break;
      default:
        throw std::runtime_error("Unexpected value");
    }
  });

  std::uint64_t ii{};
  while (ii + 32 <= length) {
    auto word{src[ii / 32]};
    for (std::uint8_t jj{}; jj < 32; ++jj) {
      word = std::rotl(word, 2);
      append_character(word);
    }

    ii += 32;
  }

  auto word{src[ii / 32]};
  while (ii < length) {
    word = std::rotl(word, 2);
    append_character(word);
    ++ii;
  }
}
}  // namespace sf
