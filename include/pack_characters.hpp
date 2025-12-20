#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>


namespace sf {
typedef std::vector<std::uint64_t> packed_word_vector;

// The non-lenient version returns UINT64_MAX if an unexpected character is encountered.
[[nodiscard]] std::uint64_t pack_characters(std::string_view sv,
                                            std::vector<std::uint64_t>& dst,
                                            std::uint64_t dst_pos = 0);
std::uint64_t pack_characters_lenient(std::string_view sv,
                                      std::vector<std::uint64_t>& dst,
                                      std::uint64_t dst_pos = 0);
void unpack_characters(std::vector<std::uint64_t> const& src,
                       std::uint64_t length, std::string& dst);

void print_read(std::span<std::uint64_t const> span, std::uint64_t len);
}  // namespace sf
