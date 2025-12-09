#pragma once

#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>


namespace sf {
typedef std::vector<std::uint64_t> packed_word_vector;

[[nodiscard]] std::uint64_t pack_characters(std::string_view sv,
                                            std::vector<std::uint64_t>& dst,
                                            std::uint64_t dst_pos);
void unpack_characters(std::vector<std::uint64_t> const& src,
                       std::uint64_t length, std::string& dst);

void print_read(std::span<std::uint64_t const> span, std::uint64_t len);
}  // namespace sf
