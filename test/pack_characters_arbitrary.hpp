#include <cstdint>
#include <range/v3/view/all.hpp>
#include <range/v3/view/for_each.hpp>
#include <range/v3/view/zip.hpp>
#include <string>
#include <vector>

#include "../include/libbio_reader_adapter.hpp"
#include "gtest/gtest.h"
#include "nucleotide.hpp"
#include "test.hpp"


RC_GTEST_PROP(pack_characters_arbitrary, PackAndUnpack,
              (std::vector<std::vector<nucleotide>> const& text)) {
  {
    auto const total_length{[&] {
      std::uint64_t retval{};
      for (auto const& vec : text) retval += vec.size();
      return retval;
    }()};
    RC_TAG(text.size(), total_length);
  }

  std::vector<std::uint64_t> packed;
  std::string unpacked;

  std::string buffer;
  std::uint64_t pos{};

        // Pack the characters.
  for (auto const& vec : text) {
    buffer.clear();
    for (auto const nt : vec) buffer.push_back(char(nt));

    pos = sf::pack_characters(buffer, packed, pos);
    EXPECT_NE(UINT64_MAX, pos);
  }

        // Unpack.
  sf::unpack_characters(packed, pos, unpacked);

        // Flatten the range.
  auto rng{ranges::views::for_each(
      text, [](auto const& vec) { return ranges::views::all(vec); })};

  for (auto const& [lhs, rhs] : ranges::views::zip(rng, unpacked))
    EXPECT_EQ(lhs, rhs);
}
