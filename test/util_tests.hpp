#pragma once

namespace sf {
TEST(Util, ByteRC1) { ASSERT_EQ(rc_byte[0b00011011], 0b00011011); }

TEST(Util, ByteRC2) { ASSERT_EQ(rc_byte[0b11010111], 0b00101000); }

TEST(Util, ByteRC3) { ASSERT_EQ(rc_byte[0b11111111], 0b00000000); }

TEST(Util, ByteRC4) { ASSERT_EQ(rc_byte[0b00000000], 0b11111111); }
}  // namespace sf