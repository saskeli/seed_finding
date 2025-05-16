/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <cstdint>
#include "test.hpp"
#include "../include/bits.hpp"

namespace sf {

	SF_RC_TEMPLATE_TEST(bit_arbitrary_pext, PEXTWorksAsExpected, (TypeParam const value, TypeParam const mask), uint32_t, uint64_t) {
		if constexpr (sf::bits::detail::pext_intrinsic_available)
		{
			auto const res(sf::bits::detail::pext(value, mask));
			auto const expected(sf::bits::detail::pext_intrinsic(value, mask));
			RC_ASSERT(res == expected);
		}
		else
		{
			GTEST_SKIP() << "Test skipped: Unable to test the PEXT implementation when the intrinsic is not available.";
		}
	}
}
