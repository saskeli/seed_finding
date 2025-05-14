/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <exception> // Needed by RapidCheck.

#include <cstdint>
#include <gtest/gtest.h>
#include <rapidcheck/gtest.h>
#include "../include/bits.hpp"

namespace sf {

	// Needed by the template test below.
	template <typename> struct bit_arbitrary_pext_fixture : public testing::Test {};

	typedef ::testing::Types <uint32_t, uint64_t> bit_arbitrary_pext_test_types;
	TYPED_TEST_SUITE(bit_arbitrary_pext_fixture, bit_arbitrary_pext_test_types);

	RC_GTEST_TYPED_FIXTURE_PROP(bit_arbitrary_pext_fixture, pextWorkAsExpected, (TypeParam const value, TypeParam const mask)) {
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
