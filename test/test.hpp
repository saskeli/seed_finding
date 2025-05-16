/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <exception>			// IWYU pragma: keep // Needed by RapidCheck.

#include <gtest/gtest.h>		// IWYU pragma: export
#include <rapidcheck/gtest.h>	// IWYU pragma: export

// GoogleTest seems to require a fixture for each type-paramterised test,
// so we have a macro to help with that.
#define SF_RC_TEMPLATE_TEST(FIXTURE_PREFIX, TEST, TEST_PARAMS, ...) \
	template <typename> struct FIXTURE_PREFIX##_fixture : public testing::Test{}; \
	typedef ::testing::Types <__VA_ARGS__> FIXTURE_PREFIX##_test_types; \
	TYPED_TEST_SUITE(FIXTURE_PREFIX##_fixture, FIXTURE_PREFIX##_test_types); \
	RC_GTEST_TYPED_FIXTURE_PROP(FIXTURE_PREFIX##_fixture, TEST, TEST_PARAMS)
