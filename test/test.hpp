/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <exception>			// IWYU pragma: keep // Needed by RapidCheck.

#include <gtest/gtest.h>		// IWYU pragma: export
#include <rapidcheck/gtest.h>	// IWYU pragma: export

#if defined(SF_REPORT_TESTED_VALUE_DISTRIBUTION) && SF_REPORT_TESTED_VALUE_DISTRIBUTION
#	define SF_RC_TAG(...) RC_TAG(__VA_ARGS__)
#else
#	define SF_RC_TAG(...)
#endif

// GoogleTest seems to require a fixture for each type-paramterised test,
// so we have a macro to help with that.
#define SF_RC_TEMPLATE_TEST(FIXTURE_PREFIX, TEST, TEST_PARAMS, ...) \
	template <typename> struct FIXTURE_PREFIX##_rc_fixture : public testing::Test{}; \
	typedef ::testing::Types <__VA_ARGS__> FIXTURE_PREFIX##_test_types; \
	TYPED_TEST_SUITE(FIXTURE_PREFIX##_rc_fixture, FIXTURE_PREFIX##_test_types); \
	RC_GTEST_TYPED_FIXTURE_PROP(FIXTURE_PREFIX##_rc_fixture, TEST, TEST_PARAMS)

// RapidCheck seems to repeatedly test with empty values in some cases.
// We would like to do so but only once, so we use GoogleTest for the empty case.
// Since we are not in RapidCheck's test case, we would like to conditionally use
// GoogleTest's assertions.
#define SF_ASSERT(CONDITION) (([&](){ if (sf_test_uses_rapidcheck) RC_ASSERT(CONDITION) else ASSERT_TRUE(CONDITION))())

// Mainly for checking if empty parentheses are passed in place of a macro parameter.
// We could detect this with __VA_OPT__ but since it is fairly new, we use the GCC extension
// which causes the comma in , ##__VA_ARGS__ to be deleted if __VA_ARGS__ is empty.
#define SF_NARGS_(_0, _16, _15, _14, _13, _12, _11, _10, _9, _8, _7, _6, _5, _4, _3, _2, _1, N, ...) N
#define SF_NARGS(...) SF_NARGS_(0, ##__VA_ARGS__, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

// We also provide this helpful macro for instantiating both the base case and the RC test.
#define SF_RC_TEST_WITH_BASE_CASE(SUITE_NAME, TEST_NAME, BASE_CASE_PARAM, ARBITRARY_PARAM, FN_PARAM, TAGS) \
void test_ ## SUITE_NAME ## _ ## TEST_NAME ## _impl(FN_PARAM, bool); /* Fwd. */ \
_Pragma("GCC diagnostic push") \
_Pragma("GCC diagnostic ignored \"-Wunused-parameter\"") /* Ignore -Wunused-parameter in case SF_RC_TAG is a no-op. */ \
inline void test_ ## SUITE_NAME ## _ ## TEST_NAME ## _tag(FN_PARAM) { /* Tag here since we don't know the variable name in FN_PARAM. */ \
	if (SF_NARGS TAGS) { \
		SF_RC_TAG TAGS; \
	} \
} \
_Pragma("GCC diagnostic pop") \
TEST(SUITE_NAME, TEST_NAME ## BaseCase) { /* GoogleTest test */ \
	BASE_CASE_PARAM val{}; \
	test_ ## SUITE_NAME ## _ ## TEST_NAME ## _impl(val, false); \
} \
RC_GTEST_PROP(SUITE_NAME, TEST_NAME, (ARBITRARY_PARAM val)) { /* RapidCheck test */ \
	test_ ## SUITE_NAME ## _ ## TEST_NAME ## _tag(val); \
	test_ ## SUITE_NAME ## _ ## TEST_NAME ## _impl(val, true); \
} \
void test_ ## SUITE_NAME ## _ ## TEST_NAME ## _impl(FN_PARAM, [[maybe_unused]] bool const sf_test_uses_rapidcheck) /* Test function body follows. */
