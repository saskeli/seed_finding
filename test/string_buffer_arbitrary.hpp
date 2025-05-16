/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include "../include/string_buffer.hpp"
#include "test.hpp"


namespace {
	struct non_empty_string
	{
		std::string value;

		operator std::string const &() const { return value; }
	};
}


namespace rc {
	template <>
	struct Arbitrary <non_empty_string>
	{
		static Gen <non_empty_string> arbitrary()
		{
			return gen::construct <non_empty_string>(rc::gen::nonEmpty <std::string>());
		}
	};
}


namespace sf { namespace tests {

	typedef string_buffer <uint64_t> string_buffer_type; // FIXME: test uint16_t, uint32_t?


	TEST(string_buffer_arbitrary, Constructor)
	{
		string_buffer_type sb;
		RC_ASSERT(sb.empty());
		RC_ASSERT(0 == sb.size());
	}

	SF_RC_TEST_WITH_BASE_CASE(string_buffer_arbitrary, Assign, std::string, non_empty_string, std::string const &str, (str.size())) {
		string_buffer_type sb;
		sb = str;
		RC_ASSERT(std::string_view(sb) == str);
	}

	SF_RC_TEST_WITH_BASE_CASE(string_buffer_arbitrary, Clear, std::string, non_empty_string, std::string const &str, (str.size())) {
		string_buffer_type sb;
		sb = str;
		sb.clear();
		RC_ASSERT(sb.empty());
		RC_ASSERT(0 == sb.size());
	}

	SF_RC_TEST_WITH_BASE_CASE(string_buffer_arbitrary, Append, std::string, non_empty_string, std::string const &str, (str.size())) {
		string_buffer_type sb;
		sb += "test";
		sb += str;

		std::string expected("test");
		expected += str;

		RC_ASSERT(std::string_view(sb) == expected);
	}

	RC_GTEST_PROP(string_buffer_arbitrary, subscriptOperatorWorksAsExpected, (non_empty_string const &str_)) {
		std::string const &str(str_);
		SF_RC_TAG(str.length());

		string_buffer_type sb;
		sb = str;

		for (std::size_t i{}; i < str.size(); ++i)
			RC_ASSERT(sb[i] == str[i]);
	}
}}
