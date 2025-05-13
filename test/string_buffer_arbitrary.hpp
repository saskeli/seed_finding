/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <cstddef>
#include <cstdint>
#include <gtest/gtest.h>
#include <rapidcheck/gtest.h>
#include <string>
#include <string_view>
#include "../include/string_buffer.hpp"


namespace sf { namespace tests {

	typedef string_buffer <uint64_t> string_buffer_type; // FIXME: test uint16_t, uint32_t?

	RC_GTEST_PROP(string_buffer_arbitrary, constructorsWorkAsExpected, ()) {
		string_buffer_type sb;
		RC_ASSERT(sb.empty());
		RC_ASSERT(0 == sb.size());
	}

	RC_GTEST_PROP(string_buffer_arbitrary, assignWorksAsExpected, (std::string const &str)) {
		string_buffer_type sb;
		sb = str;
		RC_ASSERT(std::string_view(sb) == str);
	}

	RC_GTEST_PROP(string_buffer_arbitrary, clearWorksAsExpected, (std::string const &str)) {
		string_buffer_type sb;
		sb = str;
		sb.clear();
		RC_ASSERT(sb.empty());
		RC_ASSERT(0 == sb.size());
	}

	RC_GTEST_PROP(string_buffer_arbitrary, appendWorksAsExpected, (std::string const &str)) {
		string_buffer_type sb;
		sb += "test";
		sb += str;

		std::string expected("test");
		expected += str;

		RC_ASSERT(std::string_view(sb) == expected);
	}

	RC_GTEST_PROP(string_buffer_arbitrary, subscriptOperatorWorksAsExpected, (std::string const &str)) {
		string_buffer_type sb;
		sb = str;

		for (std::size_t i{}; i < str.size(); ++i)
			RC_ASSERT(sb[i] == str[i]);
	}
}}
