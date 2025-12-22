/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#include "../include/huddinge_distance.hpp"
#include "../include/pack_characters.hpp"
#include "gtest/gtest.h"
#include "test.hpp"

namespace {

	template <std::size_t t_size>
	using span_t = std::span <std::uint64_t const, t_size>;

	typedef std::vector <std::uint64_t> buffer_type;

	constexpr static std::array const dna4_huddinge_characters{'n', 'A', 'C', 'G', 'T'};


	struct dna4_huddinge_string_length
	{
		std::size_t value{};

		operator std::size_t() const { return value; }
	};


	struct dna4_huddinge_string_pair_length
	{
		dna4_huddinge_string_length lhs;
		dna4_huddinge_string_length rhs;
	};


	struct dna4_huddinge_string_pair
	{
		constexpr static auto const &characters{dna4_huddinge_characters};

		std::string lhs;
		std::string rhs;
	};


	std::ostream &operator<<(std::ostream &os, dna4_huddinge_string_pair const &sp)
	{
		os << "lhs: " << sp.lhs << " rhs: " << sp.rhs;
		return os;
	}


	std::size_t count_defined_characters(std::string_view sv)
	{
		return std::count_if(sv.begin(), sv.end(), [](auto const cc){ return 'n' != cc; });
	}


	template <typename t_eq_fn>
	sf::huddinge_distance_return_value huddinge_distance_naive_(std::string_view lhs, std::string_view rhs, t_eq_fn &&eq_fn)
	{
		sf::huddinge_distance_return_value retval;

		auto const lhs_defined_characters{count_defined_characters(lhs)};
		auto const rhs_defined_characters{count_defined_characters(rhs)};
		auto const max_defined_characters{std::max(lhs_defined_characters, rhs_defined_characters)};

		// Test the positions where lhs has been shifted to the left w.r.t. rhs.
		for (std::size_t ii{1}; ii < lhs.size(); ++ii)
		{
			std::uint64_t current_matches{};
			for (std::size_t jj{}; jj < rhs.size() and ii + jj < lhs.size(); ++jj)
				current_matches += eq_fn(lhs[ii + jj], rhs[jj]) ? 1 : 0;

			auto const current_distance{max_defined_characters - current_matches};
			if (current_distance < retval.distance)
			{
				retval.distance = current_distance;
				retval.position = ii;
			}
		}

		// Test the positions where lhs has been shifted to the right or is at position zero w.r.t. rhs.
		for (std::size_t ii{}; ii < rhs.size(); ++ii)
		{
			std::uint64_t current_matches{};
			for (std::size_t jj{}; jj < lhs.size() and ii + jj < rhs.size(); ++jj)
				current_matches += eq_fn(lhs[jj], rhs[ii + jj]) ? 1 : 0;

			auto const current_distance{max_defined_characters - current_matches};
			if (current_distance < retval.distance)
			{
				retval.distance = current_distance;
				retval.position = -ii;
			}
		}

		return retval;
	}


	sf::huddinge_distance_return_value huddinge_distance_naive(std::string_view lhs, std::string_view rhs)
	{
		// Hamming distance with equals.
		return huddinge_distance_naive_(lhs, rhs, [](auto const ll, auto const rr) -> bool {
			if (ll == rr)
				return 'n' != ll;
			return false;
		});
	}


	std::vector <uint64_t> make_gap_mask(std::string_view sv)
	{
		std::vector <uint64_t> retval;
		std::size_t pos{};

		for (auto const cc : sv)
		{
			if (0 == pos % 32)
				retval.push_back(0);

			if ('n' != cc)
			{
				uint64_t const shift{62 - 2 * (pos % 32)};
				retval.back() |= UINT64_C(0x1) << shift;
			}

			++pos;
		}

		return retval;
	}


	sf::huddinge_distance_return_value huddinge_distance(std::string_view lhs, std::string_view rhs)
	{
		auto const lhs_mask{make_gap_mask(lhs)};
		auto const rhs_mask{make_gap_mask(rhs)};

		buffer_type lhs_buffer;
		buffer_type rhs_buffer;
		sf::pack_characters_lenient(lhs, lhs_buffer);
		sf::pack_characters_lenient(rhs, rhs_buffer);

		auto const forward_spans([&](std::string_view sv, buffer_type const &buffer, buffer_type const &mask, auto &&cb){
			if (1 <= sv.size() && sv.size() <= 32)
				return cb(span_t <1>(buffer.data(), 1), span_t <1>(mask.data(), 1));
			else if (33 <= sv.size() && sv.size() <= 64)
				return cb(span_t <2>(buffer.data(), 2), span_t <2>(mask.data(), 2));
			else
				throw std::runtime_error("Unexpected input size");
		});

		return forward_spans(lhs, lhs_buffer, lhs_mask, [&](auto const lhs_buffer_, auto const lhs_mask_){
			return forward_spans(rhs, rhs_buffer, rhs_mask, [&](auto const rhs_buffer_, auto const rhs_mask_){
				return sf::huddinge_distance(lhs_buffer_, rhs_buffer_, lhs_mask_, rhs_mask_, lhs.size(), rhs.size());
			});
		});
	}
}


namespace rc {

	template <>
	struct Arbitrary <dna4_huddinge_string_length>
	{
		static Gen <dna4_huddinge_string_length> arbitrary()
		{
			// Generate non-empty strings with 64 characters at most.
			return gen::construct <dna4_huddinge_string_length>(
				gen::inRange(std::size_t{1}, std::size_t{64})
			);
		}
	};


	template <>
	struct Arbitrary <dna4_huddinge_string_pair_length>
	{
		static Gen <dna4_huddinge_string_pair_length> arbitrary()
		{
			// Generate non-empty strings with 64 characters at most.
			return gen::construct <dna4_huddinge_string_pair_length>(
				gen::arbitrary <dna4_huddinge_string_length>(),
				gen::arbitrary <dna4_huddinge_string_length>()
			);
		}
	};


	template <>
	struct Arbitrary <dna4_huddinge_string_pair>
	{
		static Gen <dna4_huddinge_string_pair> arbitrary()
		{
			return gen::mapcat(gen::arbitrary <dna4_huddinge_string_pair_length>(), [](auto const lengths){
				return gen::construct <dna4_huddinge_string_pair>(
					gen::container <std::string>(lengths.lhs, gen::elementOf(dna4_huddinge_string_pair::characters)),
					gen::container <std::string>(lengths.rhs, gen::elementOf(dna4_huddinge_string_pair::characters))
				);
			});
		}
	};
}


RC_GTEST_PROP(huddinge_distance_arbitrary, CalculateDistance, (dna4_huddinge_string_pair sp)) {
	auto const expected{huddinge_distance_naive(sp.lhs, sp.rhs)};
	auto const actual{huddinge_distance(sp.lhs, sp.rhs)};
	ASSERT_EQ(expected.distance, actual.distance) << "  lhs:      " << sp.lhs << "\n  rhs:      " << sp.rhs << "\n  expected: " << expected << "\n  actual:   " << actual;
}
