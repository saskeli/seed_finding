/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <bit>
#include <cassert>
#include <concepts>
#include <climits>
#include <csignal>
#include <cstddef>
#include <cstdint>
#include <span>

#if defined(__linux__)
#	include <byteswap.h>
#endif

#if defined(__i386__) || defined(__x86_64__)
#	include <immintrin.h>
#endif


namespace sf::bits::detail {

#if defined(__BMI2__) // PEXT is part of BMI2.
	constexpr static inline bool const pext_intrinsic_available{true};
	inline uint32_t pext_intrinsic(uint32_t source, uint32_t mask) { return _pext_u32(source, mask); }
	inline uint64_t pext_intrinsic(uint64_t source, uint64_t mask) { return _pext_u64(source, mask); }
#else
	constexpr static inline bool const pext_intrinsic_available{false};
	inline uint64_t pext_intrinsic(uint64_t, uint64_t) { return 0; } // FIXME: I don’t remember how to use if constexpr in such a way that the discarded statement is not checked. The compiler currently checks both branches in sf::bits::pext().
#endif


#if defined(__linux__)
	inline uint16_t byteswap(uint16_t source) { return bswap_16(source); }
	inline uint32_t byteswap(uint32_t source) { return bswap_32(source); }
	inline uint64_t byteswap(uint64_t source) { return bswap_64(source); }
#endif


	template <std::unsigned_integral t_unsigned>
	constexpr inline t_unsigned pext(t_unsigned src, t_unsigned mask)
	{
		// A linear in the number of runs of set bits in the mask implementation of PEXT.

		t_unsigned retval{};
		t_unsigned dst_idx{};
		t_unsigned clear_mask{};

		clear_mask = ~clear_mask;

		while (mask)
		{
			t_unsigned const trailing_zeros(std::countr_zero(mask)); // Find the next set bit in the mask.
			src >>= trailing_zeros;                                  // Shift the source and the mask to this position.
			mask >>= trailing_zeros;
			retval |= src << dst_idx;                              // Extract.

			t_unsigned const trailing_ones(std::countr_one(mask)); // Find the next unset bit in the mask.
			src >>= trailing_ones;                                 // Skip to the end of the current run.
			mask >>= trailing_ones;
			dst_idx += trailing_ones;                              // Set up the next target position.
			clear_mask <<= trailing_ones;                          // Clear the bits outside the current run.
			retval &= ~clear_mask;
		}

		return retval;
	}
}


namespace sf::bits {

	template <std::unsigned_integral t_unsigned>
	constexpr inline t_unsigned byteswap(t_unsigned source)
	{
		if consteval
		{
#if defined(__cpp_lib_byteswap)
			return std::byteswap(source);
#else
#	error "std::byteswap not available."
#endif
		}
		else
		{
#if defined(__cpp_lib_byteswap)
				return std::byteswap(source);
#elif defined(__linux__)
				static_assert(2 <= sizeof(t_unsigned) && sizeof(t_unsigned) <= 8);
				return detail::byteswap(source);
#else
#	error "std::byteswap or equivalent not available."
#endif
		}
	}


	template <std::unsigned_integral t_unsigned>
	constexpr inline t_unsigned pext(t_unsigned source, t_unsigned mask)
	{
		if consteval
		{
			return detail::pext(source, mask);
		}
		else
		{
			if constexpr (detail::pext_intrinsic_available && (std::is_same_v <t_unsigned, std::uint32_t> || std::is_same_v <t_unsigned, std::uint64_t>))
				return detail::pext_intrinsic(source, mask);
			else
				return detail::pext(source, mask);
		}
	}

	template <std::unsigned_integral t_value, std::size_t t_n, std::unsigned_integral t_count>
	constexpr inline void shift_left(std::span <t_value, t_n> span, t_count const count)
	{
		if (!count) return;

		auto const value_bits{sizeof(t_value) * CHAR_BIT};
		assert(count < value_bits); // The general case has not been implemented.
		auto const lower_mask([&]() constexpr {
			t_value mask{};
			mask = ~mask;
			mask >>= value_bits - count;
			return mask;
		}());
		auto const higher_mask(~lower_mask);

		t_value lower{};
		for (auto &word : span)
		{
			auto word_(std::rotl(word, count));
			word = word_ & higher_mask;
			word |= lower;
			lower = word_ & lower_mask;
		}
	}

	template <std::unsigned_integral t_value, std::size_t t_n, std::unsigned_integral t_count>
	constexpr inline void shift_right(std::span <t_value, t_n> span, t_count const count)
	{
		if (!count) return;

		auto const value_bits{sizeof(t_value) * CHAR_BIT};
		assert(count < value_bits); // The general case has not been implemented.
		auto const higher_mask([&]() constexpr {
			t_value mask{};
			mask = ~mask;
			mask <<= value_bits - count;
			return mask;
		}());
		auto const lower_mask(~higher_mask);

		t_value higher{};
		for (auto it(span.rbegin()); it != span.rend(); ++it)
		{
			auto &word(*it);
			auto word_(std::rotr(word, count));
			word = word_ & lower_mask;
			word |= higher;
			higher = word_ & higher_mask;
		}
	}
}
