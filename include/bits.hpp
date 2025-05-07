/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <bit>
#include <concepts>

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

	template <std::unsigned_integral t_unsigned>
	constexpr inline t_unsigned pext(t_unsigned src, t_unsigned mask)
	{
		// A linear in the number of bits set in the mask implementation of PEXT.
		// FIXME: Needs to be tested.

		t_unsigned retval{};
		t_unsigned dst_idx{};

		while (mask)
		{
			t_unsigned const trailing_zeros(std::countr_zero(mask)); // Find the next set bit in the mask.
			src >>= trailing_zeros; // Shift the source and the mask to this position.
			mask >>= trailing_zeros;
			retval |= (src & t_unsigned(1)) << dst_idx; // Extract.
			src >>= t_unsigned(1); // For src, needs to be done after reading the masked bit.
			mask >>= t_unsigned(1); // This avoids UB (w.r.t. mask >>= trailing_zeros + 1) when only the most significant bit of the mask is set.
			++dst_idx;
		}

		return retval;
	}
}


namespace sf::bits {

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
}
