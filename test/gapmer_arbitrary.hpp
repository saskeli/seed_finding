/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <bit>
#include <cstdint>
#include <gtest/gtest.h>
#include <rapidcheck/gtest.h>
#include <vector>
#include "../include/gapmer.hpp"
#include "../include/util.hpp"


namespace {

	struct gapmer_data
	{
		uint64_t sequence{}; // Half-nibbles in reverse order like in gapmer.
		uint8_t length{};
		uint8_t suffix_length{};
		uint8_t gap_length{};

		void write_to_buffer(std::vector <uint64_t> &buffer) const;

		uint8_t gap_start() const { return suffix_length ? length - suffix_length : 0; }
	};


	void gapmer_data::write_to_buffer(std::vector <uint64_t> &buffer) const
	{
		if (!length)
			return;

		buffer.resize((length + gap_length - 1) / 8 + 1, 0);
		auto *dst(reinterpret_cast <char *>(buffer.data()));
		auto sequence_{std::rotr(sequence, 2 * (length - 1))};
		uint8_t i{};

		while (i < gap_start())
		{
			dst[i] = sf::v_to_nuc[sequence_ & 0x3];
			sequence_ = std::rotl(sequence_, 2);
			++i;
		}

		for (uint8_t j{}; j < gap_length; ++j)
		{
			dst[i] = '.';
			++i;
		}

		while (i < length + gap_length)
		{
			dst[i] = sf::v_to_nuc[sequence_ & 0x3];
			sequence_ = std::rotl(sequence_, 2);
			++i;
		}
	}
}


namespace rc {

	template <>
	struct Arbitrary <gapmer_data>
	{
		static Gen <gapmer_data> arbitrary()
		{
			// We first determine the k-mer length. If it is at least two, we allow a non-zero suffix length.
			// If we got a non-empty suffix, we determine a non-zero gap length.

			// FIXME: Test different instantiations.
			typedef sf::gapmer <> gapmer_type;

			return gen::mapcat(gen::inRange <uint8_t>(1, gapmer_type::max_k + 1), [](auto const length){
				auto suffix_length_gen(length ? gen::inRange <uint8_t>(0, length) : gen::just(uint8_t(0)));
				return gen::mapcat(suffix_length_gen, [length](auto const suffix_length){ // Suffix length depends on total length.
					auto gap_length_gen(suffix_length ? gen::inRange <uint8_t>(1, gapmer_type::max_gap + 1) : gen::just(uint8_t(0)));
					return gen::mapcat(gap_length_gen, [length, suffix_length](auto const gap_length){ // Gap length depends on suffix length.
						auto const value_limit(uint64_t(1) << (2 * length));
						return gen::construct <gapmer_data>(
							gen::inRange(uint64_t(0), value_limit),
							gen::just(length),
							gen::just(suffix_length),
							gen::just(gap_length)
						);
					});
				});
			});
		}
	};
}


namespace sf {

	RC_GTEST_PROP(gapmer_arbitrary, constructorsWorkAsExpected, (gapmer_data const &gd)) {

		std::vector <uint64_t> seq;
		gd.write_to_buffer(seq);

		gapmer const gg1(gd.sequence, gd.length, gd.gap_start(), gd.gap_length);

		// Compare gg1 to gd.
		RC_ASSERT(gd.sequence == gg1.value());
		RC_ASSERT(gd.length == gg1.length());
		RC_ASSERT(gd.gap_start() == gg1.gap_start());
		RC_ASSERT(gd.gap_length == gg1.gap_length());

		// FIXME: Currently this constructor only works for non-empty sequences.
		if (gd.length)
		{
			gapmer const gg2(seq.data(), gd.length, gd.gap_start(), gd.gap_length);
			RC_ASSERT(gg1 == gg2);
		}

		if (!gd.suffix_length)
		{
			gapmer const gg3(gd.sequence, gd.length);
			RC_ASSERT(gg1 == gg3);

			// FIXME: Currently this constructor only works for non-empty sequences.
			if (gd.length)
			{
				gapmer const gg4(seq.data(), gd.length);
				RC_ASSERT(gg1 == gg4);
			}
		}
	}
}
