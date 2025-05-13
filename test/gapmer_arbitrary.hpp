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

		sf::gapmer <> to_gapmer() const { return sf::gapmer <>(sequence, length, gap_start(), gap_length); }
		operator sf::gapmer <>() const { return to_gapmer(); }
		void write_to_buffer(std::vector <uint64_t> &buffer) const;
		uint8_t gap_start() const { return suffix_length ? length - suffix_length : 0; }
	};


	struct gapmer_pair
	{
		gapmer_data source{};
		gapmer_data target{};
	};


	struct aligning_gapmer_pair : public gapmer_pair {};
	struct non_aligning_gapmer_pair : public gapmer_pair {};


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


	std::ostream &operator<<(std::ostream &os, gapmer_data const gd)
	{
		os << gd.to_gapmer().to_string() << " (" << +gd.length << ", " << +gd.gap_start() << ", " << +gd.gap_length << ')';
		return os;
	}


	std::ostream &operator<<(std::ostream &os, gapmer_pair const &pp)
	{
		os << "source: " << pp.source << " target: " << pp.target;
		return os;
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


	template <>
	struct Arbitrary <aligning_gapmer_pair>
	{
		static Gen <aligning_gapmer_pair> arbitrary()
		{
			return gen::mapcat(gen::arbitrary <gapmer_data>(), [](auto const gd){
				return gen::mapcat(gen::inRange(uint8_t{}, gd.length), [gd](auto const pos){ // Determine the source offset in target.
					return gen::map(gen::inRange(1, gd.length - pos + 1), [gd, pos](uint8_t const length){ // Determine the source length.

						// Copy the sequence to get an aligning pair.
						auto const source_end(pos + length);
						auto sequence(gd.sequence);
						sequence >>= 2 * (gd.length - source_end);

						// Zero the non-sequence characters.
						{
							uint64_t mask{};
							mask = ~mask;
							mask >>= 32 - length;
							sequence &= mask;
						}

						// Determine the suffix length and the gap length of the source.
						auto const target_gap_start(gd.gap_start());
						uint8_t const suffix_length(pos < target_gap_start && target_gap_start < source_end ? source_end - target_gap_start : 0);
						uint8_t const gap_length(suffix_length ? gd.gap_length : 0);
						
						return aligning_gapmer_pair{
							gapmer_data{sequence, length, suffix_length, gap_length},
							gd
						};
					});
				});
			});
		}
	};
}


namespace sf {

	template <bool middle_gap_only, uint16_t t_max_gap>
	void showValue(gapmer <middle_gap_only, t_max_gap> const gg, std::ostream &os)
	{
		os << gg.to_string();
	}


	RC_GTEST_PROP(gapmer_arbitrary, constructorsWorkAsExpected, (gapmer_data const gd)) {

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


	RC_GTEST_PROP(gapmer_arbitrary, alignEmptyGapmer, (gapmer_data const gd)) {
		gapmer <> const source{};
		gapmer <> const target(gd);
		RC_ASSERT((source.aligns_to <false>(target))); // FIXME: also reverse complement (if needed).
	}


	RC_GTEST_PROP(gapmer_arbitrary, alignNonemptyGapmer, (aligning_gapmer_pair const &pair)) {
		gapmer <> const source{pair.source};
		gapmer <> const target{pair.target};
		RC_ASSERT((source.aligns_to <false>(target))); // FIXME: also reverse complement (if needed).
	}


	// FIXME: Also non-aligned. (Flip some bit of the value or adjust the gap length or position to make the operation fail.)
}
