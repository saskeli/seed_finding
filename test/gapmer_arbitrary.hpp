/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <array>
#include <bit>
#include <cstdint>
#include <iostream>
#include <ostream>
#include <stdexcept>
#include <vector>
#include "test.hpp"
#include "../include/gapmer.hpp"
#include "../include/string_buffer.hpp"
#include "../include/util.hpp"


namespace {

	struct nucleotide
	{
		char value{};

		/* implicit */ operator char() const { return value; }
	};


	struct gapmer_data
	{
		// FIXME: Test different instantiations.
		typedef sf::gapmer <> gapmer_type;

		uint64_t sequence{}; // Half-nibbles in reverse order like in gapmer.
		uint8_t length{};
		uint8_t suffix_length{};
		uint8_t gap_length{};

		gapmer_type to_gapmer() const;
		operator gapmer_type() const { return to_gapmer(); }
		void write_to_buffer(std::vector <uint64_t> &buffer) const;
		uint8_t gap_start() const { return suffix_length ? length - suffix_length : 0; }
		uint8_t suffix_start() const { return suffix_length ? gap_start() + gap_length : 0; }
	};


	constexpr static inline uint8_t const MAX_SUFFIX_LENGTH_DEFAULT{UINT8_MAX - 1};


	// Parametrised for rc::Arbitrary.
	template <uint8_t t_min_length = 0, uint8_t t_min_suffix_length = 0, uint8_t t_max_suffix_length = MAX_SUFFIX_LENGTH_DEFAULT>
	struct gapmer_data_ : public gapmer_data
	{
		static_assert(0 == t_min_suffix_length || t_min_suffix_length < t_min_length);
		static_assert(MAX_SUFFIX_LENGTH_DEFAULT == t_max_suffix_length || ((0 == t_max_suffix_length || t_max_suffix_length < t_min_length) && t_min_suffix_length <= t_max_suffix_length));
	};


	struct gapmer_pair
	{
		typedef gapmer_data::gapmer_type gapmer_type;
		gapmer_data source{};
		gapmer_data target{};
	};


	struct aligning_gapmer_pair : public gapmer_pair {};

	// Tried with a template with an enum parameter, which resulted in less-than-useful log messages;
	// hence we have the following class hierarchy.
	struct non_aligning_gapmer_pair : public gapmer_pair {};
	struct non_aligning_gapmer_pair_value_mismatch : public non_aligning_gapmer_pair {};
	struct non_aligning_gapmer_pair_source_length_greater : public non_aligning_gapmer_pair {};
	struct non_aligning_gapmer_pair_gap_length_mismatch : public non_aligning_gapmer_pair {};
	struct non_aligning_gapmer_pair_gap_position_mismatch : public non_aligning_gapmer_pair {};


	template <typename t_gapmer = sf::gapmer <>>
	struct huddinge_neighbourhood
	{
		typedef t_gapmer gapmer_type;

		struct cmp
		{
			bool operator()(gapmer_type const lhs, gapmer_type const rhs) const { return lhs.data() < rhs.data(); }
		};

		gapmer_data gd{};

		std::set <gapmer_type, cmp> values;
	};


	// Parametrised for rc::Arbitrary.
	template <typename t_gapmer_data>
	struct huddinge_neighbourhood_ : public huddinge_neighbourhood <typename t_gapmer_data::gapmer_type>
	{
		typedef t_gapmer_data gapmer_data_type;
	};


	auto gapmer_data::to_gapmer() const -> gapmer_type
	{
		auto const retval(gapmer_type(sequence, length, gap_start(), gap_length));
		if (!retval.is_valid())
		{
			std::cerr << "sequence: " << sequence << " length: " << +length << " suffix_length: " << +suffix_length << " gap_length: " << +gap_length << '\n';
			RC_FAIL("Invalid gapmer");
		}
		return retval;
	}


	void gapmer_data::write_to_buffer(std::vector <uint64_t> &buffer) const
	{
		if (!length)
			return;

		buffer.resize((length + gap_length - 1) / 8 + 1, 0);
		auto *dst(reinterpret_cast <char *>(buffer.data()));
		auto sequence_{std::rotr(sequence, 2 * (length - 1))};
		uint8_t i{};

		auto const output_to_limit([&](uint8_t const limit){
			while (i < limit)
			{
				dst[i] = sf::v_to_nuc[sequence_ & 0x3];
				sequence_ = std::rotl(sequence_, 2);
				++i;
			}
		});

		output_to_limit(gap_start());

		for (uint8_t j{}; j < gap_length; ++j)
		{
			dst[i] = '.';
			++i;
		}

		output_to_limit(length + gap_length);
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


	template <typename t_gapmer>
	std::ostream &operator<<(std::ostream &os, huddinge_neighbourhood <t_gapmer> const &hn)
	{
		os << "source: '" << hn.gd << "' values.size(): " << hn.values.size();
		return os;
	}


	char complement_nt(char const cc)
	{
		switch (cc)
		{
			case 'A': return 'T';
			case 'C': return 'G';
			case 'G': return 'C';
			case 'T': return 'A';
			case '.': return '.';
			default:
				throw std::runtime_error("Unexpected character");
		}
	}
}


namespace rc {

	template <>
	struct Arbitrary <nucleotide>
	{
		static Gen <nucleotide> arbitrary()
		{
			constexpr static std::array const values{'A', 'C', 'G', 'T'};
			return gen::map(gen::elementOf(values), [](auto const cc) -> nucleotide {
				return {cc};
			});
		}
	};


	template <>
	struct Arbitrary <gapmer_data>
	{
		static Gen <gapmer_data> arbitrary()
		{
			return gen::just(gapmer_data{});
		}
	};


	template <uint8_t t_min_length, uint8_t t_min_suffix_length, uint8_t t_max_suffix_length>
	struct Arbitrary <gapmer_data_ <t_min_length, t_min_suffix_length, t_max_suffix_length>>
	{
		static Gen <gapmer_data_ <t_min_length, t_min_suffix_length, t_max_suffix_length>> arbitrary()
		{
			// We first determine the k-mer length. If it is at least two, we allow a non-zero suffix length.
			// If we got a non-empty suffix, we determine a non-zero gap length.

			typedef gapmer_data_ <t_min_length, t_min_suffix_length, t_max_suffix_length> return_type;
			typedef typename return_type::gapmer_type gapmer_type;

			return gen::mapcat(gen::inRange <uint8_t>(t_min_length, gapmer_type::max_k + 1), [](auto const length){
				auto suffix_length_gen(length ? gen::inRange(+t_min_suffix_length, std::min(+length, 1 + t_max_suffix_length)) : gen::just(0));
				return gen::mapcat(suffix_length_gen, [length](uint8_t const suffix_length){ // Suffix length depends on total length.
					auto gap_length_gen(suffix_length ? gen::inRange <uint8_t>(1, gapmer_type::max_gap + 1) : gen::just(uint8_t{}));
					return gen::mapcat(gap_length_gen, [length, suffix_length](auto const gap_length){ // Gap length depends on suffix length.
						auto const value_limit(uint64_t(1) << (2 * length));
						return gen::map(gen::inRange(uint64_t{}, value_limit), [length, suffix_length, gap_length](auto const sequence){
							return_type const retval{sequence, length, suffix_length, gap_length};
							return retval;
						});
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
			return gen::mapcat(gen::arbitrary <gapmer_data_ <>>(), [](auto const gd){
				return gen::mapcat(gd.length ? gen::inRange(uint8_t{}, gd.length) : gen::just(uint8_t{}), [gd](auto const pos){ // Determine the source offset in target.
					return gen::map(gen::inRange(0, gd.length - pos + 1), [gd, pos](uint8_t const length){ // Determine the source length.
						if (!length)
							return aligning_gapmer_pair{gapmer_data{}, gd};

						// Copy the sequence to get an aligning pair.
						auto const source_end(pos + length);
						auto sequence(gd.sequence);
						sequence >>= 2 * (gd.length - source_end);

						// Zero the non-sequence characters.
						{
							uint64_t mask{};
							mask = ~mask;
							mask >>= 2 * (32 - length);
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


	template <>
	struct Arbitrary <non_aligning_gapmer_pair_value_mismatch>
	{
		static Gen <non_aligning_gapmer_pair_value_mismatch> arbitrary()
		{
			return gen::mapcat(gen::arbitrary <gapmer_data_ <1>>(), [](auto const gd){
				return gen::map(gen::inRange(0, 2 * gd.length), [gd](auto const pos) -> non_aligning_gapmer_pair_value_mismatch {
					// Make a new value by flipping some bit in the value.
					// (Currently all bit patterns are valid in the range specified by length.)
					auto gd_(gd);
					uint64_t mask{1};
					mask <<= pos;
					gd_.sequence ^= mask;
					RC_ASSERT(gd.sequence != gd_.sequence);
					return {gd, gd_};
				});
			});
		}
	};


	template <>
	struct Arbitrary <non_aligning_gapmer_pair_source_length_greater>
	{
		static Gen <non_aligning_gapmer_pair_source_length_greater> arbitrary()
		{
			return gen::map(gen::arbitrary <gapmer_data_ <>>(), [](auto gd) -> non_aligning_gapmer_pair_source_length_greater {
				// Make the source longer than the target.
				auto gd_(gd);
				if (gd.length < gapmer_data::gapmer_type::max_k)
					++gd.length;
				else
					--gd_.length;

				return {gd, gd_};
			});
		}
	};


	template <>
	struct Arbitrary <non_aligning_gapmer_pair_gap_length_mismatch>
	{
		static Gen <non_aligning_gapmer_pair_gap_length_mismatch> arbitrary()
		{
			return gen::map(gen::arbitrary <gapmer_data_ <2, 1>>(), []<typename t_gapmer_data>(t_gapmer_data const gd) -> non_aligning_gapmer_pair_gap_length_mismatch {
				typedef typename t_gapmer_data::gapmer_type gapmer_type;
				static_assert(1 <= gapmer_type::max_gap);

				auto gd_(gd);
				if (gd_.gap_length < gapmer_type::max_gap)
					++gd_.gap_length;
				else
					--gd_.gap_length;

				return {gd, gd_};
			});
		}
	};


	template <>
	struct Arbitrary <non_aligning_gapmer_pair_gap_position_mismatch>
	{
		static Gen <non_aligning_gapmer_pair_gap_position_mismatch> arbitrary()
		{
			return gen::map(gen::arbitrary <gapmer_data_ <3, 1>>(), []<typename t_gapmer_data>(t_gapmer_data const gd) -> non_aligning_gapmer_pair_gap_position_mismatch {
				auto gd_(gd);
				if (1 == gd_.suffix_length)
					++gd_.suffix_length;
				else
					--gd_.suffix_length;
				return {gd, gd_};
			});
		}
	};


	template <typename t_gapmer_data>
	struct Arbitrary <huddinge_neighbourhood_ <t_gapmer_data>>
	{
		static Gen <huddinge_neighbourhood_ <t_gapmer_data>> arbitrary()
		{
			return gen::map(gen::arbitrary <t_gapmer_data>(), [](auto const gd) -> huddinge_neighbourhood_ <t_gapmer_data> {
				// Given a gapmer_data, we generate its Huddinge neighbourhood by replacing its defined
				// characters and the first and last middle gap characters one by one. We then proceed to
				// prepend and append one character.
				// We try to ensure correctness at the cost of speed.
				typedef typename huddinge_neighbourhood_ <t_gapmer_data>::gapmer_type gapmer_type;
				constexpr static std::array const characters{'A', 'C', 'G', 'T'};

				huddinge_neighbourhood_ <t_gapmer_data> retval;
				retval.gd = gd;
				
				sf::string_buffer <uint64_t> str;
				auto const gg(gd.to_gapmer());
				auto const gap_start(gg.gap_start());
				auto const gap_length(gg.gap_length());
				str = gg.to_string();
				auto span(str.to_span());

				auto const modify_([&](auto it){
					auto &cc(*it);
					auto const orig_cc(cc);
					for (auto const cc_ : characters)
					{
						cc = cc_;
						gapmer_type const gg_(str.data(), span.size(), gap_start, gap_length);
						retval.values.insert(gg_);
					}
					cc = orig_cc;
				});

				auto const modify([&](auto it, auto end){
					while (it != end)
					{
						modify_(it);
						++it;
					}
				});

				switch (gap_length)
				{
					case 0:
					case 1:
						modify(span.begin(), span.end());
						break;

					default:
						modify(span.begin(), span.begin() + gap_start + 1);
						modify(span.begin() + gap_start + gap_length - 1, span.end());
						break;
				}

				// Implicit gaps at each end.
				str.append(" ");
				span = str.to_span();
				modify_(span.rbegin());

				if (gd.length)
				{
					str >>= 1;
					modify_(span.begin());
				}

				return retval;
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


	SF_RC_TEST_WITH_BASE_CASE(gapmer_arbitrary, Constructors, gapmer_data, gapmer_data_ <1>, gapmer_data const gd, (gd.length, gd.gap_length)) {
		std::vector <uint64_t> seq;
		gd.write_to_buffer(seq);

		gapmer const gg1(gd.sequence, gd.length, gd.gap_start(), gd.gap_length);

		// Compare gg1 to gd.
		SF_ASSERT(gd.sequence == gg1.value());
		SF_ASSERT(gd.length == gg1.length());
		SF_ASSERT(gd.gap_start() == gg1.gap_start());
		SF_ASSERT(gd.gap_length == gg1.gap_length());

		// FIXME: Currently this constructor only works for non-empty sequences.
		if (gd.length)
		{
			gapmer const gg2(seq.data(), gd.length, gd.gap_start(), gd.gap_length);
			SF_ASSERT(gg1 == gg2);
		}

		if (!gd.suffix_length)
		{
			gapmer const gg3(gd.sequence, gd.length);
			SF_ASSERT(gg1 == gg3);

			// FIXME: Currently this constructor only works for non-empty sequences.
			if (gd.length)
			{
				gapmer const gg4(seq.data(), gd.length);
				SF_ASSERT(gg1 == gg4);
			}
		}
	}


	SF_RC_TEST_WITH_BASE_CASE(gapmer_arbitrary, ReverseComplementGapmer, gapmer_data, gapmer_data_ <1>, gapmer_data const gd, (gd.length, gd.gap_length)) {
		typedef gapmer_data::gapmer_type gapmer_type;
		gapmer_type const source{gd};
		auto const target(source.reverse_complement());
		auto const source_(source.to_string());
		auto const target_(target.to_string());

		auto const length(source_.size());
		SF_ASSERT(length == target_.size());
		for (std::size_t i{}; i < length; ++i)
			SF_ASSERT(source_[i] == complement_nt(target_[length - i - 1]));
	}


	RC_GTEST_PROP(gapmer_arbitrary, NextGapmer, (gapmer_data_ <0, 0, 0> const gd, nucleotide const nt)) {
		SF_RC_TAG(gd.length, gd.gap_length);

		typedef gapmer_data::gapmer_type gapmer_type;
		gapmer_type const source{gd};
		auto const target{source.next(nt)};
		auto source_{source.to_string()};
		auto const target_{target.to_string()};
		source_.push_back(nt);
		source_.erase(0, 1);
		RC_ASSERT(source_ == target_);
	}


	RC_GTEST_PROP(gapmer_arbitrary, NextGappedGapmer, (gapmer_data_ <2, 1> const gd, nucleotide const nt1, nucleotide const nt2)) {
		SF_RC_TAG(gd.length, gd.gap_length);

		typedef gapmer_data::gapmer_type gapmer_type;
		gapmer_type const source{gd};
		auto const target{source.next(nt1, nt2)};
		auto source_{source.to_string()};
		auto const target_{target.to_string()};
		source_.insert(gd.gap_start(), 1, nt1);
		source_.push_back(nt2);
		source_.erase(0, 1);
		source_.erase(gd.suffix_start(), 1);
		RC_ASSERT(source_ == target_);
	}


	SF_RC_TEST_WITH_BASE_CASE(gapmer_arbitrary, AlignEmptyGapmer, gapmer_data, gapmer_data_ <1>, gapmer_data const gd, (gd.length, gd.gap_length)) {
		typedef gapmer_data::gapmer_type gapmer_type;
		gapmer_type const source{};
		gapmer_type const target(gd);
		SF_ASSERT((source.aligns_to <false>(target))); // FIXME: also reverse complement (if needed).
	}


	RC_GTEST_PROP(gapmer_arbitrary, AlignNonemptyGapmer, (aligning_gapmer_pair const &pair)) {
		SF_RC_TAG(pair.target.length, pair.target.gap_length);

		typedef aligning_gapmer_pair::gapmer_type gapmer_type;
		gapmer_type const source{pair.source};
		gapmer_type const target{pair.target};
		RC_ASSERT((source.aligns_to <false>(target))); // FIXME: also reverse complement (if needed).
	}


	SF_RC_TEMPLATE_TEST(
		gapmer_arbitrary_aligns_to,
		AlignNonMatchingGapmer,
		(TypeParam const &pair),
		non_aligning_gapmer_pair_value_mismatch,
		non_aligning_gapmer_pair_source_length_greater
	) {
		SF_RC_TAG(pair.target.length, pair.target.gap_length);
		typedef typename TypeParam::gapmer_type gapmer_type;
		gapmer_type const source{pair.source};
		gapmer_type const target{pair.target};
		RC_ASSERT((!source.template aligns_to <false>(target)));
	}


	SF_RC_TEMPLATE_TEST(
		gapmer_arbitrary_aligns_to_reversible,
		AlignNonMatchingGapmerReversible,
		(TypeParam const &pair),
		non_aligning_gapmer_pair_gap_length_mismatch,
		non_aligning_gapmer_pair_gap_position_mismatch
	) {
		SF_RC_TAG(pair.target.length, pair.target.gap_length);
		typedef typename TypeParam::gapmer_type gapmer_type;
		gapmer_type const gg1{pair.source};
		gapmer_type const gg2{pair.target};
		RC_ASSERT((!gg1.template aligns_to <false>(gg2)));
		RC_ASSERT((!gg2.template aligns_to <false>(gg1)));
	}


	// FIXME: Test other gapmer types.
#if 0
	SF_RC_TEST_WITH_BASE_CASE(
		gapmer_arbitrary,
		GenerateHuddingeNeighbourhood,
		huddinge_neighbourhood_ <gapmer_data>,
		huddinge_neighbourhood_ <gapmer_data_ <1>>,
		huddinge_neighbourhood <> const &hn,
		(hn.gd.length, hn.gd.gap_length)
	)
#else
	SF_RC_TEST_PROP(
		gapmer_arbitrary,
		GenerateHuddingeNeighbourhood,
		huddinge_neighbourhood_ <gapmer_data_ <5>>,
		huddinge_neighbourhood <> const &hn,
		(hn.gd.length, hn.gd.gap_length)
	)
#endif
	{
		typedef huddinge_neighbourhood <>::gapmer_type gapmer_type; // FIXME: Use a type parameter here.

		auto const gd(hn.gd);
		gapmer_type const gg{gd};

		std::set <gapmer_type, huddinge_neighbourhood <>::cmp> actual; // FIXME: Use a type parameter here.
		gg.middle_gap_neighbours <false, false, false>([&](gapmer_type const gg_) {
			SF_EXPECT(gg_.is_valid());
			actual.insert(gg_);
		});
		SF_ASSERT(hn.values == actual);
	}
}
