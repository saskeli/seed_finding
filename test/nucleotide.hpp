/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include "test.hpp"


namespace {

	struct nucleotide
	{
		char value{};

		/* implicit */ operator char() const { return value; }
	};


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
}
