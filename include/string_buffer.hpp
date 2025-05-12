/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

 // The purpose of this struct to make it easy to align strings to a buffer that is aligned to some number of bytes.
 // An option would have been to use a custom allocator for std::string but unfortunately the rules concerning
 // accessing values of some type that is not char or std::byte are somewhat unclear to me.

#pragma once

#include <algorithm>
#include <cstddef>
#include <string_view>
#include <vector>


namespace sf {

	template <typename t_type>
	struct string_buffer
	{
		std::vector <t_type>	data;
		std::size_t				size{};

		void append(std::string_view sv);
		std::string_view to_string_view() const { return {data.data(), size}; }
		/* implicit */ operator std::string_view() const { return to_string_view(); }
	};


	template <typename t_type>
	void string_buffer <t_type>::append(std::string_view sv)
	{
		auto const added_length(sv.size());
		if (!added_length)
			return;

		auto const new_buffer_size((size + added_length - 1) / sizeof(t_type) + 1);
		data.resize(new_buffer_size, 0);

		char *dst(data.data());
		dst += size;
		std::copy_n(sv.data(), added_length, reinterpret_cast <char *>(dst));

		size += added_length;
	}
}
