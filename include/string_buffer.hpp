/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

 // The purpose of this struct to make it easy to align strings to a buffer that is aligned to some number of bytes.
 // An option would have been to use a custom allocator for std::string but unfortunately the rules concerning
 // accessing values of some type that is not char or std::byte are somewhat unclear to me.

#pragma once

#include <algorithm>
#include <climits>
#include <cstddef>
#include <ostream>
#include <string_view>
#include <vector>


namespace sf {

	template <typename t_type>
	class string_buffer
	{
		std::vector <t_type>	data_;
		std::size_t				size_{};

	public:
		t_type const *data() const { return data_.data(); }
		std::size_t size() const { return size_; }
		char operator[](std::size_t idx) const { return 0xff & (data_[idx / sizeof(t_type)] >> (idx % sizeof(t_type) * CHAR_BIT)); }
		void append(std::string_view sv);
		std::string_view to_string_view() const { return {reinterpret_cast <char const *>(data_.data()), size_}; }
		/* implicit */ operator std::string_view() const { return to_string_view(); }
		string_buffer &operator=(std::string_view sv);
	};


	template <typename t_type>
	auto string_buffer <t_type>::operator=(std::string_view sv) -> string_buffer &
	{
		data_.clear();
		size_ = 0;
		append(sv);
		return *this;
	}


	template <typename t_type>
	void string_buffer <t_type>::append(std::string_view sv)
	{
		auto const added_length(sv.size());
		if (!added_length)
			return;

		auto const new_buffer_size((size_ + added_length - 1) / sizeof(t_type) + 1);
		data_.resize(new_buffer_size, 0);

		auto *dst(reinterpret_cast <char *>(data_.data()));
		dst += size_;
		std::copy_n(sv.data(), added_length, reinterpret_cast <char *>(dst));

		size_ += added_length;
	}


	template <typename t_type>
	std::ostream &operator<<(std::ostream &os, string_buffer <t_type> const &sb)
	{
		os << sb.to_string_view();
		return os;
	}
}
