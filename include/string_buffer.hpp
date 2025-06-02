/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

 #pragma once

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <ostream>
#include <span>
#include <string_view>
#include <vector>


namespace sf {

	// The purpose of this class template is to make it easy to align strings to a buffer that is aligned to some number of bytes.
	// An option would have been to use a custom allocator for std::string but unfortunately the rules concerning
	// accessing values of some type that is not char or std::byte are somewhat unclear to me.
	template <typename t_type>
	class string_buffer
	{
		std::vector <t_type>	data_;
		std::size_t				size_{};

	public:
		string_buffer() = default;
		explicit string_buffer(std::string_view sv) { assign(sv); }

		t_type const *data() const { return data_.data(); }
		bool empty() const { return 0 == size_; }
		std::size_t size() const { return size_; }
		void resize(std::size_t new_size);
		void clear();
		void assign(std::string_view sv);
		void append(std::string_view sv);
		void append(char cc) { append(std::string_view{&cc, 1}); }
		void shift_left(uint64_t const amt);
		void shift_right(uint64_t const amt);
		std::span <char> to_span() { return {reinterpret_cast <char *>(data_.data()), size_}; }
		std::span <char const> to_span() const { return {reinterpret_cast <char const *>(data_.data()), size_}; }
		std::string_view to_string_view() const { return {reinterpret_cast <char const *>(data_.data()), size_}; }

		char operator[](std::size_t idx) const { return 0xff & (data_[idx / sizeof(t_type)] >> (idx % sizeof(t_type) * CHAR_BIT)); }
		char &operator[](std::size_t idx) { return to_span()[idx]; }
		/* implicit */ operator std::string_view() const { return to_string_view(); }
		string_buffer &operator=(std::string_view sv) { assign(sv); return *this; }
		string_buffer &operator+=(std::string_view sv) { append(sv); return *this; }
		string_buffer &operator<<=(uint64_t const amt) { shift_left(amt); return *this; }
		string_buffer &operator>>=(uint64_t const amt) { shift_right(amt); return *this; }
		bool operator==(std::string_view sv) const { return to_string_view() == sv; }
	};


	template <typename t_type>
	void string_buffer <t_type>::resize(std::size_t new_size)
	{
		data_.resize(new_size);
		size_ = new_size;
	}


	template <typename t_type>
	void string_buffer <t_type>::clear()
	{
		data_.clear();
		size_ = 0;
	}


	template <typename t_type>
	void string_buffer <t_type>::assign(std::string_view sv)
	{
		clear();
		append(sv);
	}


	template <typename t_type>
	void string_buffer <t_type>::append(std::string_view sv)
	{
		// Determine and set the new buffer size.
		auto const added_length(sv.size());
		auto const new_buffer_size((size_ + added_length + 7) / sizeof(t_type));
		data_.resize(new_buffer_size, 0);

		// Copy the bytes.
		auto *dst(reinterpret_cast <char *>(data_.data()));
		dst += size_;
		std::copy_n(sv.data(), added_length, dst);

		size_ += added_length;
	}


	template <typename t_type>
	void string_buffer <t_type>::shift_left(uint64_t const amt)
	{
		// Deletes characters from the left end, does not reset the characters at the right end.
		if (size_ <= amt)
		{
			size_ = 0;
			return;
		}

		auto span(to_span());
		auto const begin(span.begin());
		auto const end(span.end());
		auto const src_begin(begin + amt);
		std::copy(src_begin, end, begin);
	}


	template <typename t_type>
	void string_buffer <t_type>::shift_right(uint64_t const amt)
	{
		// Deletes characters from the right end, does not reset the characters at the left end.
		if (size_ <= amt)
		{
			size_ = 0;
			return;
		}

		auto span(to_span());
		auto const begin(span.begin());
		auto const end(span.end());
		auto const src_end(end - amt);
		std::copy_backward(begin, src_end, end);
	}


	template <typename t_type>
	std::ostream &operator<<(std::ostream &os, string_buffer <t_type> const &sb)
	{
		os << sb.to_string_view();
		return os;
	}
}
