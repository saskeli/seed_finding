#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include "pack_characters.hpp"
#include "packed_character_iteration.hpp"


namespace sf {

	class reader_adapter; // Fwd.


	struct reader_adapter_delegate
	{
		virtual ~reader_adapter_delegate() {}
		virtual bool should_report_errors_for_path(reader_adapter &adapter, std::string_view path) = 0;
		virtual void found_first_read_with_unexpected_character(reader_adapter &adapter, std::string_view path, std::uint64_t lineno) = 0;
		virtual void found_total_reads_with_unexpected_characters(reader_adapter &adapter, std::string_view path, std::uint64_t count) = 0;
	};


	class reader_adapter
	{
	public:
		typedef packed_word_vector read_buffer_type;

	protected:
		read_buffer_type			m_read_buffer;
		reader_adapter_delegate		*m_delegate{};
		std::uint64_t				m_read_length{};

	public:
		constexpr reader_adapter() = default;

		explicit reader_adapter(reader_adapter_delegate &delegate):
			m_delegate(&delegate)
		{
		}

		virtual ~reader_adapter() {}

		virtual void read_from_path(char const *path) = 0;
		void read_from_path(std::string const &path) { read_from_path(path.c_str()); }
		virtual bool retrieve_next_read() = 0;
		virtual void finish() = 0;

		std::uint64_t read_length() const { return m_read_length; }
		read_buffer_type const &read_buffer() const { return m_read_buffer; }

		template <typename t_cb>
		void iterate_characters(std::uint64_t start_pos, t_cb &&cb) const { iterate_packed_characters(m_read_buffer, m_read_length, start_pos, cb); }

		template <typename t_cb>
		void iterate_character_pairs(std::uint64_t lhs_start, std::uint64_t rhs_start, t_cb &&cb) const { iterate_packed_character_pairs(m_read_buffer, m_read_length, lhs_start, rhs_start, cb); }
	};
}
