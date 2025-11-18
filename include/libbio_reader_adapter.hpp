#pragma once

#include <cstddef>
#include <cstdint>
#include <libbio/assert.hh>
#include <libbio/fasta_reader.hh>
#include <libbio/fastq_reader.hh>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/gzip_read_handle.hh>
#include <libbio/sequence_reader.hh>
#include <span>
#include <string>
#include <string_view>
#include "packed_character_iteration.hpp"


namespace sf {

	class libbio_reader_adapter;


	struct libbio_reader_adapter_delegate
	{
		virtual ~libbio_reader_adapter_delegate() {}
		virtual bool should_report_errors_for_path(libbio_reader_adapter &adapter, std::string_view path) = 0;
		virtual void found_first_read_with_unexpected_character(libbio_reader_adapter &adapter, std::string_view path, std::uint64_t lineno) = 0;
		virtual void found_total_reads_with_unexpected_characters(libbio_reader_adapter &adapter, std::string_view path, std::uint64_t count) = 0;
	};


	class libbio_reader_adapter final : public libbio::fasta_reader_delegate, libbio::fastq_reader_delegate
	{
	public:
		typedef packed_word_vector read_buffer_type;
		typedef libbio_reader_adapter_delegate delegate_type;

	private:
		libbio::fasta_reader		m_fasta_reader;
		libbio::fastq_reader		m_fastq_reader;

		libbio::file_handle			m_handle;
		libbio::gzip_reading_handle	m_gzip_handle;
		std::string					m_current_input_path;

		read_buffer_type			m_read_buffer;
		std::uint64_t				m_read_length{};

		libbio::sequence_reader		*m_reader{};
		libbio::reading_handle		*m_reading_handle{};
		delegate_type				*m_delegate{};

		std::size_t					m_parsing_block_size{};
		std::uint64_t				m_skipped_reads{};
		bool						m_should_report_errors{};
		bool						m_read_is_valid{};
		bool						m_next_read_is_reverse_complement{};

	public:
		constexpr libbio_reader_adapter() = default;

		explicit libbio_reader_adapter(libbio_reader_adapter_delegate &delegate):
			m_delegate(&delegate)
		{
			m_fasta_reader.set_delegate(*this);
			m_fastq_reader.set_delegate(*this);
			m_handle.prepare();
			m_gzip_handle.prepare();
		}

		libbio_reader_adapter(libbio_reader_adapter const &) = delete;
		libbio_reader_adapter(libbio_reader_adapter &&) = default;
		libbio_reader_adapter &operator=(libbio_reader_adapter const &) = delete;
		libbio_reader_adapter &operator=(libbio_reader_adapter &&) & = default;

		void read_from_path(std::string const &path) { read_from_path(path.c_str()); }
		void read_from_path(char const *path);
		bool retrieve_next_read();
		void finish();

		std::uint64_t read_length() const { return m_read_length; }
		read_buffer_type const &read_buffer() const { return m_read_buffer; }

		template <typename t_cb>
		void iterate_characters(std::uint64_t start_pos, t_cb &&cb) const { iterate_packed_characters(m_read_buffer, m_read_length, start_pos, cb); }

		template <typename t_cb>
		void iterate_character_pairs(std::uint64_t lhs_start, std::uint64_t rhs_start, t_cb &&cb) const { iterate_packed_character_pairs(m_read_buffer, m_read_length, lhs_start, rhs_start, cb); }

	private:
		bool handle_identifier(libbio::fasta_reader_base &, std::string_view, std::span <std::string_view const>) override { return true; }
		bool handle_sequence_chunk(libbio::fasta_reader_base &reader, std::string_view sv, bool has_newline) override; // The string view does not have the newline character.
		bool handle_sequence_end(libbio::fasta_reader_base &reader) override { return false; }

		bool handle_identifier(libbio::fastq_reader_base &, std::string_view) override { return true; }
		bool handle_sequence_chunk(libbio::fastq_reader_base &reader, std::string_view sv, bool has_newline) override; // The string view does not have the newline character.
		bool handle_sequence_end(libbio::fastq_reader_base &reader) override { return true; }
		bool handle_quality_chunk(libbio::fastq_reader_base &, std::string_view, bool) override { return true; }
		bool handle_quality_end(libbio::fastq_reader_base &) override { return false; }
	};
}
