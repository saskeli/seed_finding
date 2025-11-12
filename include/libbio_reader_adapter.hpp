#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>
#include <libbio/assert.hh>
#include <libbio/fasta_reader.hh>
#include <libbio/fastq_reader.hh>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/gzip_read_handle.hh>
#include <span>
#include <string>
#include <string_view>
#include <vector>
#include "libbio/sequence_reader.hh"

namespace sf::detail {

	void reverse_complement_packed_scalar_64b(std::span <std::uint64_t> packed_input, std::uint64_t const length);
}


namespace sf {

	class libbio_reader_adapter final : public libbio::fasta_reader_delegate, libbio::fastq_reader_delegate
	{
	public:
		typedef std::vector <std::uint64_t> read_buffer_type;

	private:
		struct iteration_context
		{
			std::uint64_t const start_idx{};
			std::uint64_t const shift_amt{};
			std::uint64_t const mask{};

			std::uint64_t const word_count{};
			std::uint64_t const last_word_idx{};

			std::uint64_t current_word{};
			std::uint64_t next_word{};

		public:
			iteration_context(
				read_buffer_type const &read_buffer,
				std::uint64_t read_length,
				std::uint64_t start_pos
			):
				start_idx{start_pos / 32U},
				shift_amt{2U * (start_pos % 32U)},
				mask{~(UINT64_C(0xFFFF'FFFF'FFFF'FFFF) << shift_amt)},
				word_count{(read_length - start_pos + 31U) / 32U},
				last_word_idx{(read_length + 31U) / 32U - start_pos / 32U},
				next_word{read_buffer[start_idx]}
			{
			}

			iteration_context(
				read_buffer_type const &read_buffer,
				std::uint64_t read_length,
				std::uint64_t start_pos,
				std::uint64_t limit
			):
				start_idx{start_pos / 32U},
				shift_amt{2U * (start_pos % 32U)},
				mask{~(UINT64_C(0xFFFF'FFFF'FFFF'FFFF) << shift_amt)},
				word_count{(limit - start_pos + 31U) / 32U},
				last_word_idx{(read_length + 31U) / 32U - (limit - 1) / 32U},
				next_word{read_buffer[start_idx]}
			{
			}

			inline void update(std::uint64_t ww);
			void rotate_current_word() { current_word = std::rotl(current_word, 2); }
		};

	private:
		libbio::fasta_reader		m_fasta_reader;
		libbio::fastq_reader		m_fastq_reader;

		libbio::file_handle			m_handle;
		libbio::gzip_reading_handle	m_gzip_handle;

		read_buffer_type			m_read_buffer;
		std::uint64_t				m_read_length{};

		libbio::sequence_reader		*m_reader{};
		libbio::reading_handle		*m_reading_handle{};

		bool						m_read_is_valid{};
		bool						m_next_read_is_reverse_complement{};

	public:
		libbio_reader_adapter()
		{
			m_fasta_reader.set_delegate(*this);
			m_fastq_reader.set_delegate(*this);
			m_handle.prepare();
			m_gzip_handle.prepare();
		}

		// We pass this to m_fasta_reader and m_fastq_reader, so moving and copying are disabled for now.
		libbio_reader_adapter(libbio_reader_adapter const &) = delete;
		libbio_reader_adapter(libbio_reader_adapter &&) = delete;
		libbio_reader_adapter &operator=(libbio_reader_adapter const &) = delete;
		libbio_reader_adapter &operator=(libbio_reader_adapter &&) = delete;

		void read_from_path(std::string const &path) { read_from_path(path.c_str()); }
		void read_from_path(char const *path);
		bool retrieve_next_read();

		std::uint64_t read_length() const { return m_read_length; }
		read_buffer_type const &read_buffer() const { return m_read_buffer; }

		template <typename t_cb>
		void iterate_characters(std::uint64_t start_pos, t_cb &&cb) const;

		template <typename t_cb>
		void iterate_character_pairs(std::uint64_t start_lhs, std::uint64_t const limit_lhs, std::uint64_t start_rhs, t_cb &&cb) const;

	private:
		bool handle_identifier(libbio::fasta_reader_base &, std::string_view, std::span <std::string_view const>) override { return true; }
		bool handle_sequence_chunk(libbio::fasta_reader_base &reader, std::string_view sv, bool has_newline) override; // The string view does not have the newline character.
		bool handle_sequence_end(libbio::fasta_reader_base &reader) override;

		bool handle_identifier(libbio::fastq_reader_base &, std::string_view) override { return true; }
		bool handle_sequence_chunk(libbio::fastq_reader_base &reader, std::string_view sv, bool has_newline) override; // The string view does not have the newline character.
		bool handle_sequence_end(libbio::fastq_reader_base &reader) override;
		bool handle_quality_chunk(libbio::fastq_reader_base &, std::string_view, bool) override { return true; }
		bool handle_quality_end(libbio::fastq_reader_base &) override { return true; }
	};


	void libbio_reader_adapter::iteration_context::update(std::uint64_t word_)
	{
		current_word = next_word << shift_amt;
		next_word = std::rotl(word_, shift_amt);
		current_word |= next_word & mask;
	}


	template <typename t_cb>
	void libbio_reader_adapter::iterate_characters(std::uint64_t const start_pos, t_cb &&cb) const
	{
		if (m_read_length <= start_pos)
			return;

		iteration_context ctx{m_read_buffer, m_read_length, start_pos};
		std::uint8_t const remaining_characters((m_read_length - start_pos) % 32U ?: 32U);

		for (std::uint64_t ii{1}; ii < ctx.word_count; ++ii)
		{
			ctx.update(m_read_buffer[ctx.start_idx + ii]);
			for (std::uint8_t jj{}; jj < 32U; ++jj)
			{
				ctx.rotate_current_word();
				cb(ctx.current_word & 0x3);
			}
		}

		ctx.update(m_read_buffer[ctx.last_word_idx]); // It does not matter if we reload the last word.
		for (std::uint8_t jj{}; jj < remaining_characters; ++jj)
		{
			ctx.rotate_current_word();
			cb(ctx.current_word & 0x3);
		}
	}


	template <typename t_cb>
	void libbio_reader_adapter::iterate_character_pairs(std::uint64_t const lhs_start_pos, std::uint64_t const lhs_limit, std::uint64_t const rhs_start_pos, t_cb &&cb) const
	{
		libbio_assert_lte(lhs_start_pos, rhs_start_pos);

		if (m_read_length <= rhs_start_pos)
			return;

		iteration_context lhs_ctx{m_read_buffer, m_read_length, lhs_start_pos, lhs_limit};
		iteration_context rhs_ctx{m_read_buffer, m_read_length, rhs_start_pos};
		std::uint8_t const rhs_remaining_characters((m_read_length - rhs_start_pos) % 32U ?: 32U);
		for (std::uint64_t ii{1}; ii < rhs_ctx.word_count; ++ii)
		{
			lhs_ctx.update(m_read_buffer[lhs_ctx.start_idx + ii]);
			rhs_ctx.update(m_read_buffer[rhs_ctx.start_idx + ii]);
			for (std::uint8_t jj{}; jj < 32U; ++jj)
			{
				lhs_ctx.rotate_current_word();
				rhs_ctx.rotate_current_word();
				cb(lhs_ctx.current_word & 0x3, rhs_ctx.current_word & 0x3);
			}
		}

		lhs_ctx.update(m_read_buffer[lhs_ctx.last_word_idx]);
		rhs_ctx.update(m_read_buffer[rhs_ctx.last_word_idx]);
		for (std::uint8_t jj{}; jj < rhs_remaining_characters; ++jj)
		{
			lhs_ctx.rotate_current_word();
			rhs_ctx.rotate_current_word();
			cb(lhs_ctx.current_word & 0x3, rhs_ctx.current_word & 0x3);
		}
	}
}
