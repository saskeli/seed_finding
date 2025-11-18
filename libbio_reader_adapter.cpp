#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/fasta_reader.hh>
#include <libbio/fastq_reader.hh>
#include <libbio/file_handle.hh>
#include <libbio/file_handling.hh>
#include <libbio/sequence_reader.hh>
#include <span>
#include <string_view>
#include <vector>

#include "libbio_reader_adapter.hpp"

namespace lb = libbio;


namespace {

	void print_read(std::span <std::uint64_t const> span, std::uint64_t len)
	{
		std::cerr << "** read: ";
		for (auto word : span)
		{
			for (std::uint8_t ii{}; ii < 32U; ++ii)
			{
				word = std::rotl(word, 2);
				auto const cc{word & 0x3};
				switch (cc)
				{
					case 0: std::cerr << 'A'; break;
					case 1: std::cerr << 'C'; break;
					case 2: std::cerr << 'G'; break;
					case 3: std::cerr << 'T'; break;
					default: std::abort();
				}

				--len;
				if (0 == len)
					break;
			}
		}
		std::cerr << '\n';
	}

	std::uint64_t read_into_packed_sequence(std::string_view sv, std::vector <std::uint64_t> &dst, std::uint64_t dst_pos)
	{
		auto const push_packed([](char const cc, std::uint64_t &dst_word, std::uint64_t const dst_pos_) -> bool {
			auto const chararcter_pos{dst_pos_ % 32U};
			std::uint64_t const shift_amt{62U - 2 * chararcter_pos};
			switch (cc)
			{
				case 'A':
				case 'a':
					break;

				case 'C':
				case 'c':
					dst_word |= UINT64_C(0x1) << shift_amt;
					break;

				case 'G':
				case 'g':
					dst_word |= UINT64_C(0x2) << shift_amt;
					break;

				case 'T':
				case 't':
					dst_word |= UINT64_C(0x3) << shift_amt;
					break;

				[[unlikely]] default:
					return false;
			}

			return true;
		});

		auto it{sv.cbegin()};
		auto const end{sv.cend()};

		if (dst_pos % 32U)
		{
			auto &dst_word{dst.back()};
			do
			{
				if (it == end)
					return dst_pos;

				if (!push_packed(*it, dst_word, dst_pos)) [[unlikely]]
					return UINT64_MAX;

				++dst_pos;
				++it;
			} while (dst_pos % 32U);
		}

		if (it == end)
			return dst_pos;

		// 0 == dst_pos % 32U.
		while (true)
		{
			auto &dst_word{dst.emplace_back(0)};
			do
			{
				if (!push_packed(*it, dst_word, dst_pos)) [[unlikely]]
					return UINT64_MAX;

				++dst_pos;
				++it;

				if (it == end)
					return dst_pos;
			} while (dst_pos % 32U);

			if (it == end)
				return dst_pos;
		}
	}
}


namespace sf::detail {

	constexpr inline unsigned char reverse_bits_8_32b(unsigned char cc)
	{
		// From https://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith32Bits
		return ((cc * UINT32_C(0x0802) & UINT32_C(0x22110)) | (cc * UINT32_C(0x8020) & UINT32_C(0x88440))) * UINT32_C(0x10101) >> 16U;
	}


	constexpr inline unsigned char reverse_bits_8_64b(unsigned char cc)
	{
		// From https://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith32Bits
		return ((cc * UINT64_C(0x80200802)) & UINT64_C(0x0884422110)) * UINT64_C(0x0101010101) >> 32U;
	}


	constexpr inline std::uint64_t reverse_bits_64b(std::uint64_t cc)
	{
		// First reverse the 4-byte words, then the 2-byte words, then
		// bytes, nibbles, groups of 2 bits and finally single bits.
		//
		// For x86-64, Clang 21.1.0 produces 20 instructions (excluding ret)
		// for the extern version of this function, GCC 15 a few more.
		// (The first two swaps get replaced by something apparently more
		// efficient.) For ARM, both compilers produce a single instruction,
		// rbit.

		cc = std::rotr(cc, 32);

		{
			auto const mask{UINT64_C(0x0000'FFFF'0000'FFFF)};
			cc = ((mask & cc) << 16U) | ((~mask & cc) >> 16U);
		}

		{
			auto const mask{UINT64_C(0x00FF'00FF'00FF'00FF)};
			cc = ((mask & cc) << 8U) | ((~mask & cc) >> 8U);
		}

		{
			auto const mask{UINT64_C(0x0F0F'0F0F'0F0F'0F0F)};
			cc = ((mask & cc) << 4U) | ((~mask & cc) >> 4U);
		}

		{
			auto const mask{UINT64_C(0x3333'3333'3333'3333)};
			cc = ((mask & cc) << 2U) | ((~mask & cc) >> 2U);
		}

		{
			auto const mask{UINT64_C(0x5555'5555'5555'5555)};
			cc = ((mask & cc) << 1U) | ((~mask & cc) >> 1U);
		}

		return cc;
	}


	constexpr inline std::uint64_t reverse_bits_8x8b(std::uint64_t word)
	{
		// This requires way more instructions than reverse_bits_64b.

		word = std::byteswap(word);

		auto *bytes{reinterpret_cast <char *>(&word)};
		bytes[0] = reverse_bits_8_64b(bytes[0]);
		bytes[1] = reverse_bits_8_64b(bytes[1]);
		bytes[2] = reverse_bits_8_64b(bytes[2]);
		bytes[3] = reverse_bits_8_64b(bytes[3]);
		bytes[4] = reverse_bits_8_64b(bytes[4]);
		bytes[5] = reverse_bits_8_64b(bytes[5]);
		bytes[6] = reverse_bits_8_64b(bytes[6]);
		bytes[7] = reverse_bits_8_64b(bytes[7]);

		return word;
	}


	constexpr inline std::uint64_t reverse_complement_packed(std::uint64_t word)
	{
		word = reverse_bits_64b(word);
		word = ~word;
		auto const b1{(word << 1U) & UINT64_C(0xAAAA'AAAA'AAAA'AAAA)};
		auto const b2{(word >> 1U) & UINT64_C(0x5555'5555'5555'5555)};
		word = b1 | b2;

		return word;
	}


	void reverse_complement_packed_scalar_64b(std::span <std::uint64_t> packed_input, std::uint64_t const length)
	{
		if (0 == length)
			return;

		// Reverse the 8-byte words.
		std::reverse(packed_input.begin(), packed_input.end());

		auto const shift_right_amt{(length % 32U) * 2U};
		auto const mask{~(UINT64_C(0xFFFF'FFFF'FFFF'FFFF) >> shift_right_amt)}; // Avoid UB by replacing shift with rotate + bitwise and.
		packed_input.front() = reverse_complement_packed(packed_input.front());
		packed_input.front() = std::rotr(packed_input.front(), shift_right_amt) & mask;
		for (std::size_t ii{1U}; ii < packed_input.size(); ++ii)
		{
			auto &word(packed_input[ii - 1]);
			auto &next_word(packed_input[ii]);
			next_word = reverse_complement_packed(next_word);
			word |= next_word >> shift_right_amt;
			next_word = std::rotr(next_word, shift_right_amt) & mask;
		}
	}
}


namespace sf {

	void libbio_reader_adapter::read_from_path(char const *path)
	{
		libbio_assert(m_delegate);

		std::string_view path_{path};

		m_handle = lb::file_handle(lb::open_file_for_reading(path));
		if (path_.ends_with(".gz"))
		{
			m_gzip_handle.set_gzip_input_handle(m_handle);
			m_reading_handle = &m_gzip_handle;
			path_.remove_suffix(3);
		}
		else
		{
			m_reading_handle = &m_handle;
		}

		if (path_.ends_with(".fasta") || path_.ends_with(".fa"))
			m_reader = &m_fasta_reader;
		else if (path_.ends_with(".fastq") || path_.ends_with(".fq"))
			m_reader = &m_fastq_reader;
		else
		{
			std::cerr << "FATAL: Unexpected suffix in path “" << path << "”.\n";
			std::exit(1);
		}

		m_current_input_path = path_;
		m_parsing_block_size = m_reading_handle->io_op_blocksize();
		m_skipped_reads = 0;
		m_should_report_errors = m_delegate->should_report_errors_for_path(*this, path_);
		m_reader->prepare();
	}


	bool libbio_reader_adapter::retrieve_next_read()
	{
		if (m_next_read_is_reverse_complement)
		{
			m_next_read_is_reverse_complement = false;
			detail::reverse_complement_packed_scalar_64b(std::span{m_read_buffer.data(), m_read_buffer.size()}, m_read_length);
		}
		else
		{
			m_next_read_is_reverse_complement = true;

			while (true)
			{
				m_read_buffer.clear();
				m_read_length = 0;
				m_read_is_valid = true;

				auto const status{m_reader->parse_(*m_reading_handle, m_parsing_block_size)};
				if (lb::sequence_reader::parsing_status::failure == status || (0 == m_read_length && m_read_is_valid))
				{
					m_reading_handle->finish();
					return false;
				}

				if (m_read_is_valid)
					break;

				// Reached if the current read is invalid. We report only the first error.
				if (0 == m_skipped_reads && m_should_report_errors)
					std::cerr << "WARNING: Skipping reads with unexpected characters in " << m_current_input_path << "; first one on line " << m_reader->line_number() << ".\n";
				++m_skipped_reads;
			}
		}

		//print_read(std::span{m_read_buffer.data(), m_read_buffer.size()}, m_read_length);
		return true;
	}


	void libbio_reader_adapter::finish()
	{
		if (m_skipped_reads && m_should_report_errors)
			m_delegate->found_total_reads_with_unexpected_characters(*this, m_current_input_path, m_skipped_reads);
	}


	bool libbio_reader_adapter::handle_sequence_chunk(lb::fasta_reader_base &reader, std::string_view sv, bool has_newline)
	{
		if (m_read_is_valid) [[likely]]
		{
			auto const chunk_length{read_into_packed_sequence(sv, m_read_buffer, m_read_length)};
			m_read_is_valid = (UINT64_MAX != chunk_length);
			if (m_read_is_valid) [[likely]]
				m_read_length += chunk_length;
		}

		return true;
	}


	bool libbio_reader_adapter::handle_sequence_chunk(lb::fastq_reader_base &reader, std::string_view sv, bool has_newline)
	{
		if (m_read_is_valid) [[likely]]
		{
			auto const chunk_length{read_into_packed_sequence(sv, m_read_buffer, m_read_length)};
			m_read_is_valid = (UINT64_MAX != chunk_length);
			if (m_read_is_valid) [[likely]]
				m_read_length += chunk_length;
		}

		return true;
	}
}
