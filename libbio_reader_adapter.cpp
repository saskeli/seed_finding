#include <algorithm>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <iostream>
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

	constexpr inline unsigned char reverse_bits_32b(unsigned char cc)
	{
		// From https://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith32Bits
		return ((cc * UINT32_C(0x0802) & UINT32_C(0x22110)) | (cc * UINT32_C(0x8020) & UINT32_C(0x88440))) * UINT32_C(0x10101) >> 16U;
	}


	constexpr inline unsigned char reverse_bits_64b(unsigned char cc)
	{
		// From https://graphics.stanford.edu/~seander/bithacks.html#ReverseByteWith32Bits
		return  ((cc * UINT64_C(0x80200802)) & UINT64_C(0x0884422110)) * UINT64_C(0x0101010101) >> 32U;
	}


	void reverse_complement_packed_scalar_64b(std::span <std::uint64_t> packed_input, std::uint64_t const length)
	{
		if (0 == length)
			return;

		// Reverse the 8-byte words.
		std::reverse(packed_input.begin(), packed_input.end() - 1);

		auto const reverse_complement([](std::uint64_t word){
			word = std::byteswap(word);

			auto *bytes{reinterpret_cast <char *>(&word)};
			bytes[0] = reverse_bits_64b(bytes[0]);
			bytes[1] = reverse_bits_64b(bytes[1]);
			bytes[2] = reverse_bits_64b(bytes[2]);
			bytes[3] = reverse_bits_64b(bytes[3]);

			word = ~word;
			auto const b1{(word << 1U) & UINT64_C(0xAAAA'AAAA'AAAA'AAAA)};
			auto const b2{(word >> 1U) & UINT64_C(0x5555'5555'5555'5555)};
			word = b1 | b2;

			return word;
		});

		auto const shift_right_amt{(length % 32U) * 2U};
		auto const mask{~(UINT64_C(0xFFFF'FFFF'FFFF'FFFF) >> shift_right_amt)}; // Avoid UB by replacing shift with rotate + bitwise and.
		packed_input.front() = reverse_complement(packed_input.front());
		packed_input.front() = std::rotr(packed_input.front(), shift_right_amt);
		for (std::size_t ii{1U}; ii < packed_input.size() - 1; ++ii)
		{
			auto &word(packed_input[ii - 1]);
			auto &next_word(packed_input[ii]);
			next_word = reverse_complement(next_word);
			word |= next_word >> shift_right_amt;
			next_word = std::rotr(next_word, shift_right_amt) & mask;
		}
	}
}


namespace sf {

	void libbio_reader_adapter::read_from_path(char const *path)
	{
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

				auto const status{m_reader->parse_(*m_reading_handle)};
				if (lb::sequence_reader::parsing_status::failure == status || 0 == m_read_length)
				{
					m_reading_handle->finish();
					return false;
				}

				if (m_read_is_valid)
					break;

				std::cerr << "Skipping read on line "  << m_reader->line_number() << " due to unexpected character.\n";
			}
		}

		return true;
	}


	bool libbio_reader_adapter::handle_sequence_chunk(lb::fasta_reader_base &reader, std::string_view sv, bool has_newline)
	{
		m_read_length = read_into_packed_sequence(sv, m_read_buffer, m_read_length);
		return (UINT64_MAX != m_read_length);
	}


	bool libbio_reader_adapter::handle_sequence_chunk(lb::fastq_reader_base &reader, std::string_view sv, bool has_newline)
	{
		m_read_length = read_into_packed_sequence(sv, m_read_buffer, m_read_length);
		return (UINT64_MAX != m_read_length);
	}


	bool libbio_reader_adapter::handle_sequence_end(lb::fasta_reader_base &reader)
	{
		return false;
	}


	bool libbio_reader_adapter::handle_sequence_end(lb::fastq_reader_base &reader)
	{
		return false;
	}
}
