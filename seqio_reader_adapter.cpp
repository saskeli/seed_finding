#include <cassert>
#include <cstdint>
#include <memory>
#include <SeqIO/SeqIO.hh>
#include "pack_characters.hpp"
#include "seqio_reader_adapter.hpp"


namespace sf {

	void seqio_reader_adapter::read_from_path(char const *path)
	{
		m_reader = std::make_unique <seq_io::Reader_x>(path); // Moving not possible.
		m_reader->enable_reverse_complements();
		m_input_path = path;
		m_should_report_errors = m_delegate->should_report_errors_for_path(*this, m_input_path);
	}


	bool seqio_reader_adapter::retrieve_next_read()
	{
		while (true)
		{
			auto const len{m_reader->get_next_read_to_buffer()};
			assert(0 <= len);
			m_read_length = len;

			if (0 == read_length())
				return false;

			m_read_buffer.clear();
			auto const res{pack_characters(m_reader->read_buf, m_read_buffer, 0)};

			if (UINT64_MAX != res)
				return true;

			// Reached if the current read is invalid. We report only the first error.
			if (0 == m_skipped_reads && m_should_report_errors)
				m_delegate->found_first_read_with_unexpected_character(*this, m_input_path, 0); // No line number information available.
			++m_skipped_reads;
		}
	}


	void seqio_reader_adapter::finish()
	{
		if (m_skipped_reads && m_should_report_errors)
			m_delegate->found_total_reads_with_unexpected_characters(*this, m_input_path, m_skipped_reads);

		m_reader.reset();
	}
}
