#pragma once

#include <cstdint>
#include <libbio/assert.hh>
#include <memory>
#include <SeqIO/SeqIO.hh>
#include <string>
#include "reader_adapter.hpp"


namespace sf {

	class seqio_reader_adapter final : public reader_adapter
	{
	private:
		std::unique_ptr <seq_io::Reader_x>	m_reader{};
		std::string							m_input_path{};
		std::uint64_t						m_skipped_reads{};
		bool								m_should_report_errors{};

	public:
		using reader_adapter::reader_adapter;
		using reader_adapter::read_from_path;

		void read_from_path(char const *path) override;
		bool retrieve_next_read() override;
		void finish() override;
	};
}
