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

#include "reader_adapter.hpp"


namespace sf {

class libbio_reader_adapter final : public reader_adapter,
                                    public libbio::fasta_reader_delegate,
                                    libbio::fastq_reader_delegate {
 private:
  libbio::fasta_reader m_fasta_reader;
  libbio::fastq_reader m_fastq_reader;

  libbio::file_handle m_handle;
  libbio::gzip_reading_handle m_gzip_handle;
  std::string m_input_path;

  libbio::sequence_reader* m_reader{};
  libbio::reading_handle* m_reading_handle{};

  std::size_t m_parsing_block_size{};
  std::uint64_t m_skipped_reads{};
  bool m_should_report_errors{};
  bool m_read_is_valid{};
  bool m_next_read_is_reverse_complement{};

 public:
  constexpr libbio_reader_adapter() = default;

  explicit libbio_reader_adapter(reader_adapter_delegate& delegate)
      : reader_adapter(delegate) {
    m_fasta_reader.set_delegate(*this);
    m_fastq_reader.set_delegate(*this);
    m_handle.prepare();
    m_gzip_handle.prepare();
  }

  libbio_reader_adapter(libbio_reader_adapter const&) = delete;
  libbio_reader_adapter(libbio_reader_adapter&&) = default;
  libbio_reader_adapter& operator=(libbio_reader_adapter const&) = delete;
  libbio_reader_adapter& operator=(libbio_reader_adapter&&) & = default;

  using reader_adapter::read_from_path;

  void read_from_path(char const* path) override;
  bool retrieve_next_read() override;
  void finish() override;

 private:
  bool handle_identifier(libbio::fasta_reader_base&, std::string_view,
                         std::span<std::string_view const>) override {
    return true;
  }
  bool handle_sequence_chunk(libbio::fasta_reader_base& reader,
                             std::string_view sv, bool has_newline)
      override; // The string view does not have the newline character.
  bool handle_sequence_end(libbio::fasta_reader_base& reader) override {
    return false;
  }

  bool handle_identifier(libbio::fastq_reader_base&,
                         std::string_view) override {
    return true;
  }
  bool handle_sequence_chunk(libbio::fastq_reader_base& reader,
                             std::string_view sv, bool has_newline)
      override; // The string view does not have the newline character.
  bool handle_sequence_end(libbio::fastq_reader_base& reader) override {
    return true;
  }
  bool handle_quality_chunk(libbio::fastq_reader_base&, std::string_view,
                            bool) override {
    return true;
  }
  bool handle_quality_end(libbio::fastq_reader_base&) override { return false; }
};
}  // namespace sf
