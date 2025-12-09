#pragma once

#include <cstdint>
#include <string>
#include <string_view>

#include "pack_characters.hpp"
#include "packed_character_iteration.hpp"


namespace sf {

class reader_adapter; // Fwd.


struct reader_adapter_delegate {
  virtual ~reader_adapter_delegate() {}
  /// Return true if the reader_adapter should call the error reporting member
  /// functions below.
  virtual bool should_report_errors_for_path(reader_adapter& adapter,
                                             std::string_view path) = 0;
  /// Called when the reader_adapter has found the first read with an unexpected
  /// (wildcard) character.
  virtual void found_first_read_with_unexpected_character(
      reader_adapter& adapter, std::string_view path, std::uint64_t lineno) = 0;
  /// Report the total number of reads with unexpected characters. Called after
  /// processing all the reads.
  virtual void found_total_reads_with_unexpected_characters(
      reader_adapter& adapter, std::string_view path, std::uint64_t count) = 0;
};


class reader_adapter {
 public:
  typedef packed_word_vector read_buffer_type;

 protected:
  read_buffer_type m_read_buffer;
  reader_adapter_delegate* m_delegate{};
  std::uint64_t m_read_length{};

 public:
  constexpr reader_adapter() = default;

  /// Construct with a delegate, needed to process reads.
  explicit reader_adapter(reader_adapter_delegate& delegate)
      : m_delegate(&delegate) {}

  virtual ~reader_adapter() {}

  /// Open the file at the given path.
  virtual void read_from_path(char const* path) = 0;
  void read_from_path(std::string const& path) { read_from_path(path.c_str()); }
  /// Retrieve the next read.
  /// \returns true if a read could be retrieved.
  virtual bool retrieve_next_read() = 0;
  /// Call after finishing with the file.
  virtual void finish() = 0;

  /// Length of the current read.
  std::uint64_t read_length() const { return m_read_length; }
  /// Access the current read.
  read_buffer_type const& read_buffer() const { return m_read_buffer; }

  void set_delegate(reader_adapter_delegate& delegate) {
    m_delegate = &delegate;
  }

  /// Iterate the characters of the current read.
  template <typename t_cb>
  void iterate_characters(std::uint64_t start_pos, t_cb&& cb) const {
    iterate_packed_characters(m_read_buffer, m_read_length, start_pos, cb);
  }

  /// Iterate the character pairs of the current read.
  template <typename t_cb>
  void iterate_character_pairs(std::uint64_t lhs_start, std::uint64_t rhs_start,
                               t_cb&& cb) const {
    iterate_packed_character_pairs(m_read_buffer, m_read_length, lhs_start,
                                   rhs_start, cb);
  }
};


// Helper for calling finish() automatically.
template <typename t_reader_adapter>
struct reader_adapter_guard {
  t_reader_adapter& adapter;

  ~reader_adapter_guard() { adapter.finish(); }
};
}  // namespace sf
