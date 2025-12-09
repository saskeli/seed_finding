#pragma once

#include <bit>
#include <cstddef>
#include <cstdint>
#include <libbio/assert.hh>

#include "pack_characters.hpp"


namespace sf::detail {

// The purpose of this class is to maintain the state while iterating
// characters or pairs of characters in a packed string.
class packed_character_iteration_context {
 private:
  std::uint64_t const m_shift_amt{};
  std::uint64_t const m_mask{};
  std::uint64_t const m_first_word_index{};
  std::uint64_t const m_last_word_index{};
  std::uint64_t const m_word_count{};

  std::uint64_t m_current_word{};
  std::uint64_t m_next_word{};

 public:
  packed_character_iteration_context(packed_word_vector const& read_buffer,
                                     std::uint64_t start_pos,
                                     std::uint64_t limit)
      : m_shift_amt{2U * (start_pos % 32U)},
        m_mask{~(UINT64_C(0xFFFF'FFFF'FFFF'FFFF) << m_shift_amt)},
        m_first_word_index(start_pos / 32U),
        m_last_word_index((limit - 1U) / 32U),
        m_word_count((limit - start_pos + 31U) / 32U),
        m_next_word{read_buffer[m_first_word_index] << m_shift_amt} {}

  inline void update(std::uint64_t ww);
  void rotate_current_word() { m_current_word = std::rotl(m_current_word, 2); }

  std::uint64_t word_count() const { return m_word_count; }
  std::uint64_t first_word_index() const { return m_first_word_index; }
  std::uint64_t last_word_index() const { return m_last_word_index; }
  std::uint64_t current_word() const { return m_current_word; }
};


void packed_character_iteration_context::update(std::uint64_t word_) {
  m_current_word = m_next_word;
  m_next_word = std::rotl(word_, m_shift_amt);
  m_current_word |= m_next_word & m_mask;
  m_next_word &= ~m_mask;
}
}  // namespace sf::detail


namespace sf {

// To iterate the characters, we have cases like the following:
// – One word encloses the character range.
// – j words enclose the character range but there are less than (j - 2) * 32
// characters in total. – j words enclose the character range and there are at
// least (j - 1) * 32 characters in total.
template <typename t_cb>
inline void iterate_packed_characters(packed_word_vector const& src,
                                      std::uint64_t length,
                                      std::uint64_t const start_pos,
                                      t_cb&& cb) {
  if (length <= start_pos) return;

  detail::packed_character_iteration_context ctx{src, start_pos, length};
  std::uint8_t const remaining_characters((length - start_pos) % 32U ?: 32U);

  for (std::uint64_t ii{1}; ii < ctx.word_count(); ++ii) {
    ctx.update(src[ctx.first_word_index() + ii]);
    for (std::uint8_t jj{}; jj < 32U; ++jj) {
      ctx.rotate_current_word();
      cb(ctx.current_word() & 0x3);
    }
  }

  // It does not matter if we reload the last word since we only handle
  // the first remaining_characters characters.
  ctx.update(src[ctx.last_word_index()]);
  for (std::uint8_t jj{}; jj < remaining_characters; ++jj) {
    ctx.rotate_current_word();
    cb(ctx.current_word() & 0x3);
  }
}


template <typename t_cb>
inline void iterate_packed_character_pairs(packed_word_vector const& src,
                                           std::uint64_t length,
                                           std::uint64_t const lhs_start_pos,
                                           std::uint64_t const rhs_start_pos,
                                           t_cb&& cb) {
  typedef detail::packed_character_iteration_context context_type;

  if (length <= rhs_start_pos) return;

  libbio_assert_lte(lhs_start_pos, rhs_start_pos);

  context_type lhs_ctx{src, lhs_start_pos, length};
  context_type rhs_ctx{src, rhs_start_pos, length};
  std::uint8_t const rhs_remaining_characters((length - rhs_start_pos) % 32U
                                                  ?: 32U);

  // Could not come up with a way to eliminate this branch.
  auto const lhs_last_word_idx{lhs_ctx.word_count() == rhs_ctx.word_count()
                                   ? rhs_ctx.last_word_index()
                                   : lhs_ctx.first_word_index() +
                                         rhs_ctx.word_count()};

  for (std::uint64_t ii{1}; ii < rhs_ctx.word_count(); ++ii) {
    lhs_ctx.update(src[lhs_ctx.first_word_index() + ii]);
    rhs_ctx.update(src[rhs_ctx.first_word_index() + ii]);
    for (std::uint8_t jj{}; jj < 32U; ++jj) {
      lhs_ctx.rotate_current_word();
      rhs_ctx.rotate_current_word();
      cb(lhs_ctx.current_word() & 0x3, rhs_ctx.current_word() & 0x3);
    }
  }

  // It does not matter if we reload the last word since we only handle
  // the first remaining_characters characters.
  lhs_ctx.update(src[lhs_last_word_idx]);
  rhs_ctx.update(src[rhs_ctx.last_word_index()]);
  for (std::uint8_t jj{}; jj < rhs_remaining_characters; ++jj) {
    lhs_ctx.rotate_current_word();
    rhs_ctx.rotate_current_word();
    cb(lhs_ctx.current_word() & 0x3, rhs_ctx.current_word() & 0x3);
  }
}
}  // namespace sf
