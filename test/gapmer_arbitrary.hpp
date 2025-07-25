/*
 * Copyright (c) 2025 Tuukka Norri
 * This code is licensed under MIT license (see LICENSE for details).
 */

#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <ostream>
#include <set>
#include <stdexcept>
#include <vector>

#include "../include/gapmer.hpp"
#include "../include/string_buffer.hpp"
#include "../include/util.hpp"
#include "nucleotide.hpp"
#include "test.hpp"


// FIXME: Some of the types used in the tests need to have typedefs for the
// associated types since the preprocessor has difficulties with commas as part
// of the type names. Some of the typedefs have very long names in order to
// encode the name of the test they have to do with. Separating the tests into
// multiple translation units could help with this.

namespace {

constexpr static inline uint8_t const GAPMER_MAX_SUFFIX_LENGTH_DEFAULT{
    UINT8_MAX - 1};
typedef sf::gapmer<> default_gapmer_type;


// The purpose of the following is to make it nicer to specify constraints for
// producing gapmer_data.
template <typename>
struct pointer_to_data_member_traits {};

template <typename t_data_member, typename t_class>
struct pointer_to_data_member_traits<t_data_member t_class::*> {
  typedef t_class class_type;
  typedef t_data_member data_member_type;
};

template <typename t_data_member_ptr>
using pointer_to_data_member_t =
    pointer_to_data_member_traits<t_data_member_ptr>::data_member_type;


struct gapmer_data_config_base {
  uint8_t min_length{};
  uint8_t max_length{sf::gapmer<>::max_k};
  uint8_t min_suffix_length{};
  uint8_t max_suffix_length{GAPMER_MAX_SUFFIX_LENGTH_DEFAULT};
};


template <typename t_type, t_type gapmer_data_config_base::*t_member>
struct gapmer_data_constraint_tpl {
  t_type value{};

  constexpr gapmer_data_constraint_tpl(t_type value_) : value{value_} {};
  constexpr void assign(gapmer_data_config_base &config) const {
    config.*t_member = value;
  }
};

template <auto t_data_member_ptr>
using gapmer_data_constraint_t = gapmer_data_constraint_tpl<
    pointer_to_data_member_t<decltype(t_data_member_ptr)>, t_data_member_ptr>;

typedef gapmer_data_constraint_t<&gapmer_data_config_base::min_length>
    gapmer_min_length_constraint;
typedef gapmer_data_constraint_t<&gapmer_data_config_base::max_length>
    gapmer_max_length_constraint;
typedef gapmer_data_constraint_t<&gapmer_data_config_base::min_suffix_length>
    gapmer_min_suffix_length_constraint;
typedef gapmer_data_constraint_t<&gapmer_data_config_base::max_suffix_length>
    gapmer_max_suffix_length_constraint;


template <auto... t_constraint>  // For now this requires C++20 (I think).
struct gapmer_data_config : gapmer_data_config_base {
  constexpr gapmer_data_config() {
    (t_constraint.assign(*this), ...);

    if (!(0 == min_suffix_length || min_suffix_length < min_length))
      throw std::logic_error{"Invalid min suffix length."};

    if (!(GAPMER_MAX_SUFFIX_LENGTH_DEFAULT == max_suffix_length ||
          ((0 == max_suffix_length || max_suffix_length < min_length) &&
           min_suffix_length <= max_suffix_length)))
      throw std::logic_error{"Invalid max suffix length."};
  }
};


template <typename t_gapmer = default_gapmer_type>
struct gapmer_data {
  // FIXME: Test different instantiations.
  typedef t_gapmer gapmer_type;

  uint64_t sequence{};  // Half-nibbles in reverse order like in gapmer.
  uint8_t length{};
  uint8_t suffix_length{};
  uint8_t gap_length{};

  gapmer_type to_gapmer() const;
  operator gapmer_type() const { return to_gapmer(); }
  void write_to_buffer(std::vector<uint64_t> &buffer) const;
  uint8_t gap_start() const {
    return suffix_length ? length - suffix_length : 0;
  }
  uint8_t suffix_start() const {
    return suffix_length ? gap_start() + gap_length : 0;
  }
  uint8_t prefix_length() const {
    return suffix_length ? length - suffix_length : length;
  }
};


// Parametrised for rc::Arbitrary.
template <typename t_gapmer = default_gapmer_type, auto... t_constraint>
struct gapmer_data_ : public gapmer_data<t_gapmer> {};


template <typename t_gapmer = default_gapmer_type>
struct gapmer_pair {
  typedef gapmer_data<t_gapmer>::gapmer_type gapmer_type;
  gapmer_data<gapmer_type> source{};
  gapmer_data<gapmer_type> target{};
};


template <typename t_gapmer = default_gapmer_type>
struct aligning_gapmer_pair : public gapmer_pair<t_gapmer> {};

// Tried with a template with an enum parameter, which resulted in
// less-than-useful log messages; hence we have the following class hierarchy.
template <typename t_gapmer = default_gapmer_type>
struct non_aligning_gapmer_pair : public gapmer_pair<t_gapmer> {};

template <typename t_gapmer = default_gapmer_type>
struct non_aligning_gapmer_pair_value_mismatch
    : public non_aligning_gapmer_pair<t_gapmer> {};

template <typename t_gapmer = default_gapmer_type>
struct non_aligning_gapmer_pair_source_length_greater
    : public non_aligning_gapmer_pair<t_gapmer> {};

template <typename t_gapmer = default_gapmer_type>
struct non_aligning_gapmer_pair_gap_length_mismatch
    : public non_aligning_gapmer_pair<t_gapmer> {};

template <typename t_gapmer = default_gapmer_type>
struct non_aligning_gapmer_pair_gap_position_mismatch
    : public non_aligning_gapmer_pair<t_gapmer> {};


template <typename t_gapmer = default_gapmer_type>
struct huddinge_neighbourhood {
  typedef t_gapmer gapmer_type;

  struct cmp {
    bool operator()(gapmer_type const lhs, gapmer_type const rhs) const {
      return lhs.data() < rhs.data();
    }
  };

  gapmer_data<gapmer_type> gd{};

  std::set<gapmer_type, cmp> values;
};


// Parametrised for rc::Arbitrary.
template <typename t_gapmer_data>
struct huddinge_neighbourhood_
    : public huddinge_neighbourhood<typename t_gapmer_data::gapmer_type> {
  typedef t_gapmer_data gapmer_data_type;
};


template <typename t_gapmer>
auto gapmer_data<t_gapmer>::to_gapmer() const -> gapmer_type {
  auto const retval(gapmer_type(sequence, length, gap_start(), gap_length));
  if (!retval.is_valid()) {
    std::cerr << "sequence: " << sequence << " length: " << +length
              << " suffix_length: " << +suffix_length
              << " gap_length: " << +gap_length << '\n';
    RC_FAIL("Invalid gapmer");
  }
  return retval;
}


template <typename t_gapmer>
void gapmer_data<t_gapmer>::write_to_buffer(
    std::vector<uint64_t> &buffer) const {
  if (!length) return;

  buffer.resize((length + gap_length - 1) / 8 + 1, 0);
  auto *dst(reinterpret_cast<char *>(buffer.data()));
  auto sequence_{std::rotr(sequence, 2 * (length - 1))};
  uint8_t i{};

  auto const output_to_limit([&](uint8_t const limit) {
    while (i < limit) {
      dst[i] = sf::v_to_nuc[sequence_ & 0x3];
      sequence_ = std::rotl(sequence_, 2);
      ++i;
    }
  });

  output_to_limit(gap_start());

  for (uint8_t j{}; j < gap_length; ++j) {
    dst[i] = '.';
    ++i;
  }

  output_to_limit(length + gap_length);
}


template <typename t_gapmer>
std::ostream &operator<<(std::ostream &os, gapmer_data<t_gapmer> const gd) {
  os << gd.to_gapmer().to_string() << " (" << +gd.length << ", "
     << +gd.gap_start() << ", " << +gd.gap_length << ')';
  return os;
}


template <typename t_gapmer>
std::ostream &operator<<(std::ostream &os, gapmer_pair<t_gapmer> const &pp) {
  os << "source: " << pp.source << " target: " << pp.target;
  return os;
}


template <typename t_gapmer>
std::ostream &operator<<(std::ostream &os,
                         huddinge_neighbourhood<t_gapmer> const &hn) {
  os << "source: '" << hn.gd << "' values.size(): " << hn.values.size();
  return os;
}
}  // namespace


namespace rc {

template <typename t_gapmer>
struct Arbitrary<gapmer_data<t_gapmer>> {
  static Gen<gapmer_data<t_gapmer>> arbitrary() {
    return gen::just(gapmer_data<t_gapmer>{});
  }
};


template <typename t_gapmer, auto... t_constraint>
struct Arbitrary<gapmer_data_<t_gapmer, t_constraint...>> {
  static Gen<gapmer_data_<t_gapmer, t_constraint...>> arbitrary() {
    // We first determine the k-mer length. If it is at least two, we allow a
    // non-zero suffix length. If we got a non-empty suffix, we determine a
    // non-zero gap length.

    typedef gapmer_data_<t_gapmer, t_constraint...> return_type;
    typedef typename return_type::gapmer_type gapmer_type;

    constexpr static gapmer_data_config<t_constraint...> const
        config;  // static b.c. I don’t want to capture config in every lambda.

    return gen::mapcat(
        gen::inRange<uint8_t>(config.min_length, config.max_length + 1),
        [](auto const length) {
          auto suffix_length_gen(
              length ? gen::inRange(
                           +config.min_suffix_length,
                           std::min(+length, 1 + config.max_suffix_length))
                     : gen::just(0));
          return gen::mapcat(
              suffix_length_gen,
              [length](uint8_t const suffix_length) {  // Suffix length depends
                                                       // on total length.
                auto gap_length_gen(
                    suffix_length
                        ? gen::inRange<uint8_t>(1, gapmer_type::max_gap + 1)
                        : gen::just(uint8_t{}));
                return gen::mapcat(
                    gap_length_gen,
                    [length, suffix_length](
                        auto const gap_length) {  // Gap length depends on
                                                  // suffix length.
                      auto const value_limit(uint64_t(1) << (2 * length));
                      return gen::map(gen::inRange(uint64_t{}, value_limit),
                                      [length, suffix_length,
                                       gap_length](auto const sequence) {
                                        return_type const retval{
                                            sequence, length, suffix_length,
                                            gap_length};
                                        return retval;
                                      });
                    });
              });
        });
  }
};


template <typename t_gapmer>
struct Arbitrary<aligning_gapmer_pair<t_gapmer>> {
  static Gen<aligning_gapmer_pair<t_gapmer>> arbitrary() {
    return gen::mapcat(
        gen::arbitrary<gapmer_data_<t_gapmer>>(), [](auto const gd) {
          return gen::mapcat(
              gd.length ? gen::inRange(uint8_t{}, gd.length)
                        : gen::just(uint8_t{}),
              [gd](auto const pos) {  // Determine the source offset in target.
                return gen::map(
                    gen::inRange(0, gd.length - pos + 1),
                    [gd, pos](
                        uint8_t const length) {  // Determine the source length.
                      if (!length)
                        return aligning_gapmer_pair<t_gapmer>{
                            gapmer_data<t_gapmer>{}, gd};

                      // Copy the sequence to get an aligning pair.
                      auto const source_end(pos + length);
                      auto sequence(gd.sequence);
                      sequence >>= 2 * (gd.length - source_end);

                      // Zero the non-sequence characters.
                      {
                        uint64_t mask{};
                        mask = ~mask;
                        mask >>= 2 * (32 - length);
                        sequence &= mask;
                      }

                      // Determine the suffix length and the gap length of the
                      // source.
                      auto const target_gap_start(gd.gap_start());
                      uint8_t const suffix_length(
                          pos < target_gap_start &&
                                  target_gap_start < source_end
                              ? source_end - target_gap_start
                              : 0);
                      uint8_t const gap_length(suffix_length ? gd.gap_length
                                                             : 0);

                      return aligning_gapmer_pair<t_gapmer>{
                          gapmer_data<t_gapmer>{sequence, length, suffix_length,
                                                gap_length},
                          gd};
                    });
              });
        });
  }
};


template <typename t_gapmer>
struct Arbitrary<non_aligning_gapmer_pair_value_mismatch<t_gapmer>> {
  static Gen<non_aligning_gapmer_pair_value_mismatch<t_gapmer>> arbitrary() {
    return gen::mapcat(
        gen::arbitrary<
            gapmer_data_<t_gapmer, gapmer_min_length_constraint{1}>>(),
        [](auto const gd) {
          return gen::map(
              gen::inRange(0, 2 * gd.length),
              [gd](auto const pos)
                  -> non_aligning_gapmer_pair_value_mismatch<t_gapmer> {
                // Make a new value by flipping some bit in the value.
                // (Currently all bit patterns are valid in the range specified
                // by length.)
                auto gd_(gd);
                uint64_t mask{1};
                mask <<= pos;
                gd_.sequence ^= mask;
                RC_ASSERT(gd.sequence != gd_.sequence);
                return {gd, gd_};
              });
        });
  }
};


template <typename t_gapmer>
struct Arbitrary<non_aligning_gapmer_pair_source_length_greater<t_gapmer>> {
  static Gen<non_aligning_gapmer_pair_source_length_greater<t_gapmer>>
  arbitrary() {
    return gen::map(
        gen::arbitrary<gapmer_data_<t_gapmer>>(),
        [](auto gd)
            -> non_aligning_gapmer_pair_source_length_greater<t_gapmer> {
          // Make the source longer than the target.
          auto gd_(gd);
          if (gd.length < t_gapmer::max_k)
            ++gd.length;
          else if (2 <= gd.prefix_length()) {
            // Make sure that the resulting gapmer is valid by removing the
            // extra character.
            --gd_.length;
            gd_.sequence >>= 2;
          } else {
            RC_ASSERT(2 <= gd.suffix_length);
            --gd_.length;
            --gd_.suffix_length;

            // Remove the extra character.
            uint64_t mask{};
            mask = ~mask;
            mask >>= 64 - (2 * gd_.length);
            gd_.sequence &= mask;
          }

          return {gd, gd_};
        });
  }
};


template <typename t_gapmer>
struct Arbitrary<non_aligning_gapmer_pair_gap_length_mismatch<t_gapmer>> {
  static Gen<non_aligning_gapmer_pair_gap_length_mismatch<t_gapmer>>
  arbitrary() {
    typedef gapmer_data_<t_gapmer, gapmer_min_length_constraint{2},
                         gapmer_min_suffix_length_constraint{1}>
        gapmer_data_type;
    return gen::map(
        gen::arbitrary<gapmer_data_type>(),
        []<typename t_gapmer_data>(t_gapmer_data const gd)
            -> non_aligning_gapmer_pair_gap_length_mismatch<t_gapmer> {
          typedef typename t_gapmer_data::gapmer_type gapmer_type;
          static_assert(1 <= gapmer_type::max_gap);

          auto gd_(gd);
          if (gd_.gap_length < gapmer_type::max_gap)
            ++gd_.gap_length;
          else
            --gd_.gap_length;

          return {gd, gd_};
        });
  }
};


template <typename t_gapmer>
struct Arbitrary<non_aligning_gapmer_pair_gap_position_mismatch<t_gapmer>> {
  static Gen<non_aligning_gapmer_pair_gap_position_mismatch<t_gapmer>>
  arbitrary() {
    typedef gapmer_data_<t_gapmer, gapmer_min_length_constraint{3},
                         gapmer_min_suffix_length_constraint{1}>
        gapmer_data_type;
    return gen::map(
        gen::arbitrary<gapmer_data_type>(),
        []<typename t_gapmer_data>(t_gapmer_data const gd)
            -> non_aligning_gapmer_pair_gap_position_mismatch<t_gapmer> {
          auto gd_(gd);
          if (1 == gd_.suffix_length)
            ++gd_.suffix_length;
          else
            --gd_.suffix_length;
          return {gd, gd_};
        });
  }
};


template <typename t_gapmer_data>
struct Arbitrary<huddinge_neighbourhood_<t_gapmer_data>> {
  static Gen<huddinge_neighbourhood_<t_gapmer_data>> arbitrary() {
    return gen::map(
        gen::arbitrary<t_gapmer_data>(),
        [](auto const gd) -> huddinge_neighbourhood_<t_gapmer_data> {
          // Given a gapmer_data, we generate its Huddinge neighbourhood by
          // replacing its defined characters and the first and last middle gap
          // characters one by one. We then proceed to prepend and append one
          // character.
          //
          // One modification may be done as follows.
          // 1. By replacing one defined character with another,
          //    i.e. ACGT ~> ACCT
          // 2. By replacing a non-boundary character with a gap, i.e.
          //    (a) ACGT ~> A-GT (no initial gap run)
          //    (b) in case 4b
          // 3. By replacing a defined character on gap boundary with a gap,
          //    i.e.
          //    (a) GAT-TACA ~> GA--TACA
          //    (b) A-GT ~> GT
          // 4. By replacing a gap character on gap boundary with a defined
          //    character, i.e.
          //    (a) GAT---TACA ~> GATC--TACA
          //    (b) AC-T ~> ACGT
          //    (c) in case 5b (either end)
          //    (d) in case 3a (opposite end)
          // 5. By contracting the string by one defined character, i.e.
          //    (a) ACGT ~> CGT
          //    (b) GAT--TACA ~> GAT--TAC
          // 6. By extending the string by one character
          //    (a) in the base case
          //    (b) in case 2a
          //    (c) in case 5 (opposite end only)
          //    (d) in case 3
          // 7. By prepending (appending) a string with a defined character
          //    followed (preceded) by an arbitrary number of gaps
          //    (a) in the base case without initial gap run
          //    (b) in case 3b
          //    (c) in case 5a (this covers inserting arbitrary number of gaps
          //        to a string without an initial gap run)
          // 8. If the string has a defined character followed (preceded) by a
          //    run of gap characters, arbitrarily extending
          //    the gap run in the base case.
          //
          // We try to ensure correctness at the cost of speed.
          typedef typename huddinge_neighbourhood_<t_gapmer_data>::gapmer_type
              gapmer_type;
          typedef sf::string_buffer<uint64_t> string_buffer_type;

          constexpr static std::array const characters{'A', 'C', 'G', 'T'};
          constexpr static uint8_t const max_edge_gap_length{5};

          huddinge_neighbourhood_<t_gapmer_data> retval;

          string_buffer_type str, str_;
          auto const gg(gd.to_gapmer());
          auto const gap_start(gg.gap_start());
          auto const gap_length(gg.gap_length());
          str = gg.to_string();
          auto span(str.to_span());

          auto const add_gapmer([&retval](string_buffer_type const &sb,
                                          uint8_t const gap_start_,
                                          uint8_t const gap_length_) {
            // The constructor takes the number of defined characters as the
            // second parameter.
            gapmer_type const gg_(sb.data(), sb.size() - gap_length_,
                                  gap_start_, gap_length_);
            auto const expected(sb.to_string_view());
            auto const actual(gg_.to_string());
            RC_ASSERT(expected == actual);
            RC_ASSERT(gg_.is_valid());

            retval.values.insert(gg_);
          });

          auto const modify_([&add_gapmer](string_buffer_type &sb,
                                           uint8_t const gap_start_,
                                           uint8_t const gap_length_, auto it) {
            auto &cc(*it);
            auto const orig_cc(cc);
            for (auto const cc_ : characters) {
              cc = cc_;
              add_gapmer(sb, gap_start_, gap_length_);
            }
            cc = orig_cc;
          });

          auto const modify([&](string_buffer_type &sb,
                                uint8_t const gap_start_,
                                uint8_t const gap_length_, auto it, auto end) {
            while (it != end) {
              modify_(sb, gap_start_, gap_length_, it);
              ++it;
            }
          });

          auto const extend_both([&modify_](string_buffer_type &sb,
                                            uint8_t const gap_start_,
                                            uint8_t const gap_length_) {
            // Case 6.
            // Restores sb.

            // Tail.
            sb.append('.');
            auto span_(sb.to_span());
            modify_(sb, gap_start_, gap_length_, span_.rbegin());

            // Head.
            sb >>= 1;
            span_[0] = '.';
            modify_(sb, gap_start_ ? gap_start_ + 1 : 0, gap_length_,
                    span_.begin());

            // Restore.
            sb <<= 1;
            sb.resize(sb.size() - 1);
          });

          auto const extend_with_edge_gap_run(
              [&modify_](string_buffer_type &sb) {
                // Case 7.
                // Does not restore sb.
                auto const orig_size(sb.size());

                // Tail.
                {
                  sb.append('.');
                  for (uint8_t i{}; i < max_edge_gap_length; ++i) {
                    sb.append('.');
                    auto span_(sb.to_span());
                    modify_(sb, orig_size, i + 1, span_.rbegin());
                  }
                }

                // Head.
                sb >>= 1;

                {
                  auto span_(sb.to_span());
                  span_[0] = '.';
                }

                for (uint8_t i{}; i < max_edge_gap_length; ++i) {
                  sb.resize(orig_size + i + 2);
                  sb >>= 1;
                  auto span_(sb.to_span());
                  span_[0] = '.';
                  modify_(sb, 1, i + 1, span_.begin());
                }
              });

          auto const extend_middle_gap_both_ends(
              [&add_gapmer, &modify_](string_buffer_type &sb,
                                      uint8_t const gap_start_,
                                      uint8_t const gap_length_) {
                // Case 4c.
                switch (gap_length_) {
                  case 0:
                    break;

                  case 1: {
                    auto span_(sb.to_span());
                    modify_(sb, 0, 0, span_.begin() + gap_start_);
                    break;
                  }

                  default: {
                    auto span_(sb.to_span());
                    modify_(sb, gap_start_ + 1, gap_length_ - 1,
                            span_.begin() + gap_start_);
                    modify_(sb, gap_start_, gap_length_ - 1,
                            span_.begin() + gap_start_ + gap_length_ - 1);
                    break;
                  }
                }
              });

          auto const insert_gap_if_needed([&](bool const should_extend) {
            // Cases 2a, 4b, 6b.
            if (3 <= str.size()) {
              auto const limit(str.size() - 2);
              for (std::size_t i{}; i < limit; ++i) {
                auto const cc_(span[i + 1]);
                span[i + 1] = '.';
                add_gapmer(str, i + 1, 1);

                if (should_extend) {
                  // Case 6b.
                  str_ = str;
                  extend_both(str_, i + 1, 1);
                }

                span[i + 1] = cc_;
              }
            }
          });

          // Case 6a.
          if (gd.length < gapmer_type::max_k) {
            str_ = str;
            extend_both(str_, gap_start, gap_length);
          }

          // Case 8.
          {
            if (1 == gap_start) {
              // Head.
              str_ = str;
              for (auto i{gap_length}; i < max_edge_gap_length; ++i) {
                str_.resize(str_.size() + 1);
                str_ >>= 1;
                auto span_(str_.to_span());
                span_[0] = span_[1];
                span_[1] = '.';
                add_gapmer(str_, 1, i + 1);
              }
            }

            if (gap_start + gap_length + 1U == str.size()) {
              // Tail.
              str_ = str;
              for (auto i{gap_length}; i < max_edge_gap_length; ++i) {
                str_.append('.');
                auto span_(str_.to_span());

                using std::swap;
                swap(*span_.rbegin(), *(span_.rbegin() + 1));
                add_gapmer(str_, gap_start, i + 1);
              }
            }
          }

          // Cases 1, 2a, 6b and 7a.
          if (0 == gap_length) {
            // Case 1.
            modify(str, 0, 0, span.begin(), span.end());

            // Cases 2a and 6b.
            insert_gap_if_needed(true);

            // Case 7a.
            if (gd.length < gapmer_type::max_k) {
              str_ = str;
              extend_with_edge_gap_run(str_);
            }
          } else {
            // Case 1.
            modify(str, gap_start, gap_length, span.begin(),
                   span.begin() + gap_start);
            modify(str, gap_start, gap_length,
                   span.begin() + gap_start + gap_length, span.end());
          }

          // Cases 4a and 4b (and 2b).
          switch (gap_length) {
            case 0:
              break;

            case 1: {
              // Cases 4b and 2b.
              auto &cc(*(span.begin() + gap_start));
              auto const orig_cc(cc);
              for (auto const cc_ : characters) {
                cc = cc_;
                if (gd.length < gapmer_type::max_k) add_gapmer(str, 0, 0);
                insert_gap_if_needed(false);
                cc = orig_cc;
              }
              break;
            }

            default:
              // Case 4a.
              if (gd.length < gapmer_type::max_k) {
                modify_(str, gap_start + 1, gap_length - 1,
                        span.begin() + gap_start);
                modify_(str, gap_start, gap_length - 1,
                        span.begin() + gap_start + gap_length - 1);
              }
              break;
          }

          // Cases 3a, 3b, 4d, 6d, 7b.
          if (gap_length) {
            // Head.
            if (1 == gap_start) {
              // Case 3b.
              str_ = str;
              str_ <<= gap_start + gap_length;
              str_.resize(str.size() - gap_start - gap_length);
              add_gapmer(str_, 0, 0);

              // Case 6d.
              extend_both(str_, 0, 0);

              // Case 7b.
              extend_with_edge_gap_run(str_);
            } else {
              // Case 3a.
              auto &cc(span[gap_start - 1]);
              auto const cc_(cc);
              cc = '.';
              if (gap_length < gapmer_type::max_gap)
                add_gapmer(str, gap_start - 1, gap_length + 1);

              // Case 4d.
              if (2 <= gap_length)
                modify_(str, gap_start - 1, gap_length,
                        span.begin() + gap_start + gap_length - 1);

              // Case 6d.
              if (gap_length < gapmer_type::max_gap) {
                str_ = str;
                extend_both(str_, gap_start - 1,
                            gap_length + 1);  // Same as in 3a.
              }

              cc = cc_;
            }

            // Tail.
            if (gap_start + gap_length + 1U == str.size()) {
              // Case 3b.
              str_ = str;
              str_.resize(gap_start);
              add_gapmer(str_, 0, 0);

              // Case 6d.
              extend_both(str_, 0, 0);

              // Case 7b.
              extend_with_edge_gap_run(str_);
            } else {
              // Case 3a.
              auto &cc(span[gap_start + gap_length]);
              auto const cc_(cc);
              cc = '.';
              if (gap_length < gapmer_type::max_gap)
                add_gapmer(str, gap_start, gap_length + 1);

              // Case 4d.
              if (2 <= gap_length)
                modify_(str, gap_start + 1, gap_length,
                        span.begin() + gap_start);

              // Case 6d.
              if (gap_length < gapmer_type::max_gap) {
                str_ = str;
                extend_both(str_, gap_start, gap_length + 1);  // Same as in 3a.
              }

              cc = cc_;
            }
          }

          // Cases 4c, 5, 6c and 7c.
          if (1 < str.size()) {
            // Head.
            if (0 == gap_start || 1 < gap_start) {
              // Case 5.
              str_ = str;
              str_ <<= 1;
              str_.resize(str.size() - 1);
              add_gapmer(str_, gap_length ? gap_start - 1 : 0, gap_length);

              // Case 6c.
              str_.append('.');
              auto span_(str_.to_span());
              modify_(str_, gap_length ? gap_start - 1 : 0, gap_length,
                      span_.rbegin());
              str_.resize(str_.size() - 1);

              if (gap_start) {
                // Case 4c.
                extend_middle_gap_both_ends(str_, gap_start - 1, gap_length);
              } else {
                // Case 7c.
                extend_with_edge_gap_run(str_);
              }
            }

            // Tail.
            if (gap_start + gap_length + 1U < str.size()) {
              str_ = str;

              {
                // Case 6c.
                str_ >>= 1;
                str_[0] = '.';
                auto span_(str_.to_span());
                modify_(str_, gap_length ? gap_start + 1 : 0, gap_length,
                        span_.begin());
                str_ <<= 1;
              }

              // Case 5.
              str_.resize(str_.size() - 1);
              add_gapmer(str_, gap_start, gap_length);

              if (gap_start) {
                // Case 4c.
                extend_middle_gap_both_ends(str_, gap_start, gap_length);
              } else {
                // Case 7c.
                extend_with_edge_gap_run(str_);
              }
            }
          }

          // Remove the input gapmer from the neighbourhood since its distance
          // to itsef is zero, not one.
          retval.values.erase(gg);
          retval.gd = gd;

          return retval;
        });
  }
};
}  // namespace rc


namespace sf {

template <bool middle_gap_only, uint16_t t_max_gap>
void showValue(gapmer<middle_gap_only, t_max_gap> const gg, std::ostream &os) {
  os << gg.to_string();
}


typedef gapmer_data_<default_gapmer_type, gapmer_min_length_constraint{1}>
    gapmer_arbitrary_constructor_data_type;
SF_RC_TEST_WITH_BASE_CASE(gapmer_arbitrary, Constructors, gapmer_data,
                          gapmer_arbitrary_constructor_data_type,
                          gapmer_data<> const gd, (gd.length, gd.gap_length)) {
  std::vector<uint64_t> seq;
  gd.write_to_buffer(seq);

  gapmer const gg1(gd.sequence, gd.length, gd.gap_start(), gd.gap_length);

  // Compare gg1 to gd.
  SF_ASSERT(gd.sequence == gg1.value());
  SF_ASSERT(gd.length == gg1.length());
  SF_ASSERT(gd.gap_start() == gg1.gap_start());
  SF_ASSERT(gd.gap_length == gg1.gap_length());

  // FIXME: Currently this constructor only works for non-empty sequences.
  if (gd.length) {
    gapmer const gg2(seq.data(), gd.length, gd.gap_start(), gd.gap_length);
    SF_ASSERT(gg1 == gg2);
  }

  if (!gd.suffix_length) {
    gapmer const gg3(gd.sequence, gd.length);
    SF_ASSERT(gg1 == gg3);

    // FIXME: Currently this constructor only works for non-empty sequences.
    if (gd.length) {
      gapmer const gg4(seq.data(), gd.length);
      SF_ASSERT(gg1 == gg4);
    }
  }
}


typedef gapmer_data_<default_gapmer_type, gapmer_min_length_constraint{1}>
    gapmer_arbitrary_reverse_complement_data_type;
SF_RC_TEST_WITH_BASE_CASE(gapmer_arbitrary, ReverseComplementGapmer,
                          gapmer_data,
                          gapmer_arbitrary_reverse_complement_data_type,
                          gapmer_data<> const gd, (gd.length, gd.gap_length)) {
  typedef default_gapmer_type gapmer_type;
  gapmer_type const source{gd};
  auto const target(source.reverse_complement());
  auto const source_(source.to_string());
  auto const target_(target.to_string());

  auto const length(source_.size());
  SF_ASSERT(length == target_.size());
  for (std::size_t i{}; i < length; ++i)
    SF_ASSERT(source_[i] == complement_nt(target_[length - i - 1]));
}


RC_GTEST_PROP(gapmer_arbitrary, NextGapmer,
              (gapmer_data_<default_gapmer_type,
                            gapmer_max_suffix_length_constraint{0}> const gd,
               nucleotide const nt)) {
  SF_RC_TAG(gd.length, gd.gap_length);

  typedef default_gapmer_type gapmer_type;
  gapmer_type const source{gd};
  auto const target{source.next(nt)};
  auto source_{source.to_string()};
  auto const target_{target.to_string()};
  source_.push_back(nt);
  source_.erase(0, 1);
  RC_ASSERT(source_ == target_);
}


RC_GTEST_PROP(
    gapmer_arbitrary, NextGappedGapmer,
    (gapmer_data_<default_gapmer_type, gapmer_min_length_constraint{2},
                  gapmer_min_suffix_length_constraint{1}> const gd,
     nucleotide const nt1, nucleotide const nt2)) {
  SF_RC_TAG(gd.length, gd.gap_length);

  typedef default_gapmer_type gapmer_type;
  gapmer_type const source{gd};
  auto const target{source.next(nt1, nt2)};
  auto source_{source.to_string()};
  auto const target_{target.to_string()};
  source_.insert(gd.gap_start(), 1, nt1);
  source_.push_back(nt2);
  source_.erase(0, 1);
  source_.erase(gd.suffix_start(), 1);
  RC_ASSERT(source_ == target_);
}


typedef gapmer_data_<default_gapmer_type, gapmer_min_length_constraint{1}>
    gapmer_arbitrary_align_empty_data_type;
SF_RC_TEST_WITH_BASE_CASE(gapmer_arbitrary, AlignEmptyGapmer, gapmer_data,
                          gapmer_arbitrary_align_empty_data_type,
                          gapmer_data<> const gd, (gd.length, gd.gap_length)) {
  typedef default_gapmer_type gapmer_type;
  gapmer_type const source{};
  gapmer_type const target(gd);
  SF_ASSERT((source.aligns_to<false>(
      target)));  // FIXME: also reverse complement (if needed).
}


RC_GTEST_PROP(gapmer_arbitrary, AlignNonemptyGapmer,
              (aligning_gapmer_pair<> const &pair)) {
  SF_RC_TAG(pair.target.length, pair.target.gap_length);

  typedef aligning_gapmer_pair<>::gapmer_type gapmer_type;
  gapmer_type const source{pair.source};
  gapmer_type const target{pair.target};
  RC_ASSERT((source.aligns_to<false>(
      target)));  // FIXME: also reverse complement (if needed).
}


SF_RC_TEMPLATE_TEST(gapmer_arbitrary_aligns_to, AlignNonMatchingGapmer,
                    (TypeParam const &pair),
                    non_aligning_gapmer_pair_value_mismatch<>,
                    non_aligning_gapmer_pair_source_length_greater<>) {
  SF_RC_TAG(pair.target.length, pair.target.gap_length);
  typedef typename TypeParam::gapmer_type gapmer_type;
  gapmer_type const source{pair.source};
  gapmer_type const target{pair.target};
  RC_ASSERT((!source.template aligns_to<false>(target)));
}


SF_RC_TEMPLATE_TEST(gapmer_arbitrary_aligns_to_reversible,
                    AlignNonMatchingGapmerReversible, (TypeParam const &pair),
                    non_aligning_gapmer_pair_gap_length_mismatch<>,
                    non_aligning_gapmer_pair_gap_position_mismatch<>) {
  SF_RC_TAG(pair.target.length, pair.target.gap_length);
  typedef typename TypeParam::gapmer_type gapmer_type;
  gapmer_type const gg1{pair.source};
  gapmer_type const gg2{pair.target};
  RC_ASSERT((!gg1.template aligns_to<false>(gg2)));
  RC_ASSERT((!gg2.template aligns_to<false>(gg1)));
}


// FIXME: Test other gapmer types.
typedef sf::gapmer<true, 5> gapmer_arbitrary_huddinge_neighbourhood_gapmer_type;

typedef gapmer_data_<
    gapmer_arbitrary_huddinge_neighbourhood_gapmer_type,
    gapmer_min_length_constraint{5},
    gapmer_max_length_constraint{
        gapmer_arbitrary_huddinge_neighbourhood_gapmer_type::max_k - 1}
    // FIXME: Use a type parameter here.
    >
    gapmer_arbitrary_huddinge_neighbourhood_gapmer_data_type;

typedef huddinge_neighbourhood_<
    gapmer_arbitrary_huddinge_neighbourhood_gapmer_data_type>
    gapmer_arbitrary_huddinge_neighbourhood_test_type;

typedef huddinge_neighbourhood<
    gapmer_arbitrary_huddinge_neighbourhood_test_type::gapmer_type>
    gapmer_arbitrary_huddinge_neighbourhood_result_type;

SF_RC_TEST_PROP(gapmer_arbitrary, GenerateHuddingeNeighbourhood,
                gapmer_arbitrary_huddinge_neighbourhood_test_type,
                gapmer_arbitrary_huddinge_neighbourhood_result_type const &hn,
                (hn.gd.length, hn.gd.gap_length)) {
  typedef gapmer_arbitrary_huddinge_neighbourhood_gapmer_type
      gapmer_type;  // FIXME: Use a type parameter here.

  auto const gd(hn.gd);
  gapmer_type const gg{gd};

  std::set<gapmer_type,
           gapmer_arbitrary_huddinge_neighbourhood_result_type::cmp>
      actual;  // FIXME: Use a type parameter here.
  gg.all_gap_neighbours<false, false, false>([&](gapmer_type const gg_) {
    int out{};
    SF_EXPECT(gg_.is_valid());
    SF_EXPECT(1 == gg.huddinge_distance(gg_, out));
    SF_EXPECT(1 == gg_.huddinge_distance(gg, out));
    actual.insert(gg_);
  });

  if (hn.values != actual) {
    std::cerr << "** Unexpected results in GenerateHuddingeNeighbourhood:\n";
    std::cerr << "* Tested gapmer: " << gg.to_string()
              << " (length: " << gg.length() << ")\n";

    std::vector<gapmer_type> not_found_in_actual, extra_elements_in_actual;
    std::set_difference(hn.values.begin(), hn.values.end(), actual.begin(),
                        actual.end(), std::back_inserter(not_found_in_actual));
    std::set_difference(actual.begin(), actual.end(), hn.values.begin(),
                        hn.values.end(),
                        std::back_inserter(extra_elements_in_actual));

    auto const output([gg](auto const gg_) {
      int out{};
      std::cerr << gg_.to_string() << " (H: " << gg.huddinge_distance(gg_, out)
                << ", length: " << gg_.length() << ")\n";
    });

    if (!not_found_in_actual.empty()) {
      std::cerr << "* Elements not found in actual:\n";
      for (auto const gg_ : not_found_in_actual) output(gg_);
    }

    if (!extra_elements_in_actual.empty()) {
      std::cerr << "* Extra elements in actual:\n";
      for (auto const gg_ : extra_elements_in_actual) output(gg_);
    }

    SF_FAIL();
  }
}
}  // namespace sf
