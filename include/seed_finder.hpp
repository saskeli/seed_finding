#pragma once

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <libbio/assert.hh>
#include <libbio/syncstream.hh>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "count_base.hpp"
#include "gapmer.hpp"
#include "gapmer_count.hpp"
#include "packed_read.hpp"
#include "partial_count.hpp"
#include "util.hpp"


namespace sf::detail {
// FIXME: Consider using the same value type (e.g. double) for all counts.
template <typename t_value>
struct enrichment_result {
  double ac_test_result{};
  t_value signal_count{};
  t_value background_count{};

  // FIXME: remove?
  template <typename t_value_>
  enrichment_result<t_value_> to_enrichment_result() const {
    return enrichment_result<t_value_>(ac_test_result, t_value_(signal_count),
                                       t_value_(background_count));
  }
};

enum class enrichment_check_status {
  success,
  fold_change_test_failed,
  ac_test_failed,
  discarded_earlier,
  not_canonical
};

template <typename t_value>
struct enrichment_check_result {
  enrichment_result<t_value> result{};
  enrichment_check_status status{enrichment_check_status::success};

  operator bool() const { return status == enrichment_check_status::success; }
};


enum class extension_validation_status {
  // Extending.
  ac_test_passed,
  binomial_and_ac_tests_passed,
  target_ac_test_p_value_smaller,
  // Not extending.
  source_signal_count_at_least_quadruple,
  ac_test_failed,
  binomial_test_failed,
  source_ac_test_p_value_smaller
};


enum class extension_validation_strategy {
  comparison_to_background,
  comparison_to_background_filtering,
  signal_count_only  // FIXME: AC test not needed in this case since it uses background count?
};


struct extension_validation_result {
  double binomial_test_result{};
  extension_validation_status status;
  extension_validation_strategy strategy;

  static extension_validation_result make(
      extension_validation_strategy strategy_,
      extension_validation_status status_, double binomial_test_result_ = 0) {
    return {binomial_test_result_, status_, strategy_};
  }

  operator bool() const {
    switch (status) {
      case extension_validation_status::ac_test_passed:
      case extension_validation_status::binomial_and_ac_tests_passed:
      case extension_validation_status::target_ac_test_p_value_smaller:
        return true;

      default:
        return false;
    }
  }
};


template <typename t_value>
std::ostream& operator<<(std::ostream& os, enrichment_result<t_value> er) {
  std::print(os, "(SC: {} BC: {} AC test p-value: {})", er.signal_count,
             er.background_count, er.ac_test_result);
  return os;
}


inline std::ostream& operator<<(std::ostream& os,
                                extension_validation_status status) {
  switch (status) {
    case extension_validation_status::ac_test_passed:
      os << "AC test passed";
      break;
    case extension_validation_status::binomial_and_ac_tests_passed:
      os << "binomial and AC tests passed";
      break;
    case extension_validation_status::target_ac_test_p_value_smaller:
      os << "target’s AC test p-value smaller";
      break;
    case extension_validation_status::source_signal_count_at_least_quadruple:
      os << "source signal count at least quadruple";
      break;
    case extension_validation_status::ac_test_failed:
      os << "AC test failed";
      break;
    case extension_validation_status::binomial_test_failed:
      os << "binomial test failed";
      break;
    case extension_validation_status::source_ac_test_p_value_smaller:
      os << "source’s AC test p-value was smaller";
      break;
  }

  return os;
}


inline std::ostream& operator<<(std::ostream& os,
                                extension_validation_strategy strategy) {
  switch (strategy) {
    case extension_validation_strategy::comparison_to_background:
      os << "comparison to background";
    case extension_validation_strategy::comparison_to_background_filtering:
      os << "comparison to background with filtering";
    case extension_validation_strategy::signal_count_only:
      os << "signal count only";
  }
  return os;
}
}  // namespace sf::detail


namespace sf {

template <bool t_middle_gap_only, uint8_t t_max_gap, bool t_enable_smoothing,
          bool t_filter_mers, bool t_enable_reporting_discarded_seeds>
struct seed_finder_configuration {
  constexpr static inline bool middle_gap_only{t_middle_gap_only};
  constexpr static inline uint8_t max_gap{t_max_gap};
  constexpr static inline bool enable_smoothing{t_enable_smoothing};
  constexpr static inline bool filter_mers{t_filter_mers};
  constexpr static inline bool enable_reporting_discarded_seeds{
      t_enable_reporting_discarded_seeds};
};


template <typename t_configuration>
class seed_finder {
 public:
  constexpr static inline bool middle_gap_only{
      t_configuration::middle_gap_only};
  constexpr static inline uint8_t max_gap{t_configuration::max_gap};
  constexpr static inline bool enable_smoothing{
      t_configuration::enable_smoothing};
  constexpr static inline bool filter_mers{t_configuration::filter_mers};
  constexpr static inline bool enable_reporting_discarded_seeds{
      t_configuration::enable_reporting_discarded_seeds};

  typedef gapmer<middle_gap_only, max_gap> gapmer_type;

 private:
  typedef detail::enrichment_check_status enrichment_check_status;
  typedef detail::extension_validation_status extension_validation_status;
  typedef detail::extension_validation_strategy extension_validation_strategy;
  typedef detail::extension_validation_result extension_validation_result;

  template <typename t_value>
  using enrichment_result = detail::enrichment_result<t_value>;

  template <typename t_value>
  using enrichment_check_result = detail::enrichment_check_result<t_value>;

  // FIXME: remove.
  struct seed_meta {
    double p{};
    uint64_t sig_count{};
    uint64_t bg_count{};
  };

 public:
  struct seed {
    friend seed_finder;

    gapmer_type g{};
    double p{};
    uint64_t sig_count{};
    uint64_t bg_count{};

   private:
    static seed from_seed_meta(gapmer_type g_, seed_meta mm) {
      return {g_, mm.p, mm.sig_count, mm.bg_count};
    }

    template <typename t_value>
    static seed from_enrichment_result(gapmer_type g_,
                                       enrichment_result<t_value> er) {
      return {g_, er.ac_test_result, uint64_t(er.signal_count),
              uint64_t(er.background_count)};
    }
  };

 private:
  typedef gapmer_count<gapmer_type> gapmer_count_type;
  typedef partial_count<gapmer_type> partial_count_type;
  typedef typename gapmer_count_type::value_type gapmer_count_value_type;

  typedef std::unordered_set<gapmer_type, typename gapmer_type::hash>
      gapmer_set;

  template <typename t_value>
  using gapmer_map =
      std::unordered_map<gapmer_type, t_value, typename gapmer_type::hash>;

  typedef gapmer_map<enrichment_result<double>> enrichment_result_map;

  packed_read_vector const& signal_reads_;
  packed_read_vector const& background_reads_;
  std::vector<seed> seeds_;
  std::ostream* discarded_gapmer_reporting_ostream_{};
  uint64_t sig_size_;
  uint64_t bg_size_;
  double p_;
  double p_ext_;
  double fold_lim_;
  double signal_to_total_length_ratio_;
  double memory_limit_;
  uint8_t k_lim_;
  uint8_t lookup_k_;
  bool prune_;

  // Thread safe.
  template <typename t_lhs_info, typename t_rhs_info>
  inline void report_discarded(gapmer_type lhs, gapmer_type rhs,
                               extension_validation_result res,
                               t_lhs_info const& lhs_info,
                               t_rhs_info const& rhs_info,
                               char const* test_fn) const;

  // Thread safe.
  template <typename t_value>
  inline void report_encrichment_check_failure(
      gapmer_type discarded, enrichment_check_result<t_value> res,
      char const* caller, char const* test_fn) const;

  template <typename t_value>
  [[nodiscard]] extension_validation_result validate_extension(
      gapmer_type aa, gapmer_type bb, enrichment_result<t_value> aa_er,
      enrichment_result<t_value> bb_er) const;

  template <bool calculate_ac_test_for_bb, typename t_value>
  [[nodiscard]] bool should_filter(gapmer_type aa, gapmer_type bb,
                                   enrichment_result<t_value> aa_er,
                                   enrichment_result<t_value>& bb_er) const;

  template <typename t_value>
  [[nodiscard]] enrichment_check_result<t_value> check_enrichment(
      count_pair<t_value> count) const;

  template <typename t_critical>
  [[nodiscard]] enrichment_check_result<gapmer_count_value_type>
  check_enrichment_and_filter(gapmer_type const gg, uint64_t const offset,
                              gapmer_count_type& counts,
                              t_critical&& critical) const;

  template <bool t_should_extend, typename t_critical,
            typename t_prev_critical = std::nullptr_t>
  void filter_huddinge_neighbourhood(
      gapmer_type const gg, uint64_t const offset,
      enrichment_result<gapmer_count_value_type> const& enrichment_res,
      gapmer_count_type& counts, t_critical&& critical,
      gapmer_count_type* const prev_counts = nullptr,
      t_prev_critical&& prev_critical = nullptr) const;

  enrichment_check_result<gapmer_count_value_type> check_count(
      gapmer_type const gg, uint64_t offset, gapmer_count_type& counts_k,
      gapmer_count_type& counts_k1);

  enrichment_check_result<gapmer_count_value_type> filter_count(
      gapmer_type const gg, uint64_t offset, gapmer_count_type& counts_k,
      enrichment_result_map& mm) const;

  void count_short_gapmers(gapmer_count_type& counts);

  void extend_counted(enrichment_result_map const& aa,
                      enrichment_result_map& bb,
                      partial_count_type const& p_counter) const;

  void extend(enrichment_result_map& aa, enrichment_result_map& bb,
              partial_count_type& p_counter, uint16_t k, bool prune) const;

  void filter(enrichment_result_map& mm) const;

 public:
  seed_finder(packed_read_vector const& signal_reads,
              packed_read_vector const& background_reads, double p,
              double log_fold = 0.5, uint8_t max_k = 10,
              double memory_limit = 4, double p_ext = 0.01,
              uint8_t lookup_k = 10, bool prune = false)
      : signal_reads_(signal_reads),
        background_reads_(background_reads),
        seeds_(),
        sig_size_(),
        bg_size_(),
        p_(p),
        p_ext_(p_ext),
        fold_lim_(std::pow(2, log_fold)),
        memory_limit_(memory_limit),
        k_lim_(max_k),
        lookup_k_(lookup_k),
        prune_(prune) {
    // For calculating the sum of the read lengths.
    auto const read_length_sum{[](packed_read_vector const& reads) {
      return std::accumulate(reads.begin(), reads.end(), std::uint64_t{},
                             [](auto const acc, packed_read const& rr) {
                               return acc + rr.length;
                             });
    }};

    gsl_set_error_handler_off();
    sig_size_ = read_length_sum(signal_reads_);
    bg_size_ = read_length_sum(background_reads_);
    signal_to_total_length_ratio_ = double(sig_size_) / (sig_size_ + bg_size_);
    std::cout << "# Background length " << bg_size_ << '\n';
    std::cout << "# Signal length " << sig_size_ << '\n';
    std::cout << "# Signal to total length ratio "
              << signal_to_total_length_ratio_ << '\n';
  }

  void find_seeds();

  void set_discarded_gapmer_reporting_ostream(std::ostream& stream) {
    discarded_gapmer_reporting_ostream_ = &stream;
  }

  const std::vector<seed>& get_seeds() const { return seeds_; }

  double signal_to_total_length_ratio() const {
    return signal_to_total_length_ratio_;
  }
};


// Thread safe.
template <typename t_configuration>
template <typename t_lhs_info, typename t_rhs_info>
inline void seed_finder<t_configuration>::report_discarded(
    gapmer_type lhs, gapmer_type rhs, extension_validation_result res,
    t_lhs_info const& lhs_info, t_rhs_info const& rhs_info,
    char const* test_fn) const {
  if constexpr (enable_reporting_discarded_seeds) {
    libbio_assert(discarded_gapmer_reporting_ostream_);
    // GCC 14 reports an error if std::osyncstream is used directly.
    libbio::osyncstream sync_stream{*discarded_gapmer_reporting_ostream_};
    std::ostream &stream{sync_stream};

    stream << "extension\t";

    if (res)
      stream << lhs << '\t' << rhs;
    else
      stream << rhs << '\t' << lhs;

    stream << '\t' << bool(res) << "\t\t" << test_fn << '\t';

    stream << res.strategy << "; " << res.status << "; ";
    if (res)
      stream << " (lhs: " << lhs_info << " rhs: " << rhs_info << ')';
    else
      stream << " (lhs: " << rhs_info << " rhs: " << lhs_info << ')';
  }
}


// Thread safe.
template <typename t_configuration>
template <typename t_value>
inline void seed_finder<t_configuration>::report_encrichment_check_failure(
    gapmer_type discarded, enrichment_check_result<t_value> res,
    char const* caller, char const* test_fn) const {
  if constexpr (enable_reporting_discarded_seeds) {
    libbio_assert(discarded_gapmer_reporting_ostream_);
    // GCC 14 reports an error if std::osyncstream is used directly.
    libbio::osyncstream sync_stream{*discarded_gapmer_reporting_ostream_};
    std::ostream &stream{sync_stream};

    stream << "encrichment check\t\t" << discarded << '\t' << false << '\t'
           << caller << '\t' << test_fn << '\t' << res << '\n';
  }
}


/**
 * Check if extension aa -> bb is valid
 *
 * @param aa     Mer to extend from
 * @param bb     Mer to extend to
 * @param aa_er  Enrichment result for aa
 * @param bb_er  Enrichment result for bb
 */
template <typename t_configuration>
template <typename t_value>
[[nodiscard]] auto seed_finder<t_configuration>::validate_extension(
    gapmer_type aa, gapmer_type bb, enrichment_result<t_value> aa_er,
    enrichment_result<t_value> bb_er) const -> extension_validation_result {
  if (aa_er.background_count <= 1.00001 && bb_er.background_count <= 1.00001) {
    constexpr auto strategy{extension_validation_strategy::signal_count_only};
    if (aa_er.signal_count > bb_er.signal_count * 4) {
      return extension_validation_result::make(
          strategy,
          extension_validation_status::source_signal_count_at_least_quadruple);
    }

    // Check the AC test result.
    // FIXME: Is this needed since it uses the background count?
    if (p_ <= bb_er.ac_test_result) {
      return extension_validation_result::make(
          strategy, extension_validation_status::ac_test_failed);
    }

    // Binomial test.
    double const p_ext{
        gsl_cdf_binomial_Q(bb_er.signal_count, 0.25, aa_er.signal_count)};
    if (p_ext_ <= p_ext) {
      return extension_validation_result::make(
          strategy, extension_validation_status::binomial_test_failed, p_ext);
    }

    // Success.
    return extension_validation_result::make(
        strategy, extension_validation_status::binomial_and_ac_tests_passed);
  }

  // Either a or b occurs in the background set.
  if constexpr (filter_mers) {
    constexpr auto strategy{
        extension_validation_strategy::comparison_to_background_filtering};
    if (bb_er.ac_test_result <= aa_er.ac_test_result) {
      return extension_validation_result::make(
          strategy,
          extension_validation_status::target_ac_test_p_value_smaller);
    }

    return extension_validation_result::make(
        strategy, extension_validation_status::source_ac_test_p_value_smaller);
  } else {
    constexpr auto strategy{
        extension_validation_strategy::comparison_to_background};
    if (bb_er.ac_test_result <= p_) {
      return extension_validation_result::make(
          strategy, extension_validation_status::ac_test_passed);
    }

    return extension_validation_result::make(
        strategy, extension_validation_status::ac_test_failed);
  }
}


/**
 * Check if aa should be discarded based on the neigbour bb.
 *
 * @tparam t_calculate_ac_test_for_bb  Indicates whether bb_er.ac_test_result
 * needs to be computed by this function or has been precalculated.
 *
 * @param aa         Mer to check
 * @param bb         Neighbour of aa
 * @param aa_er      Enrichment result for aa
 * @param bb_er      Enrichment result for bb
 */
template <typename t_configuration>
template <bool t_calculate_ac_test_for_bb, typename t_value>
[[nodiscard]] bool seed_finder<t_configuration>::should_filter(
    gapmer_type aa, gapmer_type bb, enrichment_result<t_value> aa_er,
    enrichment_result<t_value>& bb_er) const {
  if (aa_er.background_count <= 1.00001 && bb_er.background_count <= 1.00001) {
    if (aa_er.signal_count < bb_er.signal_count) {
      if constexpr (t_calculate_ac_test_for_bb) {
        bb_er.ac_test_result = error_suppressed_beta_inc(
            bb_er.signal_count, bb_er.background_count,
            signal_to_total_length_ratio_);
      }
      return true;
    }
    return false;
  }

  if constexpr (t_calculate_ac_test_for_bb) {
    bb_er.ac_test_result =
        error_suppressed_beta_inc(bb_er.signal_count, bb_er.background_count,
                                  signal_to_total_length_ratio_);
  }
  return bb_er.ac_test_result < aa_er.ac_test_result;
}


template <typename t_configuration>
template <typename t_value>
[[nodiscard]] auto seed_finder<t_configuration>::check_enrichment(
    count_pair<t_value> count) const -> enrichment_check_result<t_value> {
  // Add pseudocounts.
  count.add_pseudocounts();

  // Check the enrichment w.r.t. background by requiring that
  // the logarithmic fold change is greater than log_fold passed
  // to the constructor, 0.5 by default (with the added pseudocount, see the
  // definitions of aa and bb). The formula we use is
  // log_2((sig_count / sig_size) / (bg_count / bg_size)) ≤ log_fold
  // but with the logarithm and divisions removed.
  auto const& [sc, bc] = count;
  if (sc * bg_size_ <= fold_lim_ * bc * sig_size_) {
    return {{0, sc, bc}, enrichment_check_status::fold_change_test_failed};
  }

  // Check the enrichment w.r.t. background by using the Audic-Claverie test.
  // (Using the incomplete regularized beta function is explained in
  // the supplement of Jean-Michel Claverie, Thi Ngan Ta, ACDtool: a
  // web-server for the generic analysis of large data sets of counts,
  // Bioinformatics, Volume 35, Issue 1, January 2019, Pages 170–171,
  // https://doi.org/10.1093/bioinformatics/bty640)
  double const rr{
      error_suppressed_beta_inc(sc, bc, signal_to_total_length_ratio_)};
  return {{rr, sc, bc},
          rr <= p_ ? enrichment_check_status::success
                   : enrichment_check_status::ac_test_failed};
}


/**
 * Check the enrichment of gg in the given signal and background.
 *
 * @param gg The gapmer in question.
 * @param counts The signal and background counts.
 * @param critical OpenMP’s critical section wrapper.
 */
template <typename t_configuration>
template <typename t_critical>
[[nodiscard]] auto seed_finder<t_configuration>::check_enrichment_and_filter(
    gapmer_type const gg, uint64_t const offset, gapmer_count_type& counts,
    t_critical&& critical) const
    -> enrichment_check_result<gapmer_count_value_type> {
  auto const vv{gg.value()};
  if constexpr (filter_mers) {
    if (critical([&] { return counts.is_discarded_(vv, offset); })) {
      return {{0, 0, 0}, enrichment_check_status::discarded_earlier};
    }
  }

  if (not gg.is_canonical()) {
    if constexpr (filter_mers) {
      critical([&] { counts.mark_discarded_(vv, offset); });
    }
    return {{0, 0, 0}, enrichment_check_status::not_canonical};
  }

  auto const retval{check_enrichment(counts.count(gg, offset))};
  if (not retval) {
    if constexpr (filter_mers) {
      critical([&] { counts.mark_discarded_(vv, offset); });
    }
  }

  return retval;
}


/// Filter the Huddinge-1 neighbourhood of the given gapmer based on the
/// enrichment result.
/// \param gg             The gapmer in question.
/// \param offset         Precalculated offset in counts or prev_counts.
/// \param enrichment_res The enrichment result from check_enrichment().
/// \param counts         The counts for k.
/// \param prev_counts    The counts for k - 1 iff. t_should_consider_longer_kmers is true.
/// \param critical       Functor for accessing counts within a critical section.
/// \param prev_critical  Functor for accessing prev_counts within a critical section.
///
template <typename t_configuration>
template <bool t_should_extend, typename t_critical, typename t_prev_critical>
void seed_finder<t_configuration>::filter_huddinge_neighbourhood(
    gapmer_type const gg, uint64_t const offset,
    enrichment_result<gapmer_count_value_type> const& enrichment_res,
    gapmer_count_type& counts, t_critical&& critical,
    gapmer_count_type* const prev_counts,
    t_prev_critical&& prev_critical) const {
  constexpr bool skip_same_length{t_should_extend};
  constexpr bool skip_longer{not t_should_extend};

  auto const& [rr, sc, bc] = enrichment_res;
  gg.template huddinge_neighbours<true, skip_same_length, skip_longer>(
      [&](gapmer_type oo) {
        auto const o_offset{counts.offset(oo.gap_start(), oo.gap_length())};
        auto const o_counts{counts.count(oo, o_offset)};
        double const osc{o_counts.signal_count + 1};
        double const obc{o_counts.background_count + 1};
        if (osc * bg_size_ <= fold_lim_ * obc * sig_size_) {
          // FIXME: report discarded.
          critical([&] { counts.mark_discarded(oo, o_offset); });
          return;
        }

        if constexpr (t_should_extend) {
          double const o_r{error_suppressed_beta_inc(
              osc, obc, signal_to_total_length_ratio_)};

          enrichment_result other_enrichment_res{o_r, osc, obc};
          auto const res{
              validate_extension(gg, oo, enrichment_res, other_enrichment_res)};
          report_discarded(gg, oo, res, enrichment_res, other_enrichment_res,
                           "filter_huddinge_neighbourhood");

          if (res) {
            prev_critical([&] { prev_counts->mark_discarded(gg, offset); });
          } else {
            critical([&] { counts.mark_discarded(oo, o_offset); });
          }
        } else {
          // FIXME: report discarded.
          double o_r{};
          if (should_filter<true>(gg, oo, sc, bc, osc, obc, rr, o_r)) {
            critical([&] { counts.mark_discarded(gg, offset); });
          } else {
            critical([&] { counts.mark_discarded(oo, o_offset); });
          }
        }
      });
}


/**
 * Checks H1 neighbourhood with lengths k and k + 1 of the given gapmer,
 * and marks invalid candidates as discarded in counts_k and counts_k1.
 * Surviving mers get added to candidate list.
 *
 * @param gg        Gapmer in question
 * @param offset    Offset value for table access
 * @param counts_k  Length k gapmer count tables
 * @param counts_k1 Length k + 1 gapmer count tables
 */
template <typename t_configuration>
auto seed_finder<t_configuration>::check_count(gapmer_type const gg,
                                               uint64_t offset,
                                               gapmer_count_type& counts_k,
                                               gapmer_count_type& counts_k1)
    -> enrichment_check_result<gapmer_count_value_type> {
  auto res{check_enrichment_and_filter(gg, offset, counts_k, critical_a_bv{})};
  if (not res) return res;

  auto const& enrichment_res{res.result};
  auto const& [rr, sc, bc] = enrichment_res;
  if constexpr (filter_mers) {
    // FIXME: Add discarded seed reporting to this branch?
    filter_huddinge_neighbourhood<false>(gg, offset, enrichment_res, counts_k,
                                         critical_a_bv{});
    filter_huddinge_neighbourhood<true>(gg, offset, enrichment_res, counts_k1,
                                        critical_o_bv{}, &counts_k,
                                        critical_a_bv{});

    // FIXME: Not necessarily atomic. (Likely works on x86-64.)
    if (not counts_k.is_discarded(gg, offset)) {
#pragma omp critical
      seeds_.emplace_back(gg, rr, uint64_t(sc), uint64_t(bc));
    } else {
      res.status = enrichment_check_status::discarded_earlier;
    }
  } else {
#pragma omp critical
    seeds_.emplace_back(gg, rr, uint64_t(sc), uint64_t(bc));
  }

  return res;
}


/**
 * Checks H1 neighbourhood of the given gapmer
 * and marks invalid candidates as discarded in counts_k.
 * Surviving mers get added to candidate set mm.
 *
 * @param gg       Gapmer in question
 * @param offset   Offset value for table access
 * @param counts_k Length k gapmer count tables
 * @param mm       Map to add candidate seeds to.
 */
template <typename t_configuration>
auto seed_finder<t_configuration>::filter_count(gapmer_type const gg,
                                                uint64_t offset,
                                                gapmer_count_type& counts_k,
                                                enrichment_result_map& mm) const
    -> enrichment_check_result<gapmer_count_value_type> {
  auto res{check_enrichment_and_filter(gg, offset, counts_k, critical_d_bv{})};
  if (not res) return res;

  if constexpr (filter_mers) {
    filter_huddinge_neighbourhood<false>(gg, offset, res.result, counts_k,
                                         critical_d_bv{});
    // FIXME: Not necessarily atomic. (Likely works on x86-64.)
    if (not counts_k.is_discarded(gg, offset)) {
#pragma omp critical
      mm[gg] = res.result;
    } else {
      res.status = enrichment_check_status::discarded_earlier;
    }
  } else {
#pragma omp critical
    mm[gg] = res.result;
  }

  return res;
}


/**
 * Generate candidates for all k small enough to enable full k-mer counting
 *
 * @param counts  Output parameter, for storing counts for final k-length
 * mers
 */
template <typename t_configuration>
void seed_finder<t_configuration>::count_short_gapmers(
    gapmer_count_type& counts) {
  using std::swap;  // For ADL.

  std::cerr << "Lookup tables up to " << +lookup_k_ << '\n';
  // Initialize by counting 5-mers
  gapmer_count_type counts_k(signal_reads_, background_reads_, 5);
  if constexpr (enable_smoothing) {
    counts_k.smooth();
  }

  for (uint8_t k = 6; k <= lookup_k_; ++k) {
    // Initialize the k + 1 to compute extensions
    gapmer_count_type counts_k1(signal_reads_, background_reads_, k);
    if constexpr (enable_smoothing) {
      counts_k1.smooth();
    }
    uint64_t v_lim = gapmer_count_type::ONE << ((k - 1) * 2);
#pragma omp parallel for
    for (uint64_t v = 0; v < v_lim; ++v) {
      gapmer_type const gg(v, k - 1, 0, 0);
      if (auto const res{check_count(gg, 0, counts_k, counts_k1)}; not res) {
        report_encrichment_check_failure(gg, res, "count_short_gapmers",
                                         "check_count");
      }
    }
    uint8_t gap_s = middle_gap_only ? (k - 1) / 2 : 1;
    uint8_t gap_lim = middle_gap_only ? k - gap_s - 1 : k - 2;
    for (; gap_s <= gap_lim; ++gap_s) {
      for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
        uint64_t offset = counts_k.offset(gap_s, gap_l);
#pragma omp parallel for
        for (uint64_t v = 0; v < v_lim; ++v) {
          gapmer_type gg(v, k - 1, gap_s, gap_l);
          if (auto const res{check_count(gg, offset, counts_k, counts_k1)};
              not res) {
            report_encrichment_check_failure(gg, res, "count_short_gapmers",
                                             "check_count");
          }
        }
      }
    }

    // k -> k + 1
    swap(counts_k, counts_k1);
    std::cerr << +(k - 1) << " → " << seeds_.size() << " candidates.\n";
  }
  swap(counts_k, counts);
}


/**
 * Check collected mers for extension validity
 *
 * Newly found valid extensions from aa to elements found in counts are
 * added to bb
 *
 * @param aa     Shorter mers to extend
 * @param bb     Valid extended mers
 * @param counts Partial k-mer counts that may contain valid mer
 * extensions
 */
template <typename t_configuration>
void seed_finder<t_configuration>::extend_counted(
    enrichment_result_map const& aa, enrichment_result_map& bb,
    partial_count_type const& counts) const {
  for (auto kv : aa) {
    kv.first.template huddinge_neighbours<true, true, false>(
        [&](gapmer_type oo) {
          if (not oo.is_canonical()) {
            oo = oo.reverse_complement();
          }

          if (not bb.contains(oo)) {
            auto const count([&] {
              if constexpr (enable_smoothing)
                return counts.smooth_count(oo);
              else
                return counts.count(oo);
            }());

            auto const enrichment_res{check_enrichment(count)};
            if (not enrichment_res) {
              report_encrichment_check_failure(oo, enrichment_res, "extend",
                                               "extend_counted");
              return;
            }

            auto const enrichment_res_fp{
                enrichment_res.result.template to_enrichment_result<double>()};
            auto const res{
                validate_extension(kv.first, oo, kv.second, enrichment_res_fp)};
            report_discarded(kv.first, oo, res, kv.second, enrichment_res,
                             "extend_counted");
            if (res) {
              bb[oo] = enrichment_res_fp;
            }
          }
        });
  }
}


/**
 * Find k length extensions from a, and store valid extensions in b.
 *
 * @param aa     k - 1 length mers to extend.
 * @param bb     storage for valid k-mers.
 * @param counts Structure to use with partial counting.
 * @param k      length of k-mers to find.
 * @param prune  Should only one pass of extensions be done.
 */
template <typename t_configuration>
void seed_finder<t_configuration>::extend(enrichment_result_map& aa,
                                          enrichment_result_map& bb,
                                          partial_count_type& counts,
                                          uint16_t k, bool prune) const {
  std::cerr << "    Extend " << aa.size() << " mers.\n";

  const constexpr double fill_limit = 0.4;
  gapmer_set del_set;

  auto const init_counters{[&](gapmer_type oo) {
    if (not bb.contains(oo)) {
      counts.init(oo);
      if constexpr (enable_smoothing) {
        oo.hamming_neighbours([&](gapmer_type h_n) { counts.init(h_n); });
      }
    }
  }};

  if (prune) {
    std::vector<seed> prio;
    for (auto kv : aa) {
      prio.emplace_back(seed::from_enrichment_result(kv.first, kv.second));
    }
    // Sort by fold change.
    std::sort(prio.begin(), prio.end(), [](const auto& lhs, const auto& rhs) {
      return (lhs.sig_count / lhs.bg_count) > (rhs.sig_count / rhs.bg_count);
    });

    for (auto res : prio) {
      res.g.template huddinge_neighbours<true, true, false>(init_counters);
      if (counts.fill_rate() >= fill_limit) {
        break;
      }
    }
  } else {
    for (auto kv : aa) {
      kv.first.template huddinge_neighbours<true, true, false>(init_counters);
      if (counts.fill_rate() >= fill_limit) {
        std::print("\tLoad factor ≥ {} ({}). Counting…\n", fill_limit,
                   counts.fill_rate());
        counts.count_mers(signal_reads_, background_reads_, k);
        std::cerr << "\tFiltering extension…\n";
        extend_counted(aa, bb, counts);
        counts.clear();
      }
    }
  }

  std::cerr << "\tFinal load factor " << counts.fill_rate() << ". Counting…\n";
  counts.count_mers(signal_reads_, background_reads_, k);
  std::cerr << "\tFiltering extension…\n";
  extend_counted(aa, bb, counts);
  counts.clear();

  if constexpr (filter_mers) {
    std::cerr << "\tFiltering sources…\n";
    for (auto kv : aa) {
      bool keep = true;

      kv.first.template huddinge_neighbours<true, true, false>(
          [&](gapmer_type o) {
            if (not o.is_canonical()) {
              o = o.reverse_complement();
            }
            if (bb.contains(o)) {
              auto const res{validate_extension(kv.first, o, kv.second, bb[o])};
              report_discarded(kv.first, o, res, kv.second, bb[o], "extend");
              if (res) {
                keep = false;
              }
            }
          });

      if (not keep) {
        del_set.insert(kv.first);
      }
    }

    for (auto d : del_set) {
      aa.erase(d);
    }
  }  // if constexpr (filter_mers)
}


/**
 * Filter found mers of the same lengths. Keep only the best within H1
 * distance among same length gapmers.
 *
 * @param m   gapmers to filter
 */
template <typename t_configuration>
void seed_finder<t_configuration>::filter(enrichment_result_map& mm) const {
  gapmer_set del_set;
  auto e = mm.end();
  for (auto it = mm.begin(); it != e; ++it) {
    bool keep = true;
    if (it->second.p > p_) {
      del_set.insert(it->first);
      continue;
    }
    auto iit = it;
    ++iit;
    for (; iit != e; ++iit) {
      if (it->first.template is_neighbour<true>(iit->first)) {
        if (it->second.p > iit->second.p) {
          del_set.insert(iit->first);
        } else if (it->second.p < iit->second.p && keep) {
          del_set.insert(it->first);
          keep = false;
        }
      }
    }
  }
  for (auto d : del_set) {
    mm.erase(d);
  }
  std::cerr << "    filtered to " << mm.size() << " mers\n";
}


/**
 * Does the heavy lifting of counting increasingly long gapmers to find seed
 * candidates
 */
template <typename t_configuration>
void seed_finder<t_configuration>::find_seeds() {
  using std::swap;  // For ADL.

  enrichment_result_map aa;
  enrichment_result_map bb;

  // full k-mer couting as long as memory is sufficient.
  {
    gapmer_count_type counts;
    count_short_gapmers(counts);
    uint64_t const v_lim{gapmer_count_type::ONE << (lookup_k_ * 2)};
#pragma omp parallel for
    for (uint64_t v = 0; v < v_lim; ++v) {
      gapmer_type const gg(v, lookup_k_, 0, 0);
      if (auto const res{filter_count(gg, 0, counts, aa)}; not res) {
        report_encrichment_check_failure(gg, res, "find_seeds", "filter_count");
      }
    }
    uint8_t gap_s = middle_gap_only ? lookup_k_ / 2 : 1;
    uint8_t const gap_lim(middle_gap_only ? lookup_k_ - gap_s : lookup_k_ - 1);
    for (; gap_s <= gap_lim; ++gap_s) {
      for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
        uint64_t offset = counts.offset(gap_s, gap_l);
#pragma omp parallel for
        for (uint64_t v = 0; v < v_lim; ++v) {
          gapmer_type const gg(v, lookup_k_, 0, 0);
          if (auto const res{filter_count(gg, offset, counts, aa)}; not res) {
            report_encrichment_check_failure(gg, res, "find_seeds",
                                             "filter_count");
          }
        }
      }
    }
  }

  // Partial count with extensions when we can no longer count everything
  partial_count_type partial_counts;
  for (uint8_t k = lookup_k_ + 1; k <= k_lim_; ++k) {
    std::cerr << +(k - 1) << " → \n";
    extend(aa, bb, partial_counts, k, prune_);
    std::print("    {} {} candidates\n", aa.size(), +(k - 1));
    std::print("    {} {} potentials\n", bb.size(), +k);
    for (auto pp : aa) {
      seeds_.emplace_back(seed::from_enrichment_result(pp.first, pp.second));
    }
    aa.clear();
    if constexpr (filter_mers) {
      filter(bb);
    }
    swap(aa, bb);
    std::print("{} → {} candidates.\n", +(k - 1), seeds_.size());
    if (aa.empty()) {
      break;
    }
  }
  if (aa.size() > 0) {
    for (auto pp : aa) {
      seeds_.emplace_back(seed::from_enrichment_result(pp.first, pp.second));
    }
    std::print("{} → {} candidates.\n", +k_lim_, seeds_.size());
  }
}
}  // namespace sf
