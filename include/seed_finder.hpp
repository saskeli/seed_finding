#pragma once

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "fm_index.hpp"
#include "gapmer.hpp"
#include "gapmer_count.hpp"
#include "packed_read.hpp"
#include "partial_count.hpp"
#include "util.hpp"

namespace sf {

template <bool middle_gap_only, uint8_t max_gap, bool enable_smoothing = true,
          bool filter_mers = true>
class seed_finder {
 public:
  typedef gapmer<middle_gap_only, max_gap> gapmer_type;

 private:
  typedef gapmer_count<gapmer_type> gapmer_count_type;
  typedef partial_count<gapmer_type> partial_count_type;
  typedef typename gapmer_count_type::value_type gapmer_count_value_type;

  typedef std::unordered_set<gapmer_type, typename gapmer_type::hash>
      gapmer_set;

  template <typename t_value>
  using gapmer_map =
      std::unordered_map<gapmer_type, t_value, typename gapmer_type::hash>;

  struct Res {
    gapmer_type g;
    double p;
    uint64_t sig_count;
    uint64_t bg_count;
  };

  typedef gapmer_map <Res> gapmer_res_map;

  packed_read_vector const& signal_reads_;
  packed_read_vector const& background_reads_;
  std::vector<Res> seeds_;
  uint64_t sig_size_;
  uint64_t bg_size_;
  double p_;
  double p_ext_;
  double fold_lim_;
  double x_;
  double memory_limit_;
  uint8_t k_lim_;
  uint8_t lookup_k_;
  bool prune_;

  /**
   * Check if extension a -> b is valid
   *
   * @param a      Mer to extend from
   * @param b      Mer to extend to
   * @param a_sig  Number of occurrences of a in the signal set
   * @param a_bg   Number of occurrences of a in the backround set
   * @param b_sig  Number of occurrences of b in the signal set
   * @param b_bg   Number of occurrences of b in the background set
   * @param a_r    p-value for a from the incomplete beta function
   * @param b_r    p-value for b from the incomplete beta function
   *
   * @return True, if extension from a to b is valid
   */
  template <bool debug = false>
  bool validate_extension([[maybe_unused]] gapmer_type a,
                          [[maybe_unused]] gapmer_type b, double a_sig,
                          double a_bg, double b_sig, double b_bg, double a_r,
                          double b_r) const {
    if (a_bg <= 1.00001 && b_bg <= 1.00001) {
      if (a_sig > b_sig * 4) {
        return false;
      }
      double p_extend = gsl_cdf_binomial_Q(b_sig, 0.25, a_sig);
      if (p_extend < p_ext_ && b_r < p_) {
#ifdef DEBUG
        std::cerr << "        " << a.to_string() << " discarded by "
                  << b.to_string() << "\n            (" << a_sig << ", " << a_bg
                  << ") <-> (" << b_sig << ", " << b_bg
                  << ")\n            with p_extend " << p_extend << std::endl;
#endif
        return true;
      }
#ifdef DEBUG
      std::cerr << "        " << a.to_string() << " discards " << b.to_string()
                << "\n            (" << a_sig << ", " << a_bg << ") <-> ("
                << b_sig << ", " << b_bg << ")\n            with p_extend "
                << p_extend << std::endl;
#endif
      return false;
    }
    if constexpr (filter_mers) {
#ifdef DEBUG
      if (b_r < a_r) {
        std::cerr << "        " << a.to_string() << " discarded by "
                  << b.to_string() << "\n            (" << a_sig << ", " << a_bg
                  << ") <-> (" << b_sig << ", " << b_bg
                  << ")\n            with p " << b_r << std::endl;
      } else {
        std::cerr << "        " << a.to_string() << " discards "
                  << b.to_string() << "\n            (" << a_sig << ", " << a_bg
                  << ") <-> (" << b_sig << ", " << b_bg
                  << ")\n            with p " << b_r << std::endl;
      }
#endif
      return b_r <= a_r;
    } else {
      return b_r <= p_;
    }
  }

  /**
   * Check if a should be discarded based on the neigbour b.
   *
   * @tparam calc_b_r  Indicates, wether b_r needs to be computed by this
   * function, or has been precalculated.
   *
   * @param a          Mer to check
   * @param b          Neighbour of a
   * @param a_sig      Number of occurrences of a in signal set
   * @param a_bg       Number of occurrences of a in backgroud set
   * @param b_sig      Number of occurrences of a in signal set
   * @param b_bg       Number of occurrences of b in backgroud set
   * @param a_r        p-value for a from incomplete beta function
   * @param b_r        p-value for b from incomplete beta function, if calc_b_r
   * == false, else output parameter.
   */
  template <bool calc_b_r>
  bool should_filter([[maybe_unused]] gapmer_type a, [[maybe_unused]] gapmer_type b,
                 double a_sig, double a_bg, double b_sig, double b_bg,
                 double a_r, double& b_r) const {
    if (a_bg <= 1.00001 && b_bg <= 1.00001) {
      if (b_sig > a_sig) {
        if constexpr (calc_b_r) {
          b_r = error_suppressed_beta_inc(b_sig, b_bg, x_);
        }
        return true;
      }
      return false;
    }
    if constexpr (calc_b_r) {
      b_r = error_suppressed_beta_inc(b_sig, b_bg, x_);
    }
    return b_r < a_r;
  }


  struct enrichment_result {
    double ac_test_result{};
    uint64_t signal_count{};
    uint64_t background_count{};
  };

  struct check_enrichment_result {
    enrichment_result result{};
    bool did_pass{};
  };

  /**
   * Check the enrichment of gg in the given signal and background.
   *
   * @param gg The gapmer in question.
   * @param counts The signal and background counts.
   * @param critical OpenMP’s critical section wrapper.
   */
  template <typename t_critical>
  [[nodiscard]] check_enrichment_result check_enrichment(
      gapmer_type const gg, uint64_t const offset, gapmer_count_type &counts,
      t_critical&& critical) const {
    auto const vv{gg.value()};

    if constexpr (filter_mers) {
      if (counts.is_discarded_(vv, offset)) {
        return {{0, 0, 0}, false};
      }
    }

    if (not gg.is_canonical()) {
      if constexpr (filter_mers) {
        critical([&] { counts.mark_discarded_(vv, offset); });
      }
      return {{0, 0, 0}, false};
    }

    // Add pseudocounts.
    // Narrows the count value type (double).
    auto const [sc_, bc_] = counts.count(gg, offset);
    uint64_t const sc(sc_ + 1);
    uint64_t const bc(bc_ + 1);

    // Check the enrichment w.r.t. background by requiring that
    // the logarithmic fold change is greater than log_fold passed
    // to the constructor, 0.5 by default (with the added pseudocount, see the
    // definitions of aa and bb). The formula we use is
    // log_2((sig_count / sig_size) / (bg_count / bg_size)) ≤ log_fold
    // but with the logarithm and divisions removed.
    if (sc * bg_size_ <= fold_lim_ * bc * sig_size_) {
      if constexpr (filter_mers) {
        critical([&] { counts.mark_discarded_(vv, offset); });
      }
      return {{0, sc, bc}, false};
    }

    // Check the enrichment w.r.t. background by using the Audic-Claverie test.
    // (Using the incomplete regularized beta function is explained in
    // the supplement of Jean-Michel Claverie, Thi Ngan Ta, ACDtool: a
    // web-server for the generic analysis of large data sets of counts,
    // Bioinformatics, Volume 35, Issue 1, January 2019, Pages 170–171,
    // https://doi.org/10.1093/bioinformatics/bty640)
    double const rr{error_suppressed_beta_inc(sc, bc, x_)};
    if (rr > p_) {
      if constexpr (filter_mers) {
        critical([&] { counts.mark_discarded_(vv, offset); });
      }
      return {{rr, sc, bc}, false};
    }

    return {{rr, sc, bc}, true};
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
  template <bool t_should_extend, typename t_critical,
            typename t_prev_critical = std::nullptr_t>
  void filter_huddinge_neighbourhood(
      gapmer_type const gg, uint64_t const offset,
      enrichment_result const& enrichment_res, gapmer_count_type& counts,
      t_critical&& critical, gapmer_count_type* const prev_counts = nullptr,
      t_prev_critical&& prev_critical = nullptr) const {

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
            critical([&] { counts.mark_discarded(oo, o_offset); });
            return;
          }

          if constexpr (t_should_extend) {
            double const o_r{error_suppressed_beta_inc(osc, obc, x_)};
            if (validate_extension(gg, oo, sc, bc, osc, obc, rr, o_r)) {
              prev_critical([&] { prev_counts->mark_discarded(gg, offset); });
            } else {
              critical([&] { counts.mark_discarded(oo, o_offset); });
            }
          } else {
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
   * Checks H1 neighbourhood with lengths k and k + 1 of gapmer implied by k, v,
   * gap_s and gap_l, and marks invalid candidates as discarded in sig_bg_a an
   * sig_bg_b.
   *
   * Surviving mers get added to candidate list.
   *
   * @param k         Gapmer length
   * @param v         Gapmer value as uint64_t
   * @param gap_s     Start location of gap
   * @param gap_l     Length of gap
   * @param offset    Offset value for table access
   * @param sig_bg_k  Length k gapmer count tables
   * @param sig_bg_k1 Length k + 1 gapmer count tables
   */
  void check_count(gapmer_type const gg, uint64_t offset, gapmer_count_type& sig_bg_k,
                   gapmer_count_type& sig_bg_k1) {
    auto const& [enrichment_res, should_continue] =
        check_enrichment(gg, offset, sig_bg_k, critical_a_bv{});
    if (not should_continue) return;

    auto const &[rr, sc, bc] = enrichment_res;
    if constexpr (filter_mers) {
      filter_huddinge_neighbourhood<false>(gg, offset, enrichment_res, sig_bg_k,
                                           critical_a_bv{});
      filter_huddinge_neighbourhood<true>(gg, offset, enrichment_res, sig_bg_k1,
                                          critical_o_bv{}, &sig_bg_k,
                                          critical_a_bv{});

      if (not sig_bg_k.is_discarded(gg, offset)) {
#pragma omp critical
        seeds_.push_back({gg, rr, sc, bc});
      }
    } else {
#pragma omp critical
      seeds_.push_back({gg, rr, sc, bc});
    }
  }

  /**
   * Checks H1 neighbourhood with lengths k gapmer implied by k, v,
   * gap_s and gap_l, and marks invalid candidates as discarded in sig_bg_c.
   *
   * Surviving mers get added to candidate set m.
   *
   * @param k        Gapmer length
   * @param v        Gapmer value as uint64_t
   * @param gap_s    Start location of gap
   * @param gap_l    Length of gap
   * @param offset   Offset value for table access
   * @param sig_bg_k Length k gapmer count tables
   * @param mm       Map to add candidate seeds to.
   */
  void filter_count(gapmer_type const gg, uint64_t offset, gapmer_count_type &sig_bg_k,
                    gapmer_res_map &mm) const {
    auto const& [enrichment_res, should_continue] =
        check_enrichment(gg, offset, sig_bg_k, critical_d_bv{});
    if (not should_continue) return;

    auto const &[rr, sc, bc] = enrichment_res;
    if constexpr (filter_mers) {
      filter_huddinge_neighbourhood<false>(gg, offset, enrichment_res, sig_bg_k, critical_d_bv{});
      if (not sig_bg_k.is_discarded(gg, offset)) {
#pragma omp critical
        mm[gg] = {gg, rr, sc, bc};
      }
    } else {
#pragma omp critical
      mm[gg] = {gg, rr, sc, bc};
    }
  }

  /**
   * Generate candidates for all k small enough to enable full k-mer counting
   *
   * @param sig_bg_c  Ouput parameter, for storing counts for final k-length
   * mers
   */
  void counted_seeds_and_candidates(gapmer_count_type& sig_bg_c) {
    using std::swap;

    std::cerr << "Lookup tables up to " << int(lookup_k_) << std::endl;
    // Initialize by counting 5-mers
    gapmer_count_type sig_bg_a(signal_reads_, background_reads_, 5);
    if constexpr (enable_smoothing) {
      sig_bg_a.smooth();
    }

    for (uint8_t k = 6; k <= lookup_k_; ++k) {
      // Initialize the k + 1 to compute extensions
      gapmer_count_type sig_bg_b(signal_reads_, background_reads_, k);
      if constexpr (enable_smoothing) {
        sig_bg_b.smooth();
      }
      uint64_t v_lim = gapmer_count_type::ONE << ((k - 1) * 2);
#pragma omp parallel for
      for (uint64_t v = 0; v < v_lim; ++v) {
        gapmer_type const gg(v, k - 1, 0, 0);
        check_count(gg, 0, sig_bg_a, sig_bg_b);
      }
      uint8_t gap_s = middle_gap_only ? (k - 1) / 2 : 1;
      uint8_t gap_lim = middle_gap_only ? k - gap_s - 1 : k - 2;
      for (; gap_s <= gap_lim; ++gap_s) {
        for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
          uint64_t offset = sig_bg_a.offset(gap_s, gap_l);
#pragma omp parallel for
          for (uint64_t v = 0; v < v_lim; ++v) {
            gapmer_type gg(v, k - 1, gap_s, gap_l);
            check_count(gg, offset, sig_bg_a, sig_bg_b);
          }
        }
      }

      // k -> k + 1
      swap(sig_bg_a, sig_bg_b);
      std::cerr << int(k - 1) << " -> " << seeds_.size() << " candidates."
                << std::endl;
    }
    swap(sig_bg_a, sig_bg_c);
  }

  /**
   * Check collected mers for extension validity
   *
   * Newly found valid extensions from a to elements found in p_counter are
   * added to b
   *
   * @param a  Shorter mers to extend
   * @param b  Valid extended mers
   * @param p_counter  Partial k-mer counts that may contain valid mer
   * extensions
   */
  void extend_counted(gapmer_res_map const &aa, gapmer_res_map &bb, partial_count_type const &p_counter) const {
    for (auto p : aa) {
#ifdef DEBUG
      std::cerr << "        " << p.first.to_string() << ": "
                << p.second.sig_count << ", " << p.second.bg_count << ", "
                << p.second.p << std::endl;
#endif
      p.first.template huddinge_neighbours<true, true, false>([&](gapmer_type oo) {
        if (not oo.is_canonical()) {
          oo = oo.reverse_complement();
        }

        if (not bb.contains(oo)) {
          auto const counts([&]{
            if constexpr (enable_smoothing)
              return p_counter.smooth_count(oo);
            else
              return p_counter.count(oo);
          }());

          // Potentially narrowing conversion.
          double const sc(counts.first + 1);
          double const bc(counts.second + 1);

          if (sc * bg_size_ <= fold_lim_ * bc * sig_size_) {
            return;
          }

          double const rr{error_suppressed_beta_inc(sc, bc, x_)};
          if (validate_extension<false>(p.first, oo, p.second.sig_count,
                                        p.second.bg_count, sc, bc,
                                        p.second.p, rr)) {
            bb[oo] = {oo, rr, uint64_t(sc), uint64_t(bc)};
          }
        }
      });
    }
  }

  /**
   * Find k length extensions from a, and store valid extensions in b.
   *
   * @param a  k - 1 length mers to extend.
   * @param b  storage for valid k-mers.
   * @param p_counter  Structure to use with partial counting.
   * @param k  length of k-mers to find.
   * @param prune  Should only one pass of extensions be done.
   */
   void extend(gapmer_res_map &aa, gapmer_res_map &bb,
               partial_count_type &p_counter, uint16_t k, bool prune) const {
    std::cerr << "    Extend " << aa.size() << " mers." << std::endl;

    const constexpr double fill_limit = 0.4;
    gapmer_set del_set;

    auto const init_counters{[&](gapmer_type oo) {
      if (not bb.contains(oo)) {
        p_counter.init(oo);
        if constexpr (enable_smoothing) {
          oo.hamming_neighbours([&](gapmer_type h_n) { p_counter.init(h_n); });
        }
      }
    }};

    if (prune) {
      std::vector<Res> prio;
      for (auto res : aa) {
        prio.push_back(res.second);
      }
      // Sort by fold change.
      std::sort(prio.begin(), prio.end(), [](const auto& lhs, const auto& rhs) {
        return (lhs.sig_count / lhs.bg_count) > (rhs.sig_count / rhs.bg_count);
      });

      for (auto res : prio) {
        res.g.template huddinge_neighbours<true, true, false>(init_counters);
        if (p_counter.fill_rate() >= fill_limit) {
          break;
        }
      }
    } else {
      for (auto res : aa) {
        res.first.template huddinge_neighbours<true, true, false>(init_counters);
        if (p_counter.fill_rate() >= fill_limit) {
          std::cerr << "\tLoad factor >= " << fill_limit << " ("
                    << p_counter.fill_rate() << ") counting.." << std::endl;
          p_counter.count_mers(signal_reads_, background_reads_, k);
          std::cerr << "\tFiltering extension..." << std::endl;
          extend_counted(aa, bb, p_counter);
          p_counter.clear();
        }
      }
    }

    std::cerr << "\tFinal load factor " << p_counter.fill_rate()
              << " counting.." << std::endl;
    p_counter.count_mers(signal_reads_, background_reads_, k);
    std::cerr << "\tFiltering extension..." << std::endl;
    extend_counted(aa, bb, p_counter);
    p_counter.clear();

    if constexpr (filter_mers) {
      std::cerr << "\tFiltering sources..." << std::endl;
      for (auto res : aa) {
#ifdef DEBUG
        std::cerr << "        " << res.first.to_string() << ": "
                  << res.second.sig_count << ", " << res.second.bg_count << ", "
                  << res.second.p << std::endl;
#endif
        bool keep = true;

        res.first.template huddinge_neighbours<true, true, false>(
            [&](gapmer_type o) {
              if (not o.is_canonical()) {
                o = o.reverse_complement();
              }
              if (bb.contains(o)) {
                if (validate_extension(res.first, o, res.second.sig_count,
                                       res.second.bg_count, bb[o].sig_count,
                                       bb[o].bg_count, res.second.p, bb[o].p)) {
                  keep = false;
                }
              }
            });

        if (not keep) {
          del_set.insert(res.first);
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
  void filter(gapmer_res_map &mm) const {
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
    std::cerr << "    filtered to " << mm.size() << " mers" << std::endl;
  }

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
    x_ = double(sig_size_) / (sig_size_ + bg_size_);
    std::cerr << "Background length " << bg_size_ << '\n';
    std::cerr << "Signal length " << sig_size_ << '\n';
    std::cerr << "X = " << x_ << '\n';
  }


  /**
   * Does the heavy lifting of counting increasingly long gapmers to find seed
   * candidates
   */
  void find_seeds() {
    using std::swap;

    gapmer_res_map aa;
    gapmer_res_map bb;

    // full k-mer couting as long as memory is sufficient.
    {
      gapmer_count_type sig_bg_c;
      counted_seeds_and_candidates(sig_bg_c);
      uint64_t const v_lim{gapmer_count_type::ONE << (lookup_k_ * 2)};
#pragma omp parallel for
      for (uint64_t v = 0; v < v_lim; ++v) {
        gapmer_type const gg(v, lookup_k_, 0, 0);
        filter_count(gg, 0, sig_bg_c, aa);
      }
      uint8_t gap_s = middle_gap_only ? lookup_k_ / 2 : 1;
      uint8_t const gap_lim(middle_gap_only ? lookup_k_ - gap_s : lookup_k_ - 1);
      for (; gap_s <= gap_lim; ++gap_s) {
        for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
          uint64_t offset = sig_bg_c.offset(gap_s, gap_l);
#pragma omp parallel for
          for (uint64_t v = 0; v < v_lim; ++v) {
            gapmer_type const gg(v, lookup_k_, 0, 0);
            filter_count(gg, offset, sig_bg_c, aa);
          }
        }
      }
    }

    // Partial count with extensions when we can no longer count everything
    partial_count_type p_counter;
    for (uint8_t k = lookup_k_ + 1; k <= k_lim_; ++k) {
      std::cerr << int(k) - 1 << " -> " << std::endl;
      extend(aa, bb, p_counter, k, prune_);
      std::cerr << "    " << aa.size() << " " << int(k) - 1 << " candidates\n"
                << "    " << bb.size() << " " << int(k) << " potentials"
                << std::endl;
      for (auto p : aa) {
        seeds_.push_back(p.second);
      }
      aa.clear();
      if constexpr (filter_mers) {
        filter(bb);
      }
      swap(aa, bb);
      std::cerr << int(k) - 1 << " -> " << seeds_.size() << " candidates."
                << std::endl;
      if (aa.size() == 0) {
        break;
      }
    }
    if (aa.size() > 0) {
      for (auto p : aa) {
        seeds_.push_back(p.second);
      }
      std::cerr << int(k_lim_) << " -> " << seeds_.size() << " candidates."
                << std::endl;
    }
  }

  const std::vector<Res>& get_seeds() const { return seeds_; }

  double x() const { return x_; }
};
}  // namespace sf
