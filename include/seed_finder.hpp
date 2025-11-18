#pragma once

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <libbio/utility.hh>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "fm_index.hpp"
#include "gapmer.hpp"
#include "gapmer_count.hpp"
#include "libbio_reader_adapter.hpp"
#include "partial_count.hpp"
#include "util.hpp"

namespace sf {

template <bool middle_gap_only, uint8_t max_gap, bool enable_smootihing = true,
          bool filter_mers = true>
class seed_finder : public libbio_reader_adapter_delegate {
 public:
  typedef gapmer<middle_gap_only, max_gap> gapmer_type;

 private:
  typedef gapmer_count<middle_gap_only, max_gap> gapmer_count_type;
  typedef std::set<std::string, libbio::compare_strings_transparent> path_set;

  struct Res {
    gapmer_type g;
    double p;
    uint64_t sig_count;
    uint64_t bg_count;
  };

  std::string sig_path_;
  std::string bg_path_;
  std::vector<Res> seeds_;
  path_set paths_with_errors_;
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
  bool do_extend([[maybe_unused]] gapmer_type a, [[maybe_unused]] gapmer_type b,
                 double a_sig, double a_bg, double b_sig, double b_bg,
                 double a_r, double b_r) {
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
  bool do_filter([[maybe_unused]] gapmer_type a, [[maybe_unused]] gapmer_type b,
                 double a_sig, double a_bg, double b_sig, double b_bg,
                 double a_r, double& b_r) {
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

  /**
   * Checks H1 neighbourhood with lengths k and k + 1 of gapmer implied by k, v,
   * gap_s and gap_l, and marks invalid candidates as discarded in sig_bg_a an
   * sig_bg_b.
   *
   * Surviving mers get added to candidate list.
   *
   * @param k        Gapmer length
   * @param v        Gapmer value as uint64_t
   * @param gap_s    Start location of gap
   * @param gap_l    Length of gap
   * @param offset   Offset value for table access
   * @param sig_bg_a Length k gapmer count tables
   * @param sig_bg_b Length k + 1 gapmer count tables
   */
  void check_count(const uint8_t k, const uint64_t v, uint8_t gap_s,
                   uint8_t gap_l, uint64_t offset, gapmer_count_type& sig_bg_a,
                   gapmer_count_type& sig_bg_b) {
#ifdef DEBUG
    if (offset + v >= gapmer_count_type::lookup_elems(k)) {
      std::cerr << "accessing " << offset << " + " << v << " = " << offset + v
                << " of " << gapmer_count_type::lookup_elems(k)
                << " element table" << std::endl;
      exit(1);
    }
#endif
    if constexpr (filter_mers) {
      if (sig_bg_a.discarded[offset + v]) {
        return;
      }
    }
    gapmer_type g(v, k, gap_s, gap_l);
    if (not g.is_canonical()) {
      if constexpr (filter_mers) {
#pragma omp critical(a_bv)
        sig_bg_a.discarded[offset + v] = true;
      }
      return;
    }
    uint64_t a = sig_bg_a.sig_counts[offset + v] + 1;
    uint64_t b = sig_bg_a.bg_counts[offset + v] + 1;
    if (a * bg_size_ <= fold_lim_ * b * sig_size_) {
      if constexpr (filter_mers) {
#pragma omp critical(a_bv)
        sig_bg_a.discarded[offset + v] = true;
      }
      return;
    }
    double r = error_suppressed_beta_inc(a, b, x_);
    if (r > p_) {
      if constexpr (filter_mers) {
#pragma omp critical(a_bv)
        sig_bg_a.discarded[offset + v] = true;
      }
      return;
    }
    if constexpr (filter_mers) {
      // Length k mers.
      auto callback_a = [&](gapmer_type o) {
        uint64_t o_offset = sig_bg_a.offset(o.gap_start(), o.gap_length());
        uint64_t o_v = o.value();
        double o_a = sig_bg_a.sig_counts[o_offset + o_v] + 1;
        double o_b = sig_bg_a.sig_counts[o_offset + o_v] + 1;
        if (o_a * sig_size_ <= fold_lim_ * o_b * bg_size_) {
#pragma omp critical(a_bv)
          sig_bg_a.discarded[o_offset + o_v] = true;
          return;
        }
        double o_r;
        if (do_filter<true>(g, o, a, b, o_a, o_b, r, o_r)) {
#pragma omp critical(a_bv)
          sig_bg_a.discarded[offset + v] = true;
        } else {
#pragma omp critical(a_bv)
          sig_bg_a.discarded[o_offset + o_v] = true;
        }
      };
      g.template huddinge_neighbours<true, false, true>(callback_a);

      // Length k + 1 mers.
      auto callback_b = [&](gapmer_type o) {
        uint64_t o_offset = sig_bg_b.offset(o.gap_start(), o.gap_length());
        uint64_t o_v = o.value();
        double o_a = sig_bg_b.sig_counts[o_offset + o_v] + 1;
        double o_b = sig_bg_b.bg_counts[o_offset + o_v] + 1;
        if (o_a * sig_size_ <= fold_lim_ * o_b * bg_size_) {
#pragma omp critical(o_bv)
          sig_bg_b.discarded[o_offset + o_v] = true;
          return;
        }
        double o_r = error_suppressed_beta_inc(o_a, o_b, x_);
        if (do_extend(g, o, a, b, o_a, o_b, r, o_r)) {
#pragma omp critical(a_bv)
          sig_bg_a.discarded[offset + v] = true;
        } else {
#pragma omp critical(o_bv)
          sig_bg_b.discarded[o_offset + o_v] = true;
        }
      };
      g.template huddinge_neighbours<true, true, false>(callback_b);
      if (sig_bg_a.discarded[offset + v] == false) {
#pragma omp critical
        seeds_.push_back({g, r, a, b});
      }
    } else {
#pragma omp critical
      seeds_.push_back({g, r, a, b});
    }
  }

  /**
   * Checks H1 neighbourhood with lengths k gapmer implied by k, v,
   * gap_s and gap_l, and marks invalid candidates as discarded in sig_bg_c.
   *
   * Surviving mers get added to candidate set m.
   *
   * @tparam M       Type of m.
   *
   * @param k        Gapmer length
   * @param v        Gapmer value as uint64_t
   * @param gap_s    Start location of gap
   * @param gap_l    Length of gap
   * @param offset   Offset value for table access
   * @param sig_bg_c Length k gapmer count tables
   * @param m        Map to add candidate seeds to.
   */
  template <class M>
  void filter_count(const uint8_t k, const uint64_t v, uint8_t gap_s,
                    uint8_t gap_l, uint64_t offset, gapmer_count_type& sig_bg_c,
                    M& m) {
#ifdef DEBUG
    if (offset + v >= gapmer_count_type::lookup_elems(k)) {
      std::cerr << "k = " << int(k) << " & sig_bg_c.k_ = " << int(sig_bg_c.k_)
                << " :\n accessing " << offset << " + " << v << " = "
                << offset + v << " of " << gapmer_count_type::lookup_elems(k)
                << " element table" << std::endl;
      exit(1);
    }
#endif
    if constexpr (filter_mers) {
      if (sig_bg_c.discarded[offset + v]) {
        return;
      }
    }
    gapmer_type g(v, k, gap_s, gap_l);
    if (not g.is_canonical()) {
      if constexpr (filter_mers) {
#pragma omp critical(d_bv)
        sig_bg_c.discarded[offset + v] = true;
      }
      return;
    }
    uint64_t a = sig_bg_c.sig_counts[offset + v] + 1;
    uint64_t b = sig_bg_c.bg_counts[offset + v] + 1;
    if (a * sig_size_ <= fold_lim_ * b * bg_size_) {
      if constexpr (filter_mers) {
#pragma omp critical(d_bv)
        sig_bg_c.discarded[offset + v] = true;
      }
      return;
    }
    double r = error_suppressed_beta_inc(a, b, x_);
    if (r > p_) {
      if constexpr (filter_mers) {
#pragma omp critical(d_bv)
        sig_bg_c.discarded[offset + v] = true;
      }
      return;
    }
    if constexpr (filter_mers) {
      auto callback = [&](gapmer_type o) {
        uint64_t o_offset = sig_bg_c.offset(o.gap_start(), o.gap_length());
        uint64_t o_v = o.value();
        double o_a = sig_bg_c.sig_counts[o_offset + v] + 1;
        double o_b = sig_bg_c.bg_counts[o_offset + v] + 1;
        if (o_a * sig_size_ <= fold_lim_ * o_b * bg_size_) {
#pragma omp critical(d_bv)
          sig_bg_c.discarded[o_offset + o_v] = true;
          return;
        }
        double o_r;
        if (do_filter<true>(g, o, a, b, o_a, o_b, r, o_r)) {
#pragma omp critical(d_bv)
          sig_bg_c.discarded[offset + v] = true;
        } else {
#pragma omp critical(d_bv)
          sig_bg_c.discarded[o_offset + o_v] = true;
        }
      };
      g.template huddinge_neighbours<true, false, true>(callback);
      if (sig_bg_c.discarded[offset + v] == false) {
#pragma omp critical
        m[g] = {g, r, a, b};
      }
    } else {
#pragma omp critical
      m[g] = {g, r, a, b};
    }
  }

  /**
   * Generate candidates for all k small enough to enable full k-mer counting
   *
   * @param sig_bg_c  Ouput parameter, for storing counts for final k-length
   * mers
   */
  void counted_seeds_and_candidates(gapmer_count_type& sig_bg_c) {
    std::cerr << "Lookup tables up to " << int(lookup_k_) << std::endl;
    // Initialize by counting 5-mers
    gapmer_count_type sig_bg_a(sig_path_, bg_path_, 5, *this);
    if constexpr (enable_smootihing) {
      sig_bg_a.smooth();
    }

    for (uint8_t k = 6; k <= lookup_k_; ++k) {
      // Initialize the k + 1 to compute extensions
      gapmer_count_type sig_bg_b(sig_path_, bg_path_, k, *this);
      if constexpr (enable_smootihing) {
        sig_bg_b.smooth();
      }
      uint64_t v_lim = gapmer_count_type::ONE << ((k - 1) * 2);
#pragma omp parallel for
      for (uint64_t v = 0; v < v_lim; ++v) {
        check_count(k - 1, v, 0, 0, 0, sig_bg_a, sig_bg_b);
      }
      uint8_t gap_s = middle_gap_only ? (k - 1) / 2 : 1;
      uint8_t gap_lim = middle_gap_only ? k - gap_s - 1 : k - 2;
      for (; gap_s <= gap_lim; ++gap_s) {
        for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
          uint64_t offset = sig_bg_a.offset(gap_s, gap_l);
#ifdef DEBUG
          if (offset >= gapmer_count_type::lookup_elems(k + 1)) {
            std::cerr << "invalid offset " << offset << " for gap start "
                      << int(gap_s) << " and gap length " << int(gap_l)
                      << std::endl;
            exit(1);
          }
#endif
#pragma omp parallel for
          for (uint64_t v = 0; v < v_lim; ++v) {
            check_count(k - 1, v, gap_s, gap_l, offset, sig_bg_a, sig_bg_b);
          }
        }
      }

      // k -> k + 1
      std::swap(sig_bg_a, sig_bg_b);
      std::cerr << int(k - 1) << " -> " << seeds_.size() << " candidates."
                << std::endl;
    }
    std::swap(sig_bg_a, sig_bg_c);
  }

  /**
   * Check collected mers for extension validity
   *
   * Newly found valid extensions from a to elements found in p_counter are
   * added to b
   *
   * @tparam M Type of the k-mer maps
   * @tparam P Type of the partial k-mer count data structure
   *
   * @param a  Shorter mers to extend
   * @param b  Valid extended mers
   * @param p_counter  Partial k-mer counts that may contain valid mer
   * extensions
   */
  template <class M, class P>
  void extend_counted(M& a, M& b, P& p_counter) {
    for (auto p : a) {
#ifdef DEBUG
      std::cerr << "        " << p.first.to_string() << ": "
                << p.second.sig_count << ", " << p.second.bg_count << ", "
                << p.second.p << std::endl;
#endif
      auto callback = [&](gapmer_type o) {
        if (not o.is_canonical()) {
          o = o.reverse_complement();
        }
        if (not b.contains(o)) {
          double o_a, o_b;
          if constexpr (enable_smootihing) {
            auto sig_bg = p_counter.smooth_count(o);
            o_a = sig_bg.first;
            o_b = sig_bg.second;
          } else {
            auto sig_bg = p_counter.count(o);
            o_a = sig_bg.first;
            o_b = sig_bg.second;
          }
          o_a += 1;
          o_b += 1;
          if (o_a * sig_size_ <= fold_lim_ * o_b * bg_size_) {
            return;
          }
          double o_r = error_suppressed_beta_inc(o_a, o_b, x_);
          if (do_extend<false>(p.first, o, p.second.sig_count,
                               p.second.bg_count, o_a, o_b, p.second.p, o_r)) {
            b[o] = {o, o_r, uint64_t(o_a), uint64_t(o_b)};
          }
        }
      };
      p.first.template huddinge_neighbours<true, true, false>(callback);
    }
  }

  /**
   * Find k length extensions from a, and store valid extensions in b.
   *
   * @tparam M  K-mer map type.
   * @tparam P  Partial k-mer counting data structure type.
   *
   * @param a  k - 1 length mers to extend.
   * @param b  storage for valid k-mers.
   * @param p_counter  Structure to use with partial counting.
   * @param k  length of k-mers to find.
   * @param prune  Should only one pass of extensions be done.
   */
  template <class M, class P>
  void extend(M& a, M& b, P& p_counter, uint16_t k, bool prune) {
    const constexpr double fill_limit = 0.4;
    std::cerr << "    Extend " << a.size() << " mers." << std::endl;
    auto hash = [](const gapmer_type g) { return uint64_t(g); };
    std::unordered_set<gapmer_type, decltype(hash)> del_set;
    auto callback = [&](gapmer_type o) {
      if (not b.contains(o)) {
        p_counter.init(o);
        if constexpr (enable_smootihing) {
          auto ccb = [&](gapmer_type h_n) { p_counter.init(h_n); };
          o.hamming_neighbours(ccb);
        }
      }
    };
    if (prune) {
      std::vector<Res> prio;
      for (auto p : a) {
        prio.push_back(p.second);
      }
      std::sort(prio.begin(), prio.end(), [](const auto& lhs, const auto& rhs) {
        return (lhs.sig_count / lhs.bg_count) > (rhs.sig_count / rhs.bg_count);
      });
      for (auto km : prio) {
        km.g.template huddinge_neighbours<true, true, false>(callback);
        if (p_counter.fill_rate() >= fill_limit) {
          break;
        }
      }
    } else {
      for (auto p : a) {
        p.first.template huddinge_neighbours<true, true, false>(callback);
        if (p_counter.fill_rate() >= fill_limit) {
          std::cerr << "\tLoad factor >= " << fill_limit << " ("
                    << p_counter.fill_rate() << ") counting.." << std::endl;
          p_counter.template count_mers<middle_gap_only, max_gap>(sig_path_,
                                                                  bg_path_, k);
          std::cerr << "\tFiltering extension..." << std::endl;
          extend_counted(a, b, p_counter);
          p_counter.clear();
        }
      }
    }
    std::cerr << "\tFinal load factor " << p_counter.fill_rate()
              << " counting.." << std::endl;
    p_counter.template count_mers<middle_gap_only, max_gap>(sig_path_, bg_path_,
                                                            k);
    std::cerr << "\tFiltering extension..." << std::endl;
    extend_counted(a, b, p_counter);
    p_counter.clear();
    if constexpr (filter_mers) {
      std::cerr << "\tFiltering sources..." << std::endl;
      for (auto p : a) {
#ifdef DEBUG
        std::cerr << "        " << p.first.to_string() << ": "
                  << p.second.sig_count << ", " << p.second.bg_count << ", "
                  << p.second.p << std::endl;
#endif
        bool keep = true;
        auto callback_extend = [&](gapmer_type o) {
          if (not o.is_canonical()) {
            o = o.reverse_complement();
          }
          if (b.contains(o)) {
            if (do_extend(p.first, o, p.second.sig_count, p.second.bg_count,
                          b[o].sig_count, b[o].bg_count, p.second.p, b[o].p)) {
              if constexpr (filter_mers) {
                keep = false;
              }
            }
          }
        };
        p.first.template huddinge_neighbours<true, true, false>(
            callback_extend);
        if constexpr (filter_mers) {
          if (not keep) {
            del_set.insert(p.first);
          }
        }
      }
      for (auto d : del_set) {
        a.erase(d);
      }
    }
  }

  /**
   * Filter found mers of the same lengths. Keep only the best within H1
   * distance among same length gapmers.
   *
   * @tparam M  gapmer map type
   * @param m   gapmers to filter
   */
  template <class M>
  void filter(M& m) {
    auto hash = [](const gapmer_type& g) { return uint64_t(g); };
    std::unordered_set<gapmer_type, decltype(hash)> del_set;
    auto e = m.end();
    for (auto it = m.begin(); it != e; ++it) {
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
      m.erase(d);
    }
    std::cerr << "    filtered to " << m.size() << " mers" << std::endl;
  }

 public:
  seed_finder(const std::string& sig_path, const std::string& bg_path, double p,
              double log_fold = 0.5, uint8_t max_k = 10,
              double memory_limit = 4, double p_ext = 0.01,
              uint8_t lookup_k = 10, bool prune = false)
      : sig_path_(sig_path),
        bg_path_(bg_path),
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
    gsl_set_error_handler_off();
    seq_io::Reader_x sr(sig_path_);
    sr.enable_reverse_complements();
    while (true) {
      uint64_t len = sr.get_next_read_to_buffer();
      if (len == 0) {
        break;
      }
      sig_size_ += len;
    }
    seq_io::Reader_x br(bg_path_);
    br.enable_reverse_complements();
    while (true) {
      uint64_t len = br.get_next_read_to_buffer();
      if (len == 0) {
        break;
      }
      bg_size_ += len;
    }
    x_ = double(sig_size_) / (sig_size_ + bg_size_);
    std::cerr << "Background " << bg_path_ << " with length " << bg_size_
              << std::endl;
    std::cerr << "Signal " << sig_path_ << " with length " << sig_size_
              << std::endl;
    std::cerr << "X = " << x_ << std::endl;
  }

  /**
   * Does the heavy lifting of counting increasingly long gapmers to find seed
   * candidates
   */
  void find_seeds() {
    auto hash = [](const gapmer_type g) { return uint64_t(g); };
    std::unordered_map<gapmer_type, Res, decltype(hash)> a;
    std::unordered_map<gapmer_type, Res, decltype(hash)> b;
    // full k-mer couting as long as memory is sufficient.
    {
      gapmer_count_type sig_bg_c;
      counted_seeds_and_candidates(sig_bg_c);
      uint64_t v_lim = gapmer_count_type::ONE << (lookup_k_ * 2);
#pragma omp parallel for
      for (uint64_t v = 0; v < v_lim; ++v) {
        filter_count(lookup_k_, v, 0, 0, 0, sig_bg_c, a);
      }
      uint8_t gap_s = middle_gap_only ? lookup_k_ / 2 : 1;
      uint8_t gap_lim = middle_gap_only ? lookup_k_ - gap_s : lookup_k_ - 1;
      for (; gap_s <= gap_lim; ++gap_s) {
        for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
          uint64_t offset = sig_bg_c.offset(gap_s, gap_l);
#pragma omp parallel for
          for (uint64_t v = 0; v < v_lim; ++v) {
            filter_count(lookup_k_, v, gap_s, gap_l, offset, sig_bg_c, a);
          }
        }
      }
    }
    // Partial count with extensions when we can no longer count everything
    partial_count<gapmer_type> p_counter(*this);
    for (uint8_t k = lookup_k_ + 1; k <= k_lim_; ++k) {
      std::cerr << int(k) - 1 << " -> " << std::endl;
      extend(a, b, p_counter, k, prune_);
      std::cerr << "    " << a.size() << " " << int(k) - 1 << " candidates\n"
                << "    " << b.size() << " " << int(k) << " potentials"
                << std::endl;
      for (auto p : a) {
        seeds_.push_back(p.second);
      }
      a.clear();
      if constexpr (filter_mers) {
        filter(b);
      }
      a.swap(b);
      std::cerr << int(k) - 1 << " -> " << seeds_.size() << " candidates."
                << std::endl;
      if (a.size() == 0) {
        break;
      }
    }
    if (a.size() > 0) {
      for (auto p : a) {
        seeds_.push_back(p.second);
      }
      std::cerr << int(k_lim_) << " -> " << seeds_.size() << " candidates."
                << std::endl;
    }
  }

  const std::vector<Res>& get_seeds() const { return seeds_; }

  double x() const { return x_; }

 private:
  bool should_report_errors_for_path(libbio_reader_adapter&,
                                     std::string_view path) override {
    return !paths_with_errors_.contains(path);
  }

  void found_first_read_with_unexpected_character(
      libbio_reader_adapter&, std::string_view path,
      std::uint64_t lineno) override {
    std::cerr << "WARNING: Skipping reads with unexpected characters in "
              << path << "; first one on line " << lineno << ".\n";
  }

  void found_total_reads_with_unexpected_characters(
      libbio_reader_adapter&, std::string_view path,
      std::uint64_t count) override {
    paths_with_errors_.emplace(path);
    std::cerr << "WARNING: Skipped " << count << " reads in " << path << ".\n";
  }
};
}  // namespace sf
