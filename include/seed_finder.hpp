#pragma once

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "fm_index.hpp"
#include "gapmer.hpp"
#include "gapmer_count.hpp"
#include "partial_count.hpp"
#include "util.hpp"

namespace sf {

template <bool middle_gap_only, uint8_t max_gap, bool enable_smootihing = true>
class seed_finder {
 private:
  typedef gapmer<middle_gap_only, max_gap> G;
  typedef gapmer_count<middle_gap_only, max_gap> G_C;
  struct Res {
    G g;
    double p;
    uint64_t sig_count;
    uint64_t bg_count;
  };

  std::string sig_path_;
  std::string bg_path_;
  std::vector<Res> seeds_;
  uint64_t sig_size_;
  uint64_t bg_size_;
  double p_;
  double p_ext_;
  double fold_lim_;
  double x_;
  double memory_limit_;
  uint8_t k_lim_;

  template <bool debug = false>
  bool do_extend([[maybe_unused]] G a, [[maybe_unused]] G b, double a_sig,
                 double a_bg, double b_sig, double b_bg, double a_r,
                 double b_r) {
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
#ifdef DEBUG
    if (b_r < a_r) {
      std::cerr << "        " << a.to_string() << " discarded by "
                << b.to_string() << "\n            (" << a_sig << ", " << a_bg
                << ") <-> (" << b_sig << ", " << b_bg
                << ")\n            with p " << b_r << std::endl;
    } else {
      std::cerr << "        " << a.to_string() << " discards " << b.to_string()
                << "\n            (" << a_sig << ", " << a_bg << ") <-> ("
                << b_sig << ", " << b_bg << ")\n            with p " << b_r
                << std::endl;
    }
#endif
    return b_r <= a_r;
  }

  template <bool calc_b_r>
  bool do_filter([[maybe_unused]] G a, [[maybe_unused]] G b, double a_sig,
                 double a_bg, double b_sig, double b_bg, double a_r,
                 double& b_r) {
    if (a_bg == 1 && b_bg == 1) {
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

  void check_count(const uint8_t k, const uint64_t v, uint8_t gap_s,
                   uint8_t gap_l, uint64_t offset, G_C& sig_bg_a,
                   G_C& sig_bg_b) {
#ifdef DEBUG
    if (offset + v >= G_C::lookup_elems(k)) {
      std::cerr << "accessing " << offset << " + " << v << " = " << offset + v
                << " of " << G_C::lookup_elems(k) << " element table"
                << std::endl;
      exit(1);
    }
#endif
    if (sig_bg_a.discarded[offset + v]) {
      return;
    }
    G g(v, k, gap_s, gap_l);
    if (not g.is_canonical()) {
#pragma omp critical(a_bv)
      sig_bg_a.discarded[offset + v] = true;
      return;
    }
    uint64_t a = sig_bg_a.sig_counts[offset + v] + 1;
    uint64_t b = sig_bg_a.bg_counts[offset + v] + 1;
    if (a * bg_size_ <= fold_lim_ * b * sig_size_) {
#pragma omp critical(a_bv)
      sig_bg_a.discarded[offset + v] = true;
    }
    double r = error_suppressed_beta_inc(a, b, x_);
    if (r > p_) {
#pragma omp critical(a_bv)
      sig_bg_a.discarded[offset + v] = true;
      return;
    }
    auto callback_a = [&](G o) {
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
    auto callback_b = [&](G o) {
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
  }

  template <class M>
  void filter_count(const uint8_t k, const uint64_t v, uint8_t gap_s,
                    uint8_t gap_l, uint64_t offset, G_C& sig_bg_c, M& m) {
#ifdef DEBUG
    if (offset + v >= G_C::lookup_elems(k)) {
      std::cerr << "k = " << int(k) << " & sig_bg_c.k_ = " << int(sig_bg_c.k_)
                << " :\n accessing " << offset << " + " << v << " = "
                << offset + v << " of " << G_C::lookup_elems(k)
                << " element table" << std::endl;
      exit(1);
    }
#endif
    if (sig_bg_c.discarded[offset + v]) {
      return;
    }
    G g(v, k, gap_s, gap_l);
    if (not g.is_canonical()) {
      sig_bg_c.discarded[offset + v] = true;
      return;
    }
    uint64_t a = sig_bg_c.sig_counts[offset + v] + 1;
    uint64_t b = sig_bg_c.bg_counts[offset + v] + 1;
    if (a * sig_size_ <= fold_lim_ * b * bg_size_) {
#pragma omp critical(d_bv)
      sig_bg_c.discarded[offset + v] = true;
    }
    double r = error_suppressed_beta_inc(a, b, x_);
    if (r > p_) {
#pragma omp critical(d_bv)
      sig_bg_c.discarded[offset + v] = true;
      return;
    }
    auto callback = [&](G o) {
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
  }

  uint8_t counted_seeds_and_candidates(G_C& sig_bg_c) {
    uint8_t k_lim = 5;
    while (G_C::lookup_bytes(k_lim) < memory_limit_) {
      ++k_lim;
    }
    std::cerr << "Lookup tables up to " << int(k_lim - 1) << std::endl;
    G_C sig_bg_a(sig_path_, bg_path_, 5);
    if constexpr (enable_smootihing) {
      sig_bg_a.smooth();
    }

    for (uint8_t k = 6; k < k_lim; ++k) {
      G_C sig_bg_b(sig_path_, bg_path_, k);
      if constexpr (enable_smootihing) {
        sig_bg_b.smooth();
      }
      uint64_t v_lim = G_C::ONE << ((k - 1) * 2);
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
          if (offset >= G_C::lookup_elems(k + 1)) {
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

      std::swap(sig_bg_a, sig_bg_b);
      std::cerr << int(k - 1) << " -> " << seeds_.size() << " seeds."
                << std::endl;
    }
    std::swap(sig_bg_a, sig_bg_c);
    return k_lim;
  }

  template <class M, class P>
  void extend(M& a, M& b, P& p_counter, uint16_t k) {
    std::cerr << "    Extend " << a.size() << " mers." << std::endl;
    for (auto p : a) {
      auto callback = [&](G o) {
        p_counter.init(p.first);
        auto ccb = [&](G h_n) { p_counter.init(h_n); };
        o.hamming_neighbours(ccb);
      };
      p.first.template huddinge_neighbours<true, true, false>(callback);
    }
    p_counter.template count_mers<middle_gap_only, max_gap>(sig_path_, bg_path_,
                                                            k);
    auto hash = [](const G g) { return uint64_t(g); };
    std::unordered_set<G, decltype(hash)> del_set;
    for (auto p : a) {
#ifdef DEBUG
      std::cerr << "        " << p.first.to_string() << ": "
                << p.second.sig_count << ", " << p.second.bg_count << ", "
                << p.second.p << std::endl;
#endif
      bool keep = true;
      auto callback = [&](G o) {
        if (not o.is_canonical()) {
          o = o.reverse_complement();
        }
        if (b.contains(o)) {
          if (do_extend(p.first, o, p.second.sig_count, p.second.bg_count,
                        b[o].sig_count, b[o].bg_count, p.second.p, b[o].p)) {
            keep = false;
          } else {
            b[o] = {p.first, 1.0, 1, 1};
          }
        } else {
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
            keep = false;
          }
        }
      };
      p.first.template huddinge_neighbours<true, true, false>(callback);
      if (not keep) {
        del_set.insert(p.first);
      }
    }
    for (auto d : del_set) {
      a.erase(d);
    }
    p_counter.clear();
  }

  template <class M>
  void filter(M& m) {
    auto hash = [](const G& g) { return uint64_t(g); };
    std::unordered_set<G, decltype(hash)> del_set;
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
              double memory_limit = 4, double p_ext = 0.01)
      : sig_path_(sig_path),
        bg_path_(bg_path),
        seeds_(),
        sig_size_(),
        bg_size_(),
        p_(p),
        p_ext_(p_ext),
        fold_lim_(std::pow(2, log_fold)),
        memory_limit_(memory_limit),
        k_lim_(max_k) {
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

  void find_seeds() {
    auto hash = [](const G g) { return uint64_t(g); };
    std::unordered_map<G, Res, decltype(hash)> a;
    std::unordered_map<G, Res, decltype(hash)> b;
    uint8_t k;
    {
      G_C sig_bg_c;
      k = counted_seeds_and_candidates(sig_bg_c);
      uint64_t v_lim = G_C::ONE << (k * 2 - 2);
#pragma omp parallel for
      for (uint64_t v = 0; v < v_lim; ++v) {
        filter_count(k - 1, v, 0, 0, 0, sig_bg_c, a);
      }
      uint8_t local_k = k - 1;
      uint8_t gap_s = middle_gap_only ? local_k / 2 : 1;
      uint8_t gap_lim = middle_gap_only ? local_k - gap_s : local_k - 1;
      for (; gap_s <= gap_lim; ++gap_s) {
        for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
          uint64_t offset = sig_bg_c.offset(gap_s, gap_l);
#pragma omp parallel for
          for (uint64_t v = 0; v < v_lim; ++v) {
            filter_count(local_k, v, gap_s, gap_l, offset, sig_bg_c, a);
          }
        }
      }
    }
    partial_count<G> p_counter;
    for (; k <= k_lim_; ++k) {
      std::cerr << int(k) - 1 << " -> " << std::endl;
      extend(a, b, p_counter, k);
      std::cerr << "    " << a.size() << " " << int(k) - 1 << " seeds\n"
                << "    " << b.size() << " " << int(k) << " candidates"
                << std::endl;
      for (auto p : a) {
        seeds_.push_back(p.second);
      }
      a.clear();
      filter(b);
      a.swap(b);
      std::cerr << int(k) - 1 << " -> " << seeds_.size() << " seeds."
                << std::endl;
      if (a.size() == 0) {
        break;
      }
    }
    if (a.size() > 0) {
      for (auto p : a) {
        seeds_.push_back(p.second);
      }
      std::cerr << int(k_lim_) << " -> " << seeds_.size() << " seeds."
                << std::endl;
    }
  }

  const std::vector<Res>& get_seeds() const { return seeds_; }
};
}  // namespace sf
