#pragma once

#include <SeqIO/SeqIO.hh>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "util.hpp"

namespace sf {

template <bool middle_gap_only, uint16_t max_gap, class Res_vec_T>
class seed_clusterer {
 private:
  const Res_vec_T& seeds_;
  std::vector<std::pair<size_t, double>> local_optima_;
  std::string sig_path_;
  std::string bg_path_;
  size_t opt_idx_;
  double p_ext_;
  double h1_weight_;
  double x_;

  /**
   * Looks for reads that confrom to the Hamming 1 neighbourhood
   * of `seeds_[idx]`, and outputs the multiple alignment of these
   * as fasta to `prefix + seeds_[idx].g.to_string() + '.fa'`.
   *
   * @param prefix  Folder to write to
   * @param idx     Center mer index in seeds_
   */
  void output_alignment(std::string& prefix, size_t idx) {
    const auto center = seeds_[idx].g;
    const size_t margin = 6;
    std::string m_string = center.to_string();

    size_t p_len = m_string.size();
    std::vector<std::regex> patterns;
    for (size_t i = 0; i < p_len; ++i) {
      char c = m_string[i];
      m_string[i] = '.';
      patterns.push_back(std::regex(m_string));
      m_string[i] = c;
    }
    std::vector<std::vector<size_t>> pfm = {
        std::vector<size_t>(2 * margin + m_string.size()),
        std::vector<size_t>(2 * margin + m_string.size()),
        std::vector<size_t>(2 * margin + m_string.size()),
        std::vector<size_t>(2 * margin + m_string.size())};

    /* Cout match locations to the pfm */
    seq_io::Reader_x sr(sig_path_);
    sr.enable_reverse_complements();
    while (true) {
      uint64_t len = sr.get_next_read_to_buffer();
      if (len == 0) {
        break;
      }
      std::string read(sr.read_buf, len);
      for (size_t p_i = 0; p_i < patterns.size(); ++p_i) {
        if (m_string[p_i] == '.') {
          continue;
        }
        auto e = std::sregex_iterator();
        auto m = std::sregex_iterator(read.begin(), read.end(), patterns[p_i]);
        for (; m != e; ++m) {
          size_t pos = m->position();
          size_t i = pos < 6 ? margin - pos : 0;
          for (; i < margin; ++i) {
            ++pfm[nuc_to_v[read[pos - (margin - i)]]][i];
          }
          if (m_string[p_i] != '.') {
            ++pfm[nuc_to_v[read[pos + p_i]]][margin + p_i];
          }
          for (i = 1; i < m_string.size() - 1; ++i) {
            if (m_string[i] == '.') {
              ++pfm[nuc_to_v[read[pos + i]]][margin + i];
            }
          }
          size_t e_i = 2 * margin + p_len;
          for (i = margin + p_len; i < e_i; ++i) {
            if (pos + i >= read.size()) {
              break;
            }
            ++pfm[nuc_to_v[read[pos + i]]][i];
          }
        }
      }
    }
    std::string f_name = prefix + m_string + ".pfm";
    std::ofstream out_file(f_name);
    for (auto vec : pfm) {
      for (auto v : vec) {
        out_file << " " << v;
      }
      out_file << "\n";
    }
  }

 public:
  /**
   * Compute most likely seeds in Huddinge graph implied by the input seeds
   *
   * @param seeds    Vector of seeds, with counts and p-value computed by
   *                 sf::seed_finder
   * @param sig_path Path to signal `.fast(a|q)(.gz)?`.
   * @param bg_path  Path to backgroud `.fast(a|q)(.gz)?`.
   * @param pext     P-value to filter extensions with binomial tests.
   */
  seed_clusterer(const Res_vec_T& seeds, std::string_view sig_path,
                 std::string_view bg_path, double pext, double h1_weight,
                 double x)
      : seeds_(seeds),
        local_optima_(),
        sig_path_(sig_path),
        bg_path_(bg_path),
        opt_idx_(),
        p_ext_(pext),
        h1_weight_(h1_weight),
        x_(x) {
    std::unordered_map<uint64_t, size_t> seed_map;
    for (size_t i = 0; i < seeds_.size(); ++i) {
      seed_map[uint64_t(seeds_[i].g)] = i;
    }
    /* Compute the Huddinge priorities for all potential seeds */
#pragma omp parallel for
    for (size_t i = 0; i < seeds_.size(); ++i) {
      auto mer = seeds_[i].g;
      // Value used for attempt to filter false positives
      // uint32_t bogo_ratio = seeds_[i].sig_count - seeds_[i].bg_count;
      bool keep = true;
      double enrichment =
          seeds_[i].sig_count / x_ - seeds_[i].bg_count / (1 - x_);
      uint16_t len = seeds_[i].g.length();
      double val = 0;
      size_t h_n_count = 0;
      // Priority is based on enrichment of the candidate and its H1
      // neighbourhood.
      //
      // Adds `sig_count - bg_count` for neighbours in graph, `0` for others.
      auto cb = [&](const auto& o) {
        uint64_t o_mer =
            o.is_canonical() ? uint64_t(o) : uint64_t(o.reverse_complement());
        if (seed_map.contains(o_mer)) {
          size_t idx = seed_map[o_mer];
          double o_enrichment =
              seeds_[idx].sig_count / x_ - seeds_[idx].bg_count / (1 - x_);
          if (o_enrichment > enrichment) {
            keep = false;
          }
          val += o_enrichment * ((seeds_[idx].sig_count * (1 - x_)) /
                                 (x_ * seeds_[idx].bg_count));
          ++h_n_count;
        } else {
          ++h_n_count;
        }
      };
      if (len >= 24) {
        mer.template huddinge_neighbours<false, false, true>(cb);
      } else if (len <= 5) {
        mer.template huddinge_neighbours<true, false, false>(cb);
      } else {
        mer.huddinge_neighbours(cb);
      }
      val /= h_n_count;
      val *= h1_weight_;
      double f_m =
          ((seeds_[i].sig_count * (1 - x_)) / (x_ * seeds_[i].bg_count)) * len;
      val += enrichment * f_m;
      if (keep) {
#pragma omp critical
        local_optima_.push_back({i, 1 / val});
      }
    }
    std::sort(local_optima_.begin(), local_optima_.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    std::cerr << "Computed priorities for " << local_optima_.size()
              << std::endl;
    /*for (auto e : local_optima_) {
      std::cout << seeds_[e.first].g.to_string() << "\t" <<
    seeds_[e.first].sig_count << "\t" << seeds_[e.first].bg_count << std::endl;
    }*/
  }

  /**
   * `true` if `this` has local optima to left to output.
   */
  bool has_next() const { return opt_idx_ < local_optima_.size(); }

  /**
   * Outputs the 'best' seed that has not yet been output to `std::cout`
   *
   * Optionally also outputs MSA based on the seed.
   *
   * @param prefix  Folder for outputting fasta (needs to exist).
   * @param should_output_all_matches “All matches” refers to matches that are
   *                                  substrings of other matches.
   */
  void output_cluster(std::string prefix, bool should_output_all_matches) {
    auto res = seeds_[local_optima_[opt_idx_].first];
    std::cout << res.g.to_string() << "\t(" << res.sig_count << ","
              << res.bg_count << ")\t" << res.p << "\t"
              << local_optima_[opt_idx_].second << std::endl;
    if (prefix.size() > 0) {
      if (not prefix.ends_with('/')) {
        prefix.push_back('/');
      }
      prefix.append(std::to_string(opt_idx_));
      prefix.push_back('_');
      output_alignment(prefix, local_optima_[opt_idx_].first);
    }

    /* Get index of next seed */
    if (should_output_all_matches)
      ++opt_idx_;
    else {
      /* Skip seeds that are too close to something already output.
       *
       * Seeds that are too close include seeds that are
       * * within H1 of a previous seed,
       * * substring of a previous seed,
       * * supersting of a previous seed.
       *
       * Seeds to be skipped get priority 3, which is impossibly high for a true
       * seed
       */
      for (size_t i = ++opt_idx_; i < local_optima_.size(); ++i) {
        auto o_s = seeds_[local_optima_[i].first].g;
        int out = 0;
        if (o_s.aligns_to(res.g) or res.g.aligns_to(o_s) or
            o_s.huddinge_distance(res.g, out) <= 1 or
            o_s.reverse_complement().huddinge_distance(res.g, out) <= 1) {
          local_optima_[i].second = 3;
        }
      }
    }

    /* Skip over impossible seeds */
    for (; opt_idx_ < local_optima_.size(); ++opt_idx_) {
      if (local_optima_[opt_idx_].second < 3) {
        break;
      }
    }
  }
};


// Convenience function for instantiating seed_clusterer. t_result_vector is
// deduced. (C++ as of version 23 does not have partial class template argument
// deduction and this is the second-best option.)
template <bool t_middle_gap_only, uint16_t t_max_gap, typename t_result_vector>
auto make_seed_clusterer(t_result_vector const& seeds,
                         std::string_view sig_path, std::string_view bg_path,
                         double const pext, double const h1_weight,
                         double const x)
    -> seed_clusterer<t_middle_gap_only, t_max_gap,
                      std::remove_cvref_t<t_result_vector>> {
  return {seeds, sig_path, bg_path, pext, h1_weight, x};
}

}  // namespace sf
