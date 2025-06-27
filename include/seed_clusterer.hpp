#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace sf {

template <bool middle_gap_only, uint16_t max_gap, class Res_vec_T>
class seed_clusterer {
 private:
  const Res_vec_T& seeds_;
  std::vector<std::vector<size_t>> edges_;
  std::vector<std::pair<size_t, double>> local_optima_;
  std::string sig_path_;
  std::string bg_path_;
  size_t opt_idx_;
  double p_ext_;

 public:
  seed_clusterer(const Res_vec_T& seeds, const std::string sig_path,
                 const std::string bg_path, double pext)
      : seeds_(seeds),
        edges_(seeds.size()),
        local_optima_(),
        sig_path_(sig_path),
        bg_path_(bg_path),
        opt_idx_(),
        p_ext_(pext) {
    std::cerr << "\tBuilding Huddinge graph" << std::endl;
    std::unordered_map<uint64_t, size_t> seed_map;
    for (size_t i = 0; i < seeds_.size(); ++i) {
      seed_map[uint64_t(seeds_[i].g)] = i;
    }
#pragma omp parallel for
    for (size_t i = 0; i < seeds_.size(); ++i) {
      auto mer = seeds_[i].g;
      auto callback = [&](const auto& o) {
        auto o_mer = o.is_canonical() ? o : o.reverse_complement();
        if (seed_map.contains(uint64_t(o_mer))) {
          size_t idx = seed_map[uint64_t(o_mer)];
#pragma omp critical
          {
            edges_[i].push_back(idx);
            edges_[idx].push_back(i);
          }
        }
      };
      mer.huddinge_neighbours(callback);
    }
    /*
    for (size_t i = 0; i < seeds_.size(); ++i) {
      std::cout << seeds_[i].g.to_string() << "\t[";
      size_t lim = edges_[i].size();
      for (size_t j = 0; j < lim; ++j) {
        std::cout << edges_[i][j] << (j == lim  - 1 ? "" : ",");
      }
      std::cout << "]\t(" << i << "," << seeds_.size() << ")" << std::endl;
    }
    */
    std::cerr << "\tComputing priorities" << std::endl;
#pragma omp parallel for
    for (size_t i = 0; i < seeds_.size(); ++i) {
      bool is_opt = true;
      double val = 0;
      const auto& out_edges = edges_[i];
      auto g = seeds_[i].g;
      auto g_len = g.length();
      auto g_sig = seeds_[i].sig_count;
      for (auto e : out_edges) {
        if (seeds_[e].p <= seeds_[i].p) {
          is_opt = false;
          break;
        }
        auto o_len = seeds_[e].g.length();
        auto o_sig = seeds_[e].sig_count;
        if (o_len < g_len) {
          double p_b = gsl_cdf_binomial_Q(g_sig, 0.25, o_sig);
          if (p_b > p_ext_) {
            is_opt = false;
            break;
          }
          val = std::max(val, p_b);
        } else if (o_len > g_len) {
          double p_b = gsl_cdf_binomial_P(o_sig, 0.25, g_sig);
          if (p_b > p_ext_) {
            is_opt = false;
            break;
          }
          val = std::max(val, p_b);
        } else {
          if (o_sig >= g_sig) {
            double p_b = gsl_cdf_binomial_P(g_sig, 0.5, 2 * o_sig);
            if (p_b > p_ext_) {
              is_opt = false;
              break;
            }
            val = std::max(val, p_b);
          }
        }
      }
      if (is_opt) {
        size_t in_h = 0;
        size_t ne_s = 0;
        auto cb = [&](const auto &o) {
          auto o_mer = o.is_canonical() ? o : o.reverse_complement();
          if (seed_map.contains(uint64_t(o_mer))) {
            ++in_h;
          }
          ++ne_s;
        };
        g.huddinge_neighbours(cb);
        val = 1 - double(in_h) / ne_s;
#pragma omp critical
        local_optima_.push_back({i, val});
      }
    }
    std::sort(local_optima_.begin(), local_optima_.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
  }

  bool has_next() const { return opt_idx_ < local_optima_.size(); }

  // “All matches” refers to matches that are substrings of other matches.
  void output_cluster(bool should_output_all_matches) {
    auto res = seeds_[local_optima_[opt_idx_].first];
    std::cout << res.g.to_string() << "\t(" << res.sig_count << ","
              << res.bg_count << ")\t" << res.p << "\t"
              << local_optima_[opt_idx_].second << std::endl;
    if (should_output_all_matches)
      ++opt_idx_;
    else {
      for (size_t i = ++opt_idx_; i < local_optima_.size(); ++i) {
        auto o_s = seeds_[local_optima_[i].first].g;
        if (o_s.aligns_to(res.g)) {
          local_optima_[i].second = 1;
        }
      }
    }
    for (; opt_idx_ < local_optima_.size(); ++opt_idx_) {
      if (local_optima_[opt_idx_].second < 1) {
        break;
      }
    }
  }
};
}  // namespace sf
