#pragma once

#include <SeqIO/SeqIO.hh>
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <regex>
#include <unordered_map>
#include <unordered_set>

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

  void output_alignment(std::string& prefix, size_t idx) {
    auto center = seeds_[idx].g;
    std::vector<std::pair<std::string, int>> mers;
    mers.push_back({center.to_string(), 0});
    prefix.append(mers.back().first).append(".txt");
    int a_offset = 0;
    int b_offset = 0;
    int min_offset = 0;
    for (auto e : edges_[idx]) {
      auto g = seeds_[e].g;
      auto g_r = g.reverse_complement();
      size_t a_d = center.huddinge_distance(g, a_offset);
      size_t b_d = center.huddinge_distance(g_r, b_offset);
      if (a_d < b_d) {
        min_offset = std::min(min_offset, a_offset);
        mers.push_back({g.to_string(), a_offset});
      } else {
        min_offset = std::min(min_offset, b_offset);
        mers.push_back({g_r.to_string(), b_offset});
      }
    }
    std::string reg = "(";
    for (auto p : mers) {
      if (reg.size() > 2) {
        reg.push_back('|');
      }
      for (int i = p.second; i > min_offset; --i) {
        reg.push_back('.');
      }
      reg.append(p.first);
    }
    reg.push_back(')');
    std::regex rex(reg);

    std::ofstream out_file(prefix);
    seq_io::Reader_x sr(sig_path_);
    sr.enable_reverse_complements();
    std::vector<std::pair<std::string, size_t>> hits;
    while (true) {
      uint64_t len = sr.get_next_read_to_buffer();
      if (len == 0) {
        break;
      }
      std::string read(sr.read_buf, len);
      auto e = std::sregex_iterator();
      auto m = std::sregex_iterator(read.begin(), read.end(), rex);
      for (; m != e; ++m) {
        hits.push_back({read, m->position()});
      }
    }
    size_t max_pos = 0;
    for (auto p : hits) {
      max_pos = std::max(max_pos, p.second);
    }
    size_t max_len = 0;
    for (auto p : hits) {
      size_t n_count = max_pos - p.second;
      max_len = std::max(n_count + p.first.size(), max_len);
    }
    size_t c = 1;
    for (auto p : hits) {
      out_file << "> " << c++ << "\n";
      size_t n_count = max_pos - p.second;
      for (size_t i = 0; i < n_count; ++i) {
        out_file << "N";
      }
      out_file << p.first;
      n_count = max_len - p.first.size() - n_count;
      for (size_t i = 0; i < n_count; ++i) {
        out_file << "N";
      }
      out_file << std::endl;
    }
  }

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
      const auto g = seeds_[i].g;
      const size_t p_sig = seeds_[i].sig_count;
      const auto& edges = edges_[i];
      for (auto e : edges) {
        const auto o_g = seeds_[e].g;
        const auto o_sig = seeds_[e].sig_count;
        if (o_g.length() > g.length() and o_sig * 4 > p_sig) {
          double binom_p = gsl_cdf_binomial_P(p_sig, 0.25, o_sig);
          if (binom_p < pext) {
            is_opt = false;
            break;
          }
        } else if (o_g.length() < g.length() and o_sig > p_sig * 4) {
          double binom_p = gsl_cdf_beta_Q(o_sig, 0.25, p_sig);
          if (binom_p < pext) {
            is_opt = false;
            break;
          }
        } else {
          if (o_sig > p_sig) {
            is_opt = false;
            break;
          }
        }
      }
      if (not is_opt) {
        continue;
      }
      size_t ne_s = 0;
      seeds_[i].g.huddinge_neighbours(
          [&]([[maybe_unused]] const auto& o) { ++ne_s; });
      double val = 1 - double(edges.size()) / ne_s;
#pragma omp critical
      local_optima_.push_back({i, val});
    }
    std::sort(local_optima_.begin(), local_optima_.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
  }

  bool has_next() const { return opt_idx_ < local_optima_.size(); }

  void output_cluster(std::string prefix = "") {
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
    for (size_t i = ++opt_idx_; i < local_optima_.size(); ++i) {
      auto o_s = seeds_[local_optima_[i].first].g;
      int out = 0;
      if (o_s.aligns_to(res.g) or res.g.aligns_to(o_s) or
          o_s.huddinge_distance(res.g, out) <= 1 or
          o_s.reverse_complement().huddinge_distance(res.g, out) <= 1) {
        local_optima_[i].second = 1;
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
