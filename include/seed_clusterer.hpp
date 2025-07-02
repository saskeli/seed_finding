#pragma once

#include <SeqIO/SeqIO.hh>
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <regex>
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

  /**
   * Looks for reads that confrom to the Huddinge 1 neighbourhood
   * of `seeds_[idx]`, and outputs the multiple alignment of these
   * as fasta to `prefix + seeds_[idx].g.to_string() + '.fa'`.
   *
   * @param prefix  Folder to write to
   * @param idx     Center mer index in seeds_
   */
  void output_alignment(std::string& prefix, size_t idx) {
    auto center = seeds_[idx].g;
    std::vector<std::pair<std::string, int>> mers;
    mers.push_back({center.to_string(), 0});
    prefix.append(mers.back().first).append(".fa");
    int a_offset = 0;
    int b_offset = 0;
    int min_offset = 0;
    /* Collects all H1 neighbours and aligns them relative to the center */
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
    /* Build a regex from all the seeds. */
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

    /* Collect all reads that match the regex */
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
    /* Compute MSA for the reads */
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
    /* Output MSA to file */
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
  /**
   * Compute most likely seeds in Huddinge graph implied by the input seeds
   *
   * @param seeds    Vector of seeds, with counts and p-value computed by
   *                 sf::seed_finder
   * @param sig_path Path to signal `.fast(a|q)(.gz)?`.
   * @param bg_path  Path to backgroud `.fast(a|q)(.gz)?`.
   * @param pext     P-value to filter extensions with binomial tests.
   */
  seed_clusterer(const Res_vec_T& seeds, const std::string sig_path,
                 const std::string bg_path, double pext)
      : seeds_(seeds),
        edges_(seeds.size()),
        local_optima_(),
        sig_path_(sig_path),
        bg_path_(bg_path),
        opt_idx_(),
        p_ext_(pext) {
    std::unordered_map<uint64_t, size_t> seed_map;
    for (size_t i = 0; i < seeds_.size(); ++i) {
      seed_map[uint64_t(seeds_[i].g)] = i;
    }
    /* Compute the Huddinge graph and priorities for all potential seeds */
#pragma omp parallel for
    for (size_t i = 0; i < seeds_.size(); ++i) {
      auto mer = seeds_[i].g;
      size_t len = mer.length();
      double val = 0;
      size_t neigbour_count = 0;
      /*
       * Priority is based on enrichment of seed and H1 neighbourhood.
       *
       * Adds `bg_count / seed_count` for neighbours in graph, `1` for others.
       */
      auto cb = [&](const auto& o) {
        uint64_t o_mer =
            o.is_canonical() ? uint64_t(o) : uint64_t(o.reverse_complement());
        if (seed_map.contains(o_mer)) {
          size_t idx = seed_map[o_mer];
          edges_[i].push_back(idx);
          val += double(seeds_[idx].bg_count) / seeds_[idx].sig_count;
        } else {
          val += 1;
        }
        ++neigbour_count;
      };
      if (len >= 24) {
        mer.template huddinge_neighbours<false, false, true>(cb);
      } else if (len <= 5) {
        mer.template huddinge_neighbours<true, false, false>(cb);
      } else {
        mer.huddinge_neighbours(cb);
      }
      val = val / neigbour_count;
      val += double(seeds_[i].bg_count) / seeds_[i].sig_count;
#pragma omp critical
      local_optima_.push_back({i, val});
    }
    std::sort(local_optima_.begin(), local_optima_.end(),
              [](const auto& a, const auto& b) { return a.second < b.second; });
    std::cerr << "Computed priorities for " << local_optima_.size()
              << std::endl;
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
}  // namespace sf
