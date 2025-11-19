#pragma once

#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <libbio/algorithm.hh>
#include <string>
#include <utility>

#include "libbio_reader_adapter.hpp"

namespace sf {

template <class gapmer_t>
class partial_count {
 private:
  static const constexpr uint64_t MOD = 59999999;
  static const constexpr uint64_t STEP = 127;
  struct elem {
    gapmer_t mer;
    uint64_t sig_count;
    uint64_t bg_count;
  };

  libbio_reader_adapter reader_adapter;
  elem* data_;
  uint64_t size_;

  uint64_t find(gapmer_t g) const {
    uint64_t g_val(g);
    uint64_t idx = g_val % MOD;
    uint64_t val(data_[idx].mer);
    while (val > 0 && val != g_val) {
      idx += STEP;
      idx -= idx >= MOD ? MOD : 0;
      val = data_[idx].mer;
    }
    return idx;
  }

 public:
  explicit partial_count(libbio_reader_adapter_delegate& delegate)
      : reader_adapter(delegate), size_(0) {
    data_ = (elem*)calloc(MOD, sizeof(elem));
  }

  ~partial_count() {
    if (data_ != nullptr) {
      free(data_);
    }
  }

  void init(gapmer_t g) {
    uint64_t idx = find(g);
    if (data_[idx].mer != g) {
      data_[idx] = {g, 0, 0};
      ++size_;
    }
  }

  void clear() {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wclass-memaccess"
    std::memset(data_, 0, MOD * sizeof(elem));
#pragma GCC diagnostic pop
    size_ = 0;
  }

  void inc_sig(gapmer_t g) {
    uint64_t idx = find(g);
    if (data_[idx].mer == g) {
#pragma omp atomic
      data_[idx].sig_count = data_[idx].sig_count + 1;
    }
  }

  void inc_bg(gapmer_t g) {
    uint64_t idx = find(g);
    if (data_[idx].mer == g) {
#pragma omp atomic
      data_[idx].bg_count = data_[idx].bg_count + 1;
    }
  }

  template <bool middle_gap_only, uint16_t max_gap>
  void count_mers(const std::string& sig_path, const std::string& bg_path,
                  uint16_t k) {
    auto const do_count([&](std::string const& fasta_path, auto&& inc) {
      reader_adapter.read_from_path(fasta_path);

#pragma omp parallel
      while (true) {
        bool should_continue{};
#pragma omp critical
        should_continue = reader_adapter.retrieve_next_read();

        if (!should_continue) break;

        if (reader_adapter.read_length() < k) continue;

        auto const& read_buffer{reader_adapter.read_buffer()};
        gapmer_t g(read_buffer.front() >> (64U - 2U * k), k);
        inc(g);
        reader_adapter.iterate_characters(k, [&](std::uint8_t const cc) {
          g = g.next_(cc);
          inc(g);
        });

        uint8_t gap_s = middle_gap_only ? k / 2 : 1;
        auto gap_lim = libbio::min_ct(reader_adapter.read_length(),
                                      middle_gap_only ? (k + 3) / 2 : k);

        for (; gap_s < gap_lim; ++gap_s) {
          for (uint8_t gap_l = 1; gap_l <= max_gap; ++gap_l) {
            auto const tail_start{gap_s + gap_l};
            g = gapmer_t(read_buffer, k, gap_s, gap_l);
            inc(g);
            reader_adapter.iterate_character_pairs(
                gap_s, tail_start,
                [&](std::uint8_t const lhsc, std::uint8_t const rhsc) {
                  g = g.next_(lhsc, rhsc);
                  inc(g);
                });
          }
        }
      }
      reader_adapter.finish();
    });

    do_count(sig_path, [this](auto const gg) { inc_sig(gg); });
    do_count(bg_path, [this](auto const gg) { inc_bg(gg); });
  }

  std::pair<uint64_t, uint64_t> count(gapmer_t g) const {
    uint64_t idx = find(g);
    if (data_[idx].mer == g) {
      return {data_[idx].sig_count, data_[idx].bg_count};
    }
    return {0, 0};
  }

  std::pair<double, double> smooth_count(gapmer_t g) const {
    uint64_t idx = find(g);
    double sig_val = data_[idx].sig_count;
    double bg_val = data_[idx].bg_count;
    std::array<uint64_t, 10> addables{};
    std::array<uint64_t, 10> bddables{};
    uint16_t smallest = 0;
    auto callback = [&](gapmer_t o) {
      uint64_t o_idx = find(o);
      uint64_t o_sig = data_[o_idx].sig_count;
      if (addables[smallest] < o_sig) {
        addables[smallest] = o_sig;
        bddables[smallest] = data_[o_idx].bg_count;
        smallest = 0;
        for (uint16_t i = 1; i < 10; ++i) {
          smallest = addables[i] < addables[smallest] ? i : smallest;
        }
      }
    };
    g.hamming_neighbours(callback);
    for (uint16_t i = 0; i < 10; ++i) {
      sig_val += addables[i];
      bg_val += bddables[i];
    }
    return {sig_val / 11, bg_val / 11};
  }

  double fill_rate() const { return double(size_) / MOD; }
};

}  // namespace sf
