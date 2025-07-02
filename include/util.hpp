#pragma once
#include <array>
#include <cstdint>
#include <unordered_set>
#include <string>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>

namespace sf {
const static constexpr std::array<char, 4> v_to_nuc = {'A', 'C', 'G', 'T'};
const static constexpr std::array<uint8_t, 256> nuc_to_v = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 2,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

const static constexpr std::array<uint8_t, 256> rc_byte = []() constexpr {
  std::array<uint8_t, 256> ret;
  for (uint16_t i = 0; i < 256; ++i) {
    uint16_t ii = ~i;
    uint8_t v = (ii >> 6) & 0b11;
    v |= (ii >> 2) & 0b1100;
    v |= (ii << 2) & 0b110000;
    v |= (ii << 6) & 0b11000000;
    ret[i] = v;
  }
  return ret;
}();

template <class G, typename F, bool middle_gap_only = true, uint8_t max_gap = 6>
void gap_mer_generation(uint8_t k, F& callback) {
  uint64_t lim;
  if (k == 5) {
    lim = uint64_t(1) << (4 * 2);
    for (uint64_t v = 0; v < lim; ++v) {
      G o(v, 4);
      callback(o);
    }
    for (uint16_t gl = 1; gl <= max_gap; ++gl) {
      for (uint64_t v = 0; v < lim; ++v) {
        if constexpr (middle_gap_only) {
          G o(v, 4, 2, gl);
          callback(o);
        } else {
          for (uint16_t gs = 1; gs < 4; ++gs) {
            G o(v, 4, gs, gl);
            callback(o);
          }
        }
      }
    }
  }

  lim = uint64_t(1) << (5 * 2);
  for (uint64_t v = 0; v < lim; ++v) {
    G o(v, 5);
    callback(o);
  }
  for (uint16_t gl = 1; gl <= max_gap; ++gl) {
    for (uint64_t v = 0; v < lim; ++v) {
      if constexpr (middle_gap_only) {
        for (uint16_t gs = 2; gs <= 3; ++gs) {
          G o(v, 5, gs, gl);
          callback(o);
        }
      } else {
        for (uint16_t gs = 1; gs < 5; ++gs) {
          G o(v, 5, gs, gl);
          callback(o);
        }
      }
    }
  }

  lim = uint64_t(1) << (6 * 2);
  for (uint64_t v = 0; v < lim; ++v) {
    G o(v, 6);
    callback(o);
  }
  for (uint16_t gl = 1; gl <= max_gap; ++gl) {
    for (uint64_t v = 0; v < lim; ++v) {
      if constexpr (middle_gap_only) {
        G o(v, 6, 3, gl);
        callback(o);
      } else {
        for (uint16_t gs = 1; gs < 6; ++gs) {
          G o(v, 6, gs, gl);
          callback(o);
        }
      }
    }
  }

  if (k == 6) {
    lim = uint64_t(1) << (7 * 2);
    for (uint64_t v = 0; v < lim; ++v) {
      G o(v, 7);
      callback(o);
    }
    for (uint16_t gl = 1; gl <= max_gap; ++gl) {
      for (uint64_t v = 0; v < lim; ++v) {
        if constexpr (middle_gap_only) {
          for (uint16_t gs = 3; gs <= 4; ++gs) {
            G o(v, 7, gs, gl);
            callback(o);
          }
        } else {
          for (uint16_t gs = 1; gs < 7; ++gs) {
            G o(v, 7, gs, gl);
            callback(o);
          }
        }
      }
    }
  }
}

template <class G, typename F, bool middle_gap_only = true, uint8_t max_gap = 6>
void gap_mer_neighbour_generation(G g, F& callback) {
  auto second_callback = [&](G o) {
    if (g.is_neighbour(o)) {
      callback(o);
    }
  };
  gap_mer_generation<G, decltype(second_callback), middle_gap_only, max_gap>(
      g.length(), second_callback);
}

template <class G, bool middle_gap_only, uint16_t max_gap>
bool compare_generation(G g, std::unordered_set<std::string>& a, std::unordered_set<std::string>& b) {
  auto callback = [&](G o) { a.insert(o.to_string()); };
  gap_mer_neighbour_generation<G, decltype(callback), middle_gap_only,
                                   max_gap>(g, callback);

  auto o_callback = [&](G o) { b.insert(o.to_string()); };
  g.huddinge_neighbours(o_callback);

  return a == b;
}

template <class G, bool middle_gap_only, uint16_t max_gap>
bool compare_generation(G g) {
  std::unordered_set<std::string> a;
  std::unordered_set<std::string> b;
  return compare_generation<G, middle_gap_only, max_gap>(g, a, b);
}

inline double error_suppressed_beta_inc(double a, double b, double x) {
  gsl_sf_result res;
  int err = gsl_sf_beta_inc_e(a, b, x, &res);
  return err ? 0 : res.val;
}

}  // namespace sf
