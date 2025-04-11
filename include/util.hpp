#pragma once
#include <array>
#include <cstdint>

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

}  // namespace sf
