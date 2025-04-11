#include <cstdint>
#include <generator>
#include <iostream>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <vector>

#include "include/gapmer.hpp"
#include "include/util.hpp"

void help(const char* call) {
  std::cout << R"(
Explore huddinge neighbourhood

Usage: )" << call
            << R"( [OPTION]... [input_mer]

input_mer  Possibly gapped input mer consisting of \"[ACGT.]*\".
           If given, works with Huddinge 1 distance of input mer.
           Else works with all 5- and 6-mers (requires -c).
-c         Compare single gapped mer generation techinques.
-a         Single gap can be in any postions.
           Defaults to only consider middle-gapped mers.
-h         Print usage information and terminate.
           Overrides other flags.
-s         Use string-based neighbour generation.
           (Does not impact multicomparison)
           (Will ignore gapping rules when not comparing)

)" << std::endl;
  exit(0);
}

const uint16_t max_gap = 6;

bool single_gapped(const std::string& input_mer) {
  bool gap_foud = false;
  for (uint64_t i = 1; i < input_mer.size(); ++i) {
    if (input_mer[i] == '.' && input_mer[i - 1] != '.') {
      if (gap_foud) {
#ifdef DEBUG
        std::cout << "  multigapped " << input_mer << " discarded."
                  << std::endl;
#endif
        return false;
      } else {
        gap_foud = true;
      }
    }
  }
  return true;
}

bool middle_gapped(const std::string& input_mer) {
  if (not single_gapped(input_mer)) {
#ifdef DEBUG
    std::cout << "  multigapped " << input_mer << " discarded." << std::endl;
#endif
    return false;
  }
  uint16_t pk = 0;
  for (uint64_t i = 0; i < input_mer.size(); ++i) {
    if (input_mer[i] != '.') {
      ++pk;
    } else {
      break;
    }
  }
  uint16_t sk = 0;
  for (uint64_t i = input_mer.size() - 1; i < input_mer.size(); --i) {
    if (input_mer[i] != '.') {
      ++sk;
    } else {
      break;
    }
  }
  bool ret = sk == pk || sk - 1 == pk || sk + 1 == pk;
#ifdef DEBUG
  if (not ret) {
    std::cout << "  non-middle gapped " << input_mer << " discarded."
              << std::endl;
  }
#endif
  return ret;
}

bool atmost_k_gapped(const std::string& input_mer, uint16_t k) {
  uint16_t gc = 0;
  for (auto c : input_mer) {
    gc += c == '.';
  }
  bool ret = gc <= k;
#ifdef DEBUG
  if (not ret) {
    std::cout << "  overgapped " << input_mer << " discarder" << std::endl;
  }
#endif
  return ret;
}

void string_gen(const std::string& input_mer, auto callback) {
  // For defined bases + 1
  // Add one nuc to each gap (implicit gap at start and end)

  // Start:
  for (auto n : sf::v_to_nuc) {
    std::string ne = "";
    ne.push_back(n);
    ne.append(input_mer);
    callback(ne);
    for (uint16_t gl = 1; gl <= max_gap; ++gl) {
      ne = "";
      ne.push_back(n);
      ne.append(std::string(gl, '.'));
      ne.append(input_mer);
      callback(ne);
    }
  }

  // Middle:
  for (size_t i = 0; i < input_mer.size(); ++i) {
    if (input_mer[i] == '.') {
      for (auto n : sf::v_to_nuc) {
        std::string ne = input_mer.substr(0, i);
        ne.push_back(n);
        ne.append(input_mer.substr(i + 1));
        callback(ne);
      }
    }
  }

  // End:
  for (auto n : sf::v_to_nuc) {
    std::string ne = input_mer;
    ne.push_back(n);
    callback(ne);
    for (uint16_t gl = 1; gl <= max_gap; ++gl) {
      ne = input_mer;
      ne.append(std::string(gl, '.'));
      ne.push_back(n);
      callback(ne);
    }
  }

  // For defined bases + 0
  // Hamming distance in defined bases

  // Hamming
  for (auto n : sf::v_to_nuc) {
    for (size_t i = 0; i < input_mer.size(); ++i) {
      if (input_mer[i] != '.' and input_mer[i] != n) {
        std::string ne = input_mer.substr(0, i);
        ne.push_back(n);
        ne.append(input_mer.substr(i + 1));
        callback(ne);
      }
    }
  }

  // Add one to any gap and gap any existing position
  // Start:
  for (auto n : sf::v_to_nuc) {
    for (size_t i = 0; i < input_mer.size(); ++i) {
      if (input_mer[i] != '.') {
        std::string ne = "";
        ne.push_back(n);
        ne.append(input_mer.substr(0, i));
        if (i < input_mer.size() - 1) {
          ne.push_back('.');
          ne.append(input_mer.substr(i + 1));
        }
        callback(ne);
        for (uint16_t gl = 1; gl <= max_gap; ++gl) {
          ne = "";
          ne.push_back(n);
          ne.append(std::string(gl, '.'));
          ne.append(input_mer.substr(0, i));
          if (i < input_mer.size() - 1) {
            ne.push_back('.');
            ne.append(input_mer.substr(i + 1));
          }
          callback(ne);
        }
      }
    }
  }

  // Middle:
  for (size_t i = 0; i < input_mer.size(); ++i) {
    if (input_mer[i] == '.') {
      continue;
    }
    for (size_t j = 0; j < input_mer.size(); ++j) {
      if (input_mer[j] == '.') {
        for (auto n : sf::v_to_nuc) {
          std::string ne = input_mer.substr(0, j);
          ne.push_back(n);
          ne.append(input_mer.substr(j + 1));

          std::string nne = ne.substr(0, i);
          if (i != 0 and i < input_mer.size() - 1) {
            nne.push_back('.');
          }
          nne.append(ne.substr(i + 1));
          callback(nne);
        }
      }
    }
  }

  // End:
  for (auto n : sf::v_to_nuc) {
    for (size_t i = 0; i < input_mer.size(); ++i) {
      if (input_mer[i] != '.') {
        std::string ne = input_mer.substr(0, i);
        if (i > 0) {
          ne.push_back('.');
        }
        ne.append(input_mer.substr(i + 1));
        ne.push_back(n);
        callback(ne);
        for (uint16_t gl = 1; gl <= max_gap; ++gl) {
          ne = input_mer.substr(0, i);
          if (i > 0) {
            ne.push_back('.');
          }
          ne.append(input_mer.substr(i + 1));
          ne.append(std::string(gl, '.'));
          ne.push_back(n);
          callback(ne);
        }
      }
    }
  }

  // For defined bases - 1
  // gap any existing nuc
  for (size_t i = 0; i < input_mer.size(); ++i) {
    if (input_mer[i] != '.') {
      std::string ne = input_mer.substr(0, i);
      if (i != 0 and i < input_mer.size() - 1) {
        ne.push_back('.');
      }
      ne.append(input_mer.substr(i + 1));
      callback(ne);
    }
  }
}

template <class G, bool middle_gap_only>
bool compare_single(G g) {
  std::unordered_set<std::string> a;
  auto callback = [&](G o) { a.insert(o.to_string()); };
  sf::gap_mer_neighbour_generation<G, decltype(callback), middle_gap_only,
                                   max_gap>(g, callback);

  std::unordered_set<std::string> b;
  auto o_callback = [&](G o) { b.insert(o.to_string()); };
  g.huddinge_neighbours(o_callback);

  std::vector<std::string> m_g;
  for (auto a_e : a) {
    if (not b.contains(a_e)) {
      m_g.push_back(a_e);
    }
  }

  std::vector<std::string> m_b;
  for (auto b_e : b) {
    if (not a.contains(b_e)) {
      m_b.push_back(b_e);
    }
  }

  if (m_g.size() || m_b.size()) {
    std::cout << "Problem when generating " << g.to_string()
              << " neighbours: " << std::endl;
    std::cout << "Missing from neighbour generation" << std::endl;
    for (auto s : m_g) {
      std::cout << s << std::endl;
    }
    std::cout << "Missing from brute force" << std::endl;
    for (auto s : m_b) {
      std::cout << s << std::endl;
    }
    return false;
  }
  return true;
}

template <class G, bool middle_gap_only>
bool comp_str(G g, std::string& input_mer) {
  std::unordered_set<std::string> a;
  if constexpr (middle_gap_only) {
    auto callback = [&](std::string& o) {
      if (atmost_k_gapped(o, max_gap) && middle_gapped(o)) a.insert(o);
    };
    string_gen(input_mer, callback);
  } else {
    auto callback = [&](std::string& o) {
      if (atmost_k_gapped(o, max_gap) && single_gapped(o)) a.insert(o);
    };
    string_gen(input_mer, callback);
  }

  std::unordered_set<std::string> b;
  auto o_callback = [&](G o) { b.insert(o.to_string()); };
  g.huddinge_neighbours(o_callback);

  std::vector<std::string> m_g;
  for (auto a_e : a) {
    if (not b.contains(a_e)) {
      m_g.push_back(a_e);
    }
  }

  std::vector<std::string> m_b;
  for (auto b_e : b) {
    if (not a.contains(b_e)) {
      m_b.push_back(b_e);
    }
  }

  if (m_g.size() || m_b.size()) {
    std::cout << "Problem when generating " << g.to_string()
              << " neighbours: " << std::endl;
    std::cout << "Missing from neighbour generation" << std::endl;
    for (auto s : m_g) {
      std::cout << s << std::endl;
    }
    std::cout << "Missing from string generation" << std::endl;
    for (auto s : m_b) {
      std::cout << s << std::endl;
    }
    return false;
  }
  return true;
}

template <class G, bool center_gap>
bool compare_all() {
  auto callback = [](G g) {
    if (not compare_single<G, center_gap>(g)) {
      exit(1);
    }
  };
  sf::gap_mer_generation<G, decltype(callback), center_gap, max_gap>(0,
                                                                     callback);
  return true;
}

int main(int argc, char const* argv[]) {
  if (argc < 2) {
    help(argv[0]);
  }

  std::string input_mer = "";
  bool compare = false;
  bool center_gap = true;
  bool use_string_gen = false;
  for (uint16_t i = 1; i < argc; ++i) {
    std::string s = argv[i];
    if (s == "-c") {
      compare = true;
    } else if (s == "-a") {
      center_gap = false;
    } else if (s == "-h") {
      help(argv[0]);
    } else if (s == "-s") {
      use_string_gen = true;
    } else {
      input_mer = s;
    }
  }

  if (input_mer.size() == 0) {
    if (not compare) {
      help(argv[0]);
    }
    bool res;
    if (center_gap) {
      res = compare_all<sf::gapmer<true, max_gap>, true>();
    } else {
      res = compare_all<sf::gapmer<false, max_gap>, false>();
    }
    if (res) {
      std::cout << "OK!" << std::endl;
    }
    exit(0);
  }

  uint16_t defined_bases = 0;
  uint16_t gap_bases = 0;
  uint16_t gap_start = input_mer.size();
  bool has_multiple_gaps = false;
  for (uint16_t i = 0; i < input_mer.size(); ++i) {
    if (input_mer[i] == 'A' || input_mer[i] == 'C' || input_mer[i] == 'G' ||
        input_mer[i] == 'T') {
      ++defined_bases;
    } else {
      if (gap_start >= input_mer.size()) {
        gap_start = i;
      }
      if (i > gap_start + gap_bases) {
        has_multiple_gaps = true;
      }
      ++gap_bases;
    }
  }
  if (has_multiple_gaps || use_string_gen) {
    if (compare) {
      if (has_multiple_gaps) {
        std::cerr << "Can't compare multigapped mers" << std::endl;
        exit(1);
      } else {
        bool res;
        if (center_gap) {
          res = comp_str<sf::gapmer<true, max_gap>, true>(
              sf::gapmer<true, max_gap>(input_mer.c_str(), defined_bases,
                                        gap_start, gap_bases),
              input_mer);
        } else {
          res = comp_str<sf::gapmer<false, max_gap>, false>(
              sf::gapmer<false, max_gap>(input_mer.c_str(), defined_bases,
                                         gap_start, gap_bases),
              input_mer);
        }
        if (res) {
          std::cout << "str_comp OK!" << std::endl;
        }
        exit(0);
      }
    }
    std::cout << "Huddinge 1 neighbours " << input_mer
              << " using string generation" << std::endl;
    string_gen(input_mer,
               [](std::string& o_mer) { std::cout << o_mer << std::endl; });
    exit(0);
  }

  gap_start = gap_start >= input_mer.size() ? 0 : gap_start;
  std::cout << defined_bases << " defined bases, " << gap_start
            << " gap start, " << gap_bases << " gap bases." << std::endl;

  if (compare) {
    bool res;
    if (defined_bases < 5 || defined_bases > 6) {
      if (center_gap) {
        res = comp_str<sf::gapmer<true, max_gap>, true>(
            sf::gapmer<true, max_gap>(input_mer.c_str(), defined_bases,
                                      gap_start, gap_bases),
            input_mer);
      } else {
        res = comp_str<sf::gapmer<false, max_gap>, false>(
            sf::gapmer<false, max_gap>(input_mer.c_str(), defined_bases,
                                       gap_start, gap_bases),
            input_mer);
      }

    } else if (center_gap) {
      res = compare_single<sf::gapmer<true, max_gap>, true>(
          sf::gapmer<true, max_gap>(input_mer.c_str(), defined_bases, gap_start,
                                    gap_bases));

    } else {
      res = compare_single<sf::gapmer<false, max_gap>, false>(
          sf::gapmer<false, max_gap>(input_mer.c_str(), defined_bases,
                                     gap_start, gap_bases));
    }
    if (res) {
      std::cout << "OK!" << std::endl;
    }
    exit(0);
  }

  std::cout << "Huddinge 1 neighbours of " << input_mer << std::endl;
  if (center_gap) {
    typedef sf::gapmer<true, max_gap> T;
    T m(input_mer.c_str(), defined_bases, gap_start, gap_bases);
    if (not m.is_valid()) {
      std::cerr << "Validation error " << m.to_string() << " <-> " << m.bits()
                << std::endl;
      exit(1);
    }
#ifdef DEBUG
    auto callback = [&](T g) {
      std::cout << m.to_string() << " -> " << g.to_string() << std::endl;
    };
#else
    auto callback = [](T g) { std::cout << g.to_string() << std::endl; };
#endif
    m.huddinge_neighbours(callback);
  } else {
    typedef sf::gapmer<false, max_gap> T;
    T m(input_mer.c_str(), defined_bases, gap_start, gap_bases);
    if (not m.is_valid()) {
      std::cerr << "Validation error " << m.to_string() << " <-> " << m.bits()
                << std::endl;
      exit(1);
    }
#ifdef DEBUG
    auto callback = [&](T g) {
      std::cout << m.to_string() << " -> " << g.to_string() << std::endl;
    };
#else
    auto callback = [](T g) { std::cout << g.to_string() << std::endl; };
#endif
    m.huddinge_neighbours(callback);
  }

  return 0;
}
