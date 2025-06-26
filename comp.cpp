
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

typedef uint64_t packmer;

const constexpr std::array<char, 4> v_to_nuc = {'A', 'C', 'T', 'G'};

const constexpr static uint16_t max_gap = 10;

/// @brief Length of the gap in the gapped kmer
/// @param key
/// @return Length of the gap in bases
uint8_t get_mer_gap(const packmer key) {
  return int(((1 << 5) - 1) & (key >> (60 - 1)));
}

/// @brief  Get position of gap in bitpacked kmer
/// @param key bitpacked gapped kmer
/// @return Start index of the gap
uint8_t get_mer_gap_pos(const packmer key) {
  return int(((1 << 5) - 1) & (key >> (55 - 1)));
}

/// @brief Number of defined bases in the gapped kmer
/// @param key
/// @return Number of defined bases in the gapped kmer
uint8_t get_mer_len(const packmer key) {
  return int(((1 << 5) - 1) & (key >> (50 - 1)));
}

uint64_t bitExtracted(uint64_t number, int k, int p) {
  return int(((1 << k) - 1) & (number >> (p - 1)));
}

uint64_t base_char_to_idx(const char c) { return ((c >> 1) & 3); }

/** Convert string kmer "ACGT..TCA" to bin packed version with
 with s=ACGTTCA, gaps==2,  gap_start==4, len = 7
 */
packmer mer_to_key(const char *s, uint64_t gaps, uint64_t gap_start,
                   uint64_t len) {
  packmer key = 0;
  for (uint64_t j = 0; j < len; j++) {
    const char s_char = s[len - j - 1];
    const uint64_t c = base_char_to_idx(s_char);
    key |= c << (2 * j);
  }
  key |= gaps << (59);
  key |= gap_start << (54);
  key |= len << (49);

  return key;
}

packmer mer_to_key(const std::string &kmer) {
  uint64_t gaps = 0;
  uint64_t gap_start = UINT64_MAX;  // Initialize to an invalid value
  uint64_t len = 0;
  std::string bases;

  for (uint64_t i = 0; i < kmer.length(); i++) {
    if (kmer[i] == '.') {
      gaps++;

      if (UINT64_MAX == gap_start) {
        gap_start = i;
      }
    } else {
      bases.push_back(kmer[i]);
      len++;
    }
  }
  if (UINT64_MAX == gap_start) {
    gap_start = 0;
  }
  return mer_to_key(bases.c_str(), gaps, gap_start, len);
}

std::string key_to_mer(packmer key) {
  std::string s;
  int gap_lenght = get_mer_gap(key);
  int gap_pos = get_mer_gap_pos(key);
  int len = get_mer_len(key);
  for (int j = len; j > 0; j--) {
    if (len - j == gap_pos) {
      for (int g = 0; g < gap_lenght; g++) {
        s.append(1, '.');
      }
    }
    s.append(1, v_to_nuc[bitExtracted(key, 2, (j * 2 - 1))]);
  }

  return s;
}

std::string get_mer_bases(const packmer key) {
  std::string s;

  int len = get_mer_len(key);
  for (int j = len; j > 0; j--) {
    s.append(1, v_to_nuc[bitExtracted(key, 2, (j * 2 - 1))]);
  }

  return s;
}

std::vector<packmer> get_hamming_neighbors(packmer kmer) {
  std::vector<packmer> out;
  int len = get_mer_len(kmer) * 2;
  for (int i = 0; i < len; i += 2) {
    kmer ^= (1UL << i);
    out.push_back(kmer);

    kmer ^= (1UL << (i + 1));
    out.push_back(kmer);
    kmer ^= (1UL << i);
    out.push_back(kmer);
    kmer ^= (1UL << (i + 1));
  }

  return out;
}

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "input mer required" << std::endl;
    exit(1);
  }
  // e is a pair of (kmer, (signal_count, bg_count))
  std::string mer = argv[1];
  const packmer kmer = mer_to_key(mer);
  std::cout << "Hudding 1 neighbours for " << key_to_mer(kmer) << std::endl;

  int len = get_mer_len(kmer);
  std::string bases = get_mer_bases(kmer);
  int gap = get_mer_gap(kmer);
  int gap_pos = get_mer_gap_pos(kmer);

  // STEP 1: Same length k-mers (counts can be compared directly)
  //
  // 1.1: First go through the hamming neighbors
  std::cout << "*************ONE************" << std::endl;
  std::vector<packmer> hamming_neighbors = get_hamming_neighbors(kmer);

  for (auto pm : hamming_neighbors) {
    std::cout << key_to_mer(pm) << std::endl;
  }

  // 1.2 Shifted k-mers, same length
  std::cout << "*************ONE.TWO************" << std::endl;
  if (gap_pos != 1 && (gap_pos != len - 1)) {
    if (gap > 0) {
      for (auto a : v_to_nuc) {
        std::string temp = a + bases.substr(0, len - 1);  // extend k-mer left
        packmer mer_key = mer_to_key(temp.c_str(), gap, gap_pos + 1, len);
        std::cout << key_to_mer(mer_key) << std::endl;

        temp = bases.substr(1, len - 1) + a;  // extend k-mer right
        mer_key = mer_to_key(temp.c_str(), gap, gap_pos - 1, len);
        std::cout << key_to_mer(mer_key) << std::endl;
      }
    } else {
      for (auto a : v_to_nuc) {
        std::string temp = a + bases.substr(0, len - 1);  // extend k-mer left
        packmer mer_key = mer_to_key(temp.c_str(), 0, 0, len);
        std::cout << key_to_mer(mer_key) << std::endl;

        temp = bases.substr(1, len - 1) + a;  // extend k-mer right
        mer_key = mer_to_key(temp.c_str(), 0, 0, len);
        std::cout << key_to_mer(mer_key) << std::endl;
      }
    }
  }

  // 1.3: same number of defined bases, one shorter gap: "shorter gap"
  std::cout << "*************ONE.THREE************" << std::endl;
  if (gap > 0) {
    if (gap == 1) {
      for (auto a : v_to_nuc) {
        std::string temp =
            bases.substr(1, gap_pos - 1) + a +
            bases.substr(
                gap_pos,
                len - gap_pos);  // gap == 1, fill the gap shorten from start

        packmer mer_key = mer_to_key(temp.c_str(), 0, 0, len);
        std::cout << key_to_mer(mer_key) << std::endl;
      }
      for (auto a : v_to_nuc) {
        std::string temp =
            bases.substr(0, gap_pos) + a +
            bases.substr(
                gap_pos,
                len - gap_pos - 1);  // gap == 1, fill the gap shorten from end
        packmer mer_key = mer_to_key(temp.c_str(), 0, 0, len);
        std::cout << key_to_mer(mer_key) << std::endl;
      }
    } else {
      // if there is a longer gap, it must be reduced from both directions
      for (auto a : v_to_nuc) {
        // 1.shortening kmer from start, shortening gap from start
        std::string temp = bases.substr(1, gap_pos - 1) + a +
                           bases.substr(gap_pos, len - gap_pos);
        packmer mer_key = mer_to_key(temp.c_str(), gap - 1, gap_pos, len);
        std::cout << key_to_mer(mer_key) << std::endl;

        // 2. shortening kmer from end, shortening gap from start
        temp = bases.substr(0, gap_pos) + a +
               bases.substr(gap_pos, len - gap_pos - 1);
        mer_key = mer_to_key(temp.c_str(), gap - 1, gap_pos + 1, len);
        std::cout << key_to_mer(mer_key) << std::endl;

        // 3. shortening the kmer from start, shortening gap from end
        temp = bases.substr(1, gap_pos - 1) + a +
               bases.substr(gap_pos, len - gap_pos);
        mer_key = mer_to_key(temp.c_str(), gap - 1, gap_pos - 1, len);
        std::cout << key_to_mer(mer_key) << std::endl;

        // 4. shortening the kmer from end, shortening gap from end
        temp = bases.substr(0, gap_pos) + a +
               bases.substr(gap_pos, len - gap_pos - 1);
        mer_key = mer_to_key(temp.c_str(), gap - 1, gap_pos, len);
        std::cout << key_to_mer(mer_key) << std::endl;
      }
    }
  }

  std::cout << "*************ONE.ZERO?************" << std::endl;
  // 1.0: remove first/last letter of non-gapped k-mer, add arbitrary length
  // gap
  // This has no effect if only middle gapped -kmers are counted same
  // number of defined bases, counts can be directly compared
  if (gap == 0) {
    for (int g = 1; g <= max_gap; g++)  // add gap until max_gap
    {
      for (auto a : v_to_nuc) {
        // shortening the kmer from beginning
        std::string temp = bases.substr(1, len - 1) + a;
        packmer mer_key = mer_to_key(temp.c_str(), g, len - 1, len);
        std::cout << key_to_mer(mer_key) << std::endl;

        // removing last letter, adding gap to beginning
        temp = a + bases.substr(0, len - 1);
        mer_key = mer_to_key(temp.c_str(), g, 1, len);
        std::cout << key_to_mer(mer_key) << std::endl;
      }
    }
  }

  /***************************************************************************************************
   *  STEP 2: One shorter Huddinge-neighbors (counts cannot be directly
   *compared, need to test for extension)
   *
   * 2.1: longer gaps
   * testing if the current k-mer can be extended from shorter ones that have
   *a longer gap, for example: seed=ACnT, can it be extended from AnnT (so
   *ACnT must have higher fraction of counts counts than thresholds[k-1]['C'])
   *************************************************************************************************/
  std::cout << "*************TWO.ONE************" << std::endl;
  if (gap > 0 && (gap < ((len + gap) - 2))) {
    if (gap_pos > 1) {
      std::cout << bases.substr(0, gap_pos - 1);
      for (uint16_t g_i = 0; g_i <= gap; ++g_i) {
        std::cout << ".";
      }
      std::cout << bases.substr(gap_pos) << std::endl;
      for (auto a : v_to_nuc) {
        // shortening the kmer from beginning
        std::string temp = bases.substr(0, gap_pos - 1) + a +
                           bases.substr(gap_pos, len - gap_pos);
        packmer mer_key = mer_to_key(temp.c_str(), gap, gap_pos, len);
        std::cout << "** -> " << key_to_mer(mer_key) << std::endl;
      }
    }
    if (gap_pos < len - 1) {
      std::cout << bases.substr(0, gap_pos);
      for (uint16_t g_i = 0; g_i <= gap; ++g_i) {
        std::cout << ".";
      }
      std::cout << bases.substr(gap_pos + 1) << std::endl;
      for (auto a : v_to_nuc) {
        std::string temp = bases.substr(0, gap_pos) + a +
                           bases.substr(gap_pos + 1, len - gap_pos - 1);
        packmer mer_key = mer_to_key(temp.c_str(), gap, gap_pos, len);
        std::cout << "** -> " << key_to_mer(mer_key) << std::endl;
      }
    }
  }

  if (gap == 0) {
    // shorter with longer gap, when kmer does not have a gap
    for (int g_ind = 1; g_ind < len - 1; g_ind++) {
      for (auto a : v_to_nuc) {
        // shortening the kmer from beginning
        std::string temp = bases.substr(0, g_ind) + a +
                           bases.substr(g_ind + 1, len - g_ind - 1);
        packmer mer_key = mer_to_key(temp.c_str(), 0, 0, len);
        std::cout << key_to_mer(mer_key) << std::endl;
      }
    }
  }

  // 2.2: shorter Huddinge neighbors
  std::cout << "*************TWO.TWO************" << std::endl;
  if (gap_pos != 1) {
    // extension to the first position of kmer
    std::cout << key_to_mer(mer_to_key(bases.substr(1).c_str(), gap,
                                       gap_pos - 1, len - 1))
              << std::endl;
    for (auto a : v_to_nuc) {
      // shortening the kmer from beginning
      std::string temp = a + bases.substr(1, len - 1);
      packmer mer_key = mer_to_key(temp.c_str(), gap, gap_pos, len);
      std::cout << "** -> " << key_to_mer(mer_key) << std::endl;
    }
  }
  if (gap_pos != len - 1) {
    // extension to the last position of kmer
    std::cout << key_to_mer(mer_to_key(bases.substr(0, len - 1).c_str(), gap,
                                       gap_pos, len - 1))
              << std::endl;
    for (auto a : v_to_nuc) {
      // shortening the kmer from beginning
      std::string temp = bases.substr(0, len - 1) + a;
      packmer mer_key = mer_to_key(temp.c_str(), gap, gap_pos, len);
      std::cout << "** -> " << key_to_mer(mer_key) << std::endl;
    }
  }

  /*********************************************************************************************************
   * if we got this far, it means that kmer is local max within its own length
   *and k-1 Huddinge-neighbors! next testing if kmer extends to any k+1-mers,
   *unless k == kmax, as for these greedy search is used
   ********************************************************************************************************
   */
  std::cout << "*************L.ONE************" << std::endl;
  // 1: same number of defined bases but one longer gap
  if (gap == 0) {
    // if there is no gap in kmer, only gaps of length 1 are considered
    for (int g_ind = 1; g_ind < (len - 1); g_ind++) {
      // first extending kmer to left
      for (auto a : v_to_nuc) {
        // shortening the kmer from end, shortening gap from end
        std::string temp = a + bases.substr(0, g_ind) +
                           bases.substr(g_ind + 1, len - g_ind - 1);
        packmer mer_key = mer_to_key(temp.c_str(), 1, g_ind + 1, len);
        std::cout << key_to_mer(mer_key) << std::endl;

        temp = bases.substr(0, g_ind) +
               bases.substr(g_ind + 1, len - g_ind - 1) + a;
        mer_key = mer_to_key(temp.c_str(), 1, g_ind, len);
        std::cout << key_to_mer(mer_key) << std::endl;
      }
    }
  }

  if (gap > 0) {
    // 1.2:if there is a gap, we must test extending the gap to both
    std::cout << "*************L.ONE.TWO************" << std::endl;
    // directions
    int g_ind = gap_pos;
    if (g_ind > 0) {
      for (auto a : v_to_nuc) {
        // extending kmer towards beginning and gap towards beginning
        std::string temp = a + bases.substr(0, gap_pos - 1) +
                           bases.substr(gap_pos, len - gap_pos);
        packmer mer_key = mer_to_key(temp.c_str(), gap + 1, gap_pos, len);
        std::cout << key_to_mer(mer_key) << std::endl;
        // extending kmer towards end and gap towards beginning
        if (gap_pos > 1) {
          temp = bases.substr(0, gap_pos - 1) +
                 bases.substr(g_ind, len - gap_pos) + a;
          mer_key = mer_to_key(temp.c_str(), gap + 1, gap_pos - 1, len);
          std::cout << key_to_mer(mer_key) << std::endl;
        }
      }
    }

    for (auto a : v_to_nuc) {
      // 1.3:extending kmer towards end and gap towards end
      std::cout << "*************L.ONE.THREE************" << std::endl;
      std::string temp = bases.substr(0, gap_pos) +
                         bases.substr(gap_pos + 1, len - gap_pos - 1) + a;
      packmer mer_key = mer_to_key(temp.c_str(), gap + 1, gap_pos, len);
      std::cout << key_to_mer(mer_key) << std::endl;

      // 1.4:extending kmer towards beginning and gap tpwards end
      std::cout << "*************L.ONE.FOUR************" << std::endl;
      if (gap_pos < (len - 1)) {
        temp = a + bases.substr(0, gap_pos) +
               bases.substr(gap_pos + 1, len - gap_pos - 1);
        mer_key = mer_to_key(temp.c_str(), gap + 1, gap_pos + 1, len);
        std::cout << key_to_mer(mer_key) << std::endl;
      }
    }
  }

  // 2: Trying to extend seed towards either end, one more defined bases, so
  std::cout << "*************L.TWO************" << std::endl;
  // comparison using threshold
  for (auto a : v_to_nuc) {
    // extending to left
    std::string temp = a + bases;
    packmer mer_key;
    if (gap > 0) {
      mer_key = mer_to_key(temp.c_str(), gap, gap_pos + 1, len + 1);
    } else {
      mer_key = mer_to_key(temp.c_str(), gap, gap_pos, len + 1);
    }
    std::cout << key_to_mer(mer_key) << std::endl;
  }

  for (auto a : v_to_nuc) {
    // extending to right
    std::string temp = bases + a;
    packmer mer_key = mer_to_key(temp.c_str(), gap, gap_pos, len + 1);
    std::cout << key_to_mer(mer_key) << std::endl;
  }

  // 3: Extending gapped k-mer at the position of the gap, one more defined
  std::cout << "*************L.THREE************" << std::endl;
  // bases, comparison using threshold
  if (gap > 0) {
    int chances = (gap == 1) ? 1 : 2;  // if gap == 1 only fill the gap once.
                                       // If gap>1 fill start and end of gap

    for (int i = 0; i < chances; i++) {
      for (auto a : v_to_nuc) {
        std::string temp =
            bases.substr(0, gap_pos) + a + bases.substr(gap_pos, len - gap_pos);
        int temp_gap_pos = gap_pos + i;
        int temp_gap = gap - 1;

        if (temp_gap == 0) {
          temp_gap_pos = 0;
        }

        packmer mer_key =
            mer_to_key(temp.c_str(), temp_gap, temp_gap_pos, len + 1);
        std::cout << key_to_mer(mer_key) << std::endl;
      }
    }
  }
  std::cout << "***************L.FOUR*************" << std::endl;
  // 4: Extending non-gapped kmer by adding arbitrary length gap followed by
  // any nucleotide. One more defined bases, comparison using threshold. (this
  // should have no meaning if only middle gapped k-mers are counted).
  if (gap == 0) {
    for (int g = 1; g <= max_gap; g++) {
      for (auto a : v_to_nuc) {
        // extend right
        std::string temp = bases + a;
        packmer mer_key = mer_to_key(temp.c_str(), g, len, len + 1);
        std::cout << key_to_mer(mer_key) << std::endl;
      }

      for (auto a : v_to_nuc) {
        // extend left
        std::string temp = a + bases;
        packmer mer_key = mer_to_key(temp.c_str(), g, 1, len + 1);
        std::cout << key_to_mer(mer_key) << std::endl;
      }
    }
  }
}
