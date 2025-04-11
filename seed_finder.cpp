#include <iostream>

#include "include/fm_index.hpp"
#include "include/gapmer.hpp"

int main(int argc, char const *argv[]) {
  if (argc < 3) {
    std::cerr << "Input fasta and sequence is required" << std::endl;
    exit(1);
  }
  sf::fm_index<> index(argv[1]);

  std::string s = argv[2];
  auto gap_s = s.find(".");
  if (gap_s != std::string::npos) {
    uint64_t gap_l = gap_s + 1;
    for (; gap_l < s.size(); ++gap_l) {
      if (s[gap_l] != '.') {
        break;
      }
    }
    gap_l -= gap_s;
    std::cout << s << " with " << s.size() - gap_l << " defined bases, gap start " << gap_s
              << " and gap length " << gap_l << std::endl;
    sf::gapmer mer(s.c_str(), s.size() - gap_l, gap_s, gap_l);
    std::cout << index.count(mer) << " instances of " << mer.to_string() << std::endl;
  } else {
    sf::gapmer mer(s.c_str(), s.size());
    std::cout << index.count(mer) << " instances of " << mer.to_string() << std::endl;
  }

  return 0;
}
