#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>

#include "include/util.hpp"

int main(int argc, char const *argv[]) {
  if (argc < 2) {
    std::cerr << "input mer required" << std::endl;
    exit(1);
  }

  std::string input_mer = argv[1];

  std::cout << "Huddinge neighbours of\n " << input_mer << "\n" << std::endl;

  // For defined bases + 1
  // Add one nuc to each gap (implicit gap at start and end)

  // Start:
  for (auto n : sf::v_to_nuc) {
    std::string ne = "";
    ne.push_back(n);
    ne.append(input_mer);
    std::cout << ne << std::endl;
  }

  // Middle:
  for (size_t i = 0; i < input_mer.size(); ++i) {
    if (input_mer[i] == '.') {
      for (auto n : sf::v_to_nuc) {
        std::string ne = input_mer.substr(0, i);
        ne.push_back(n);
        ne.append(input_mer.substr(i + 1));
        std::cout << " " << ne << std::endl;
      }
    }
  }

  // End:
  for (auto n : sf::v_to_nuc) {
    std::string ne = input_mer;
    ne.push_back(n);
    std::cout << " " << ne << std::endl;
  }

  // For defined bases + 0
  // Hamming distance in defined bases

  //Hamming
  for (auto n : sf::v_to_nuc) {
    for (size_t i = 0; i < input_mer.size(); ++i) {
      if (input_mer[i] != '.' and input_mer[i] != n) {
        std::string ne = input_mer.substr(0, i);
        ne.push_back(n);
        ne.append(input_mer.substr(i + 1));
        std::cout << ne << std::endl;
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
        std::cout << ne << std::endl;
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
          std::cout << " " << nne << std::endl;
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
        std::cout << " " << ne << std::endl;
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
      std::cout << (i == 0 ? "  " : " ") << ne << std::endl;
    }
  }

  return 0;
}
