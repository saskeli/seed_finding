#include <cstdint>
#include <string>
#include <iostream>

#include "include/kmer.hpp"

int main() {
    std::string s = "AAAAAAAA";
    sf::kmer<15> km(s, 8);
    km.print();
    std::cout << km.value() << std::endl;
    s = "AAAAAAAC";
    km = sf::kmer<15>(s, 8);
    km.print();
    std::cout << km.value() << std::endl;
    s = "AAAAAAAG";
    km = sf::kmer<15>(s, 8);
    km.print();
    std::cout << km.value() << std::endl;
    s = "AAAAAAAT";
    km = sf::kmer<15>(s, 8);
    km.print();
    std::cout << km.value() << std::endl;
    s = "AAAAAACA";
    km = sf::kmer<15>(s, 8);
    km.print();
    std::cout << km.value() << std::endl;
    s = "AAAAAACC";
    km = sf::kmer<15>(s, 8);
    km.print();
    std::cout << km.value() << std::endl;
    s = "AAAAAACG";
    km = sf::kmer<15>(s, 8);
    km.print();
    std::cout << km.value() << std::endl;
    s = "GGGGGGCT";
    km = sf::kmer<15>(s, 8);
    km.print();
    std::cout << km.value() << std::endl;
    return 0;
}
