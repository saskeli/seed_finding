#include <cstdint>
#include <string>
#include <iostream>

#include "include/gapmer.hpp"

int main() {
    std::string s = "AAAAAAAA";
    sf::gapmer km(s, 4, 2, 2);
    std::cout << km.to_string() << std::endl;
    return 0;
}
