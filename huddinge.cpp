#include <cstdint>
#include <string>
#include <iostream>
#include <bitset>
#include <string>

#include "include/fm_index.hpp"

int main(int argc, char const *argv[]) {
    if (argc < 2) {
        std::cerr << "Fasta file required" << std::endl;
        exit(1);
    }
    
    sf::fm_index index(argv[1]);

    return 0;
}

