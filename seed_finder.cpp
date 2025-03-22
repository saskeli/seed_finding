#include <iostream>

#include "include/fm_index.hpp"

int main(int argc, char const *argv[]) {
    if (argc < 2) {
        std::cerr << "Input fasta is required" << std::endl;
    }
    sf::fm_index<> index(argv[1]);

    std::string s = "AGTTTTTAAGTTAGTAGTGGCCAGTCTC";
    std::cout << index.count(s) << " instances of " << s << std::endl;
    return 0;
}
