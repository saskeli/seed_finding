#include <iostream>

#include "include/fm_index.hpp"

int main(int argc, char const *argv[]) {
    if (argc < 3) {
        std::cerr << "Input fasta and sequence is required" << std::endl;
        exit(1);
    }
    sf::fm_index<> index(argv[1]);

    std::string s = argv[2];
    std::cout << "Instances of " << s << std::endl;
    for (auto l : index.locate(s)) {
        std::cout << "at " << l << std::endl;
    }
    return 0;
}
