#include <SeqIO/SeqIO.hh>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

int main(int argc, char const *argv[]) {
  if (argc < 3) {
    std::cerr << "Nope" << std::endl;
    exit(1);
  }
  std::string reg = "(";
  std::ifstream file(argv[1]);
  std::string str;
  while (std::getline(file, str)) {
    if (reg.size() > 1) {
      reg.push_back('|');
    }
    reg.append(str);
  }
  reg.push_back(')');
  std::regex rex(reg);
  seq_io::Reader sr(argv[2]);
  sr.enable_reverse_complements();
  std::vector<std::pair<std::string, size_t>> hits;
  while (true) {
    uint64_t len = sr.get_next_read_to_buffer();
    if (len == 0) {
      break;
    }
    std::string read(sr.read_buf, len);
    auto e = std::sregex_iterator();
    auto m = std::sregex_iterator(read.begin(), read.end(), rex);
    for (; m != e; ++m) {
      hits.push_back({read, m->position()});
    }
  }
  size_t max_pos = 0;
  for (auto p : hits) {
    max_pos = std::max(max_pos, p.second);
  }
  std::cerr << "Max pos: " << max_pos << std::endl;
  size_t max_len = 0;
  for (auto p : hits) {
    size_t n_count = max_pos - p.second;
    max_len = std::max(n_count + p.first.size(), max_len);
  }
  std::cerr << "Max len: " << max_len << std::endl;
  size_t c = 1;
  for (auto p : hits) {
    std::cout << "> " << c++ << "\n";
    size_t n_count = max_pos - p.second;
    std::cerr << max_pos << " - " << p.second << " = " << n_count << std::endl;
    if (n_count > 200) {
      exit(1);
    }
    for (size_t i = 0; i < n_count; ++i) {
      std::cout << "N";
    }
    std::cout << p.first;
    std::cerr << max_len << " - " << p.first.size() << " - " << n_count;
    n_count = max_len - p.first.size() - n_count;
    std::cerr << " = " << n_count << std::endl;
    if (n_count > 200) {
      exit(1);
    }
    for (size_t i = 0; i < n_count; ++i) {
      std::cout << "N";
    }
    std::cout << std::endl;
  }
  return 0;
}
