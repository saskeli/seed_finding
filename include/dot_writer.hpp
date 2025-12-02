#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>

namespace sf {
class Dot_Writer {
 public:
  template <class G, class T>
  static void write_dot(const std::string& dot_path, const T& node_vec,
                        uint8_t max_k) {
    std::ofstream dot_file(dot_path);
    std::unordered_map<uint64_t, size_t> node_map;
    dot_file << "graph G {\n";
    for (size_t i = 0; i < node_vec.size(); ++i) {
      const auto r = node_vec[i];
      dot_file << i << " [label=\"" << r.g.to_string()
               << "\", signal_count=" << r.sig_count
               << ", background_count=" << r.bg_count << "]\n";
      node_map[uint64_t(r.g)] = i;
    }
    size_t edges = 1;
    size_t i = 0;
    auto cb = [&](G o) {
      G actual = o.is_canonical() ? o : o.reverse_complement();
      if (node_map.contains(uint64_t(actual))) {
        size_t o_i = node_map[uint64_t(actual)];
        if (i < o_i) {
          dot_file << i << " -- " << o_i << "\n";
          ++edges;
        }
      }
    };
    for (; i < node_vec.size(); ++i) {
      G n = node_vec[i].g;
      if (n.length() <= 5) {
        n.template huddinge_neighbours<true, false, false>(cb);
      } else if (n.length() >= max_k) {
        n.template huddinge_neighbours<false, false, true>(cb);
      } else {
        n.huddinge_neighbours(cb);
      }
    }
    dot_file << "}" << std::endl;
    std::cerr << node_vec.size() << " element Huddinge graph with " << edges
              << " edges written to " << dot_path << std::endl;
  }
};
} // namespace sf
