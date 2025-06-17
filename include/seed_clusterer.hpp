#include <cstdint>
#include <iostream>

namespace sf {

template <bool middle_gap_only, uint16_t max_gap, class Res_vec_T>
class seed_clusterer {
 private:
  const Res_vec_T& seeds_;
  std::vector<std::vector<size_t>> edges_;
  std::string sig_path_;
  std::string bg_path_;

 public:
  seed_clusterer(const Res_vec_T& seeds, const std::string sig_path,
                 const std::string bg_path)
      : seeds_(seeds),
        edges_(seeds.size()),
        sig_path_(sig_path),
        bg_path_(bg_path) {
#pragma omp parallel for
    for (size_t i = 0; i < seeds_.size(); ++i) {
      auto mer = seeds_[i].g;
      for (size_t j = i + 1; j < seeds.size(); ++j) {
        if (mer.template is_neighbour<true>(seeds_[j].g)) {
#pragma omp critical
          {
            edges_[i].push_back(j);
            edges_[j].push_back(i);
          }
        }
      }
    }
    for (size_t i = 0; i < seeds_.size(); ++i) {
        std::cout << seeds_[i].g.to_string() << "\t[";
        size_t lim = edges_[i].size() - 1;
        for (size_t j = 0; j <= lim; ++j) {
            std::cout << edges_[i][j] << (j == lim ? "]" : ","); 
        }
        std::cout << std::endl;
    }
  }

  size_t size() const {
    return 0;
  }

  void output_cluster() {
    return;
  }
};
}  // namespace sf
