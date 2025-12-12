#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace sf {

template <typename t_type>
concept UInt64Constructible = requires(t_type val) {
	uint64_t(val);
};

// We require that the key type is convertible to uint64_t and, given t_key key,
// uint64_t(key) == 0 is a special empty value.
template <UInt64Constructible t_key, typename t_value, typename t_hash = std::hash <t_key>>
class lossy_hash_map {
 public:
  typedef t_key key_type;
  typedef t_value value_type;
  typedef t_hash hash_type;

  struct allocate_tag {};

  struct element {
    key_type key{};
    value_type value{};
  };

 private:
  constexpr static uint64_t MOD{UINT64_C(59999999)};
  constexpr static uint64_t STEP{UINT64_C(40960001)};

  std::vector <element> elements_;
  hash_type hash_;
  std::size_t size_{};

public:
  /// Construct an empty hash map.
  constexpr lossy_hash_map() = default;

  /// Allocate space for the values.
  explicit lossy_hash_map(allocate_tag):
  	elements_(MOD) {}

  // For now we assume that key_type is sufficiently small.
  void init(key_type key);
  /// Find the given key. Always returns an index.
  std::size_t find_(key_type key) const;
  /// Find the given key. Returns nullptr if not found.
  value_type *find(key_type key);
  value_type const *find(key_type key) const;
  void clear();

  double fill_rate() const { return double(size_) / MOD; }
};


template <UInt64Constructible t_key, typename t_value, typename t_hash>
void lossy_hash_map <t_key, t_value, t_hash>::init(key_type key) {
  auto &el{elements_[find_(key)]};
  if (not uint64_t(el.key)) {
    el.key = key;
    ++size_;
  }
}


template <UInt64Constructible t_key, typename t_value, typename t_hash>
std::size_t lossy_hash_map <t_key, t_value, t_hash>::find_(key_type key) const {
  auto const hh{hash_(key)};
  auto idx{hh % MOD};
  while (uint64_t(elements_[idx].key) && elements_[idx].key != key) [[unlikely]] {
    idx += STEP;
    idx -= idx >= MOD ? MOD : 0;
  }
  return idx;
}


template <UInt64Constructible t_key, typename t_value, typename t_hash>
auto lossy_hash_map <t_key, t_value, t_hash>::find(key_type key) -> value_type * {
  auto &el{elements_[find_(key)]};
  if (el.key != key) return nullptr;
    return &el.value;
}


template <UInt64Constructible t_key, typename t_value, typename t_hash>
auto lossy_hash_map <t_key, t_value, t_hash>::find(key_type key) const -> value_type const * {
  auto &el{elements_[find_(key)]};
  if (el.key != key) return nullptr;
    return &el.value;
}


template <UInt64Constructible t_key, typename t_value, typename t_hash>
void lossy_hash_map <t_key, t_value, t_hash>::clear() {
  auto const size{elements_.size()};
  elements_.clear();
  elements_.resize(size);
  size_ = 0;
}
}  // namespace sf
