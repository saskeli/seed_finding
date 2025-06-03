#pragma once

#include <cstdint>
#include <gtest/gtest.h>
#include <string>
#include <vector>
#include "../include/gapmer.hpp"
#include "../include/string_buffer.hpp"
#include "../include/util.hpp"

namespace sf {

std::vector<uint64_t> vec_from_string(const std::string& s) {
  if (s.empty()) return {};

  std::vector<uint64_t> vec(8U * (s.length() - 1) + 1); // Causes the string to be not NUL-terminated.
  char* ptr = reinterpret_cast<char*>(vec.data());
  for (uint16_t i = 0; i < s.size(); ++i) {
    ptr[i] = s[i];
  }
  return vec;
}

TEST(gapmer, DefaultConstructor) {
  gapmer gg;
  ASSERT_TRUE(gg.is_valid());
  ASSERT_EQ("", gg.to_string());
}

TEST(gapmer, StrConstructor1) {
  std::string s = "TTTTTTTT";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 8);
  ASSERT_TRUE(g.is_valid());
  auto sm = g.to_string();
  ASSERT_EQ(sm, s);
}

TEST(gapmer, StrConstructor2) {
  std::string s = "TTT..TTT";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 6, 3, 2);
  ASSERT_TRUE(g.is_valid());
  auto sm = g.to_string();
  ASSERT_EQ(sm, s);
}

TEST(gapmer, StrConstructor3) {
  std::string s = "CATTATTAC";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 9, 0, 0);
  ASSERT_TRUE(g.is_valid());
  auto sm = g.to_string();
  ASSERT_EQ(sm, s);
}

TEST(gapmer, StrConstructor4) {
  std::string s = "ATATATATAT...CATTATTAC";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 19, 10, 3);
  ASSERT_TRUE(g.is_valid());
  auto sm = g.to_string();
  ASSERT_EQ(sm, s);
}

TEST(gapmer, StrConstructor5) {
  std::string s = "AATGATT..........TCTGTGG";
  auto vec = vec_from_string(s);
  gapmer<true, 10> g(vec.data(), 14, 7, 10);
  ASSERT_TRUE(g.is_valid()) << g.to_string() << " " << g.bits();
  auto sm = g.to_string();
  ASSERT_EQ(sm, s);
}

TEST(gapmer, StrConstructor6) {
  sf::string_buffer <uint64_t> str;
  str = "CCCCC";
  const gapmer<> gg(str.data(), str.size(), 0, 0);
  ASSERT_TRUE(gg.is_valid()) << gg.to_string() << " " << gg.bits();
  const auto sm(gg.to_string());
  ASSERT_EQ(sm, str);
}

TEST(gapmer, StrConstructor7) {
  sf::string_buffer <uint64_t> str;
  str = "C.CCC";
  const gapmer<> gg(str.data(), str.size() - 1, 1, 1);
  ASSERT_TRUE(gg.is_valid()) << gg.to_string() << " " << gg.bits();
  const auto sm(gg.to_string());
  ASSERT_EQ(sm, str);
}

TEST(gapmer, StrConstructor8) {
  sf::string_buffer <uint64_t> str;
  str = "AAAA.....A";
  const gapmer<> gg(str.data(), 5, 4, 5);
  ASSERT_TRUE(gg.is_valid()) << gg.to_string() << " " << gg.bits();
  const auto sm(gg.to_string());
  ASSERT_EQ(sm, str);
}

TEST(gapmer, NextComp1) {
  std::string s = "CATTATAC";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 5);
  g = g.next(s[5]);
  vec = vec_from_string(s.substr(1));
  gapmer o(vec.data(), 5);
  ASSERT_TRUE(g.is_valid());
  ASSERT_TRUE(g == o) << g.to_string() << " <-> " << o.to_string();
}

TEST(gapmer, NextComp2) {
  std::string s = "CATTATAC";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 5, 2, 1);
  g = g.next(s[2], s[6]);
  vec = vec_from_string(s.substr(1));
  gapmer o(vec.data(), 5, 2, 1);
  ASSERT_TRUE(g.is_valid());
  ASSERT_TRUE(g == o) << g.to_string() << " <-> " << o.to_string();
}

TEST(gapmer, Neighbour1) {
  std::string s = "TTTTTTTT";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 8);
  ASSERT_FALSE(g.is_neighbour(g));
}

TEST(gapmer, Neighbour2) {
  std::string s = "TTTTTTTT";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 8);
  gapmer h(vec.data(), 7);
  ASSERT_TRUE(g.is_neighbour(h));
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour3) {
  std::string s = "TTTTTTTT";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 8);
  gapmer h(vec.data(), 6);
  ASSERT_FALSE(g.is_neighbour(h));
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour4) {
  std::string s = "TTTTTTTT";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 8);
  std::string t = "AAAAAAAA";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 6);
  ASSERT_FALSE(g.is_neighbour(h));
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour5) {
  std::string s = "TTTTTTTT";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 8);
  std::string t = "TTTTTTAA";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 6);
  ASSERT_FALSE(g.is_neighbour(h));
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour6) {
  std::string s = "CATTATTAC";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 9);
  std::string t = "ATTATTAC";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 8);
  ASSERT_TRUE(g.is_neighbour(h));
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour7) {
  std::string s = "CATTATTAC";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 9);
  std::string t = "CATTATTA";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 8);
  ASSERT_TRUE(g.is_neighbour(h));
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour8) {
  std::string s = "CATTATTAC";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 9);
  std::string t = "CCATTATTA";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 9);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour9) {
  std::string s = "CATTATTAC";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 9);
  gapmer h(vec.data(), 8, 4, 1);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour10) {
  std::string s = "AA.AAT";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 5, 2, 1);
  std::string t = "AAAA";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 4);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour11) {
  std::string s = "AC.ATG";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 5, 2, 1);
  std::string t = "C.ATGG";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 5, 1, 1);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour12) {
  std::string s = "ACATG";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 5);
  std::string t = "ACTGC.T";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 6, 5, 1);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour13) {
  std::string s = "ACATG";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 5);
  std::string t = "C...ACATG";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 6, 1, 3);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour14) {
  std::string s = "ACATG";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 5);
  std::string t = "T.CATG";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 5, 1, 1);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour15) {
  std::string s = "ACATG";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 5);
  std::string t = "G.....GTGTA";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 6, 1, 5);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour16) {
  std::string s = "ACATG";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 5);
  std::string t = "C.....GCCA";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 6, 1, 5);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour17) {
  std::string s = "ACATG";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 5);
  std::string t = "T....GATG";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 5, 1, 4);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour18) {
  std::string s = "AA.....AAA";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 5, 2, 5);
  std::string t = "AA......AAG";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 5, 2, 6);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour19) {
  std::string s = "AAAAAA";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 6);
  std::string t = "AGA.AAAA";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 7, 3, 1);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, RCNeighbour1) {
  std::string s = "CATTA...TTACC";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 10, 3, 5);
  gapmer o = g.reverse_complement();
  auto callback = [&](const decltype(o)& n) {
    ASSERT_TRUE(g.is_neighbour<true>(n));
  };
  o.huddinge_neighbours(callback);
}

TEST(gapmer, RCNeighbour2) {
  std::string s = "CAAT.TTAC";
  auto vec = vec_from_string(s);
  gapmer g(vec.data(), 8, 1, 4);
  std::string t = "AAAA.AAAA";
  vec = vec_from_string(t);
  gapmer h(vec.data(), 8, 1, 4);
  ASSERT_FALSE(g.is_neighbour<true>(h))
      << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour<true>(g));
}

TEST(gapmer, Huddinge1) {
  std::string s = "ACGTGT";
  auto vec = vec_from_string(s);
  typedef gapmer<false, 5> G;
  G g(vec.data(), 6);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge2) {
  std::string s = "A.GTGT";
  auto vec = vec_from_string(s);
  typedef gapmer<false, 5> G;
  G g(vec.data(), 5, 1, 1);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge3) {
  std::string s = "ACGT.T";
  auto vec = vec_from_string(s);
  typedef gapmer<false, 5> G;
  G g(vec.data(), 5, 4, 1);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge4) {
  std::string s = "AC.TCT";
  auto vec = vec_from_string(s);
  typedef gapmer<false, 5> G;
  G g(vec.data(), 5, 4, 1);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge5) {
  std::string s = "ACGTGT";
  auto vec = vec_from_string(s);
  typedef gapmer<true, 5> G;
  G g(vec.data(), 6);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge6) {
  std::string s = "ACG.GT";
  auto vec = vec_from_string(s);
  typedef gapmer<true, 5> G;
  G g(vec.data(), 5, 3, 1);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge7) {
  std::string s = "ACG.TGT";
  auto vec = vec_from_string(s);
  typedef gapmer<true, 5> G;
  G g(vec.data(), 5, 3, 1);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge8) {
  std::string s = "AC..TCT";
  auto vec = vec_from_string(s);
  typedef gapmer<true, 5> G;
  G g(vec.data(), 5, 2, 2);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge9) {
  std::string s = "GAC..TCT";
  auto vec = vec_from_string(s);
  typedef gapmer<true, 5> G;
  G g(vec.data(), 6, 3, 2);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge10) {
  std::string s = "GAC..TC";
  auto vec = vec_from_string(s);
  typedef gapmer<true, 5> G;
  G g(vec.data(), 5, 3, 2);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge11) {
  std::string s = "GACTC";
  auto vec = vec_from_string(s);
  typedef gapmer<true, 5> G;
  G g(vec.data(), 5, 0, 0);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge12) {
  std::string s = "AC.TCT";
  auto vec = vec_from_string(s);
  typedef gapmer<true, 5> G;
  G g(vec.data(), 5, 2, 1);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge13) {
  std::string s = "A...CTCT";
  auto vec = vec_from_string(s);
  typedef gapmer<false, 5> G;
  G g(vec.data(), 5, 1, 3);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge14) {
  std::string s = "AG.TCT";
  auto vec = vec_from_string(s);
  typedef gapmer<false, 5> G;
  G g(vec.data(), 5, 2, 1);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge15) {
  std::string s = "AGTC...T";
  auto vec = vec_from_string(s);
  typedef gapmer<false, 5> G;
  G g(vec.data(), 5, 4, 3);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge16) {
  std::string s = "AG...TCT";
  auto vec = vec_from_string(s);
  typedef gapmer<false, 5> G;
  G g(vec.data(), 5, 2, 3);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, IsCanonical1) {
  std::string s = "AAAA...TTT";
  auto vec = vec_from_string(s);
  // rc = AAA...TTTT so technically smaller but ignore gap
  gapmer<true, 5> g(vec.data(), 7, 4, 3);
  ASSERT_TRUE(g.is_canonical());
}

TEST(gapmer, IsCanonical2) {
  std::string s = "ACGCCGT";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 7);
  ASSERT_TRUE(g.is_canonical());
}

TEST(gapmer, IsCanonical3) {
  std::string s = "ACGGCGT";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 7);
  ASSERT_FALSE(g.is_canonical());
}

TEST(gapmer, IsCanonical4) {
  std::string s = "AAATTT";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 6);
  ASSERT_TRUE(g.is_canonical());
}

TEST(gapmer, IsCanonical5) {
  std::string s = "GGGTTT";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 6);
  ASSERT_FALSE(g.is_canonical());
}

TEST(gapmer, RevComp1) {
  std::string s = "AAAATTTT";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 8);
  g = g.reverse_complement();
  ASSERT_EQ(g.to_string(), s);
}

TEST(gapmer, RevComp2) {
  std::string s = "AAA.TTTT";
  std::string rc = "AAAA.TTT";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 7, 3, 1);
  g = g.reverse_complement();
  ASSERT_EQ(g.to_string(), rc);
}

TEST(gapmer, RevComp3) {
  std::string s = "ACGTA";
  std::string rc = "TACGT";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 5);
  g = g.reverse_complement();
  ASSERT_EQ(g.to_string(), rc);
}

TEST(gapmer, RevComp4) {
  std::string s = "ACGT";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 4);
  g = g.reverse_complement();
  ASSERT_EQ(g.to_string(), s);
}

TEST(gapmer, Align1) {
  std::string s = "AAAAA";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 5);
  bool res = g.aligns_to(g);
  ASSERT_TRUE(res);
  gapmer<true, 5> o = g.reverse_complement();
  res = g.aligns_to(o);
  ASSERT_TRUE(res);
}

TEST(gapmer, Align2) {
  std::string s = "AAAAA";
  auto vec = vec_from_string(s);
  std::string os = "GGGGAAAAAT";
  gapmer<true, 5> g(vec.data(), 5);
  vec = vec_from_string(os);
  gapmer<true, 5> o(vec.data(), 10);
  bool res = g.aligns_to(o);
  ASSERT_TRUE(res);
}

TEST(gapmer, Align3) {
  std::string s = "AAAAA";
  std::string os = "GGGG...AAAAA";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 5);
  vec = vec_from_string(os);
  gapmer<true, 5> o(vec.data(), 9, 4, 3);
  bool res = g.aligns_to(o);
  ASSERT_TRUE(res);
}

TEST(gapmer, Align4) {
  std::string s = "AAA...AA";
  std::string os = "GAAA...AATTT";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 5, 3, 3);
  vec = vec_from_string(os);
  gapmer<true, 5> o(vec.data(), 9, 4, 3);
  bool res = g.aligns_to(o);
  ASSERT_TRUE(res);
}

TEST(gapmer, Align5) {
  std::string s = "GTGAATGATTAA";
  auto vec = vec_from_string(s);
  gapmer<true, 5> g(vec.data(), 9, 4, 3);
  gapmer<true, 5> o(vec.data(), 12);
  bool res = g.template aligns_to(o);
  ASSERT_TRUE(res);
}

TEST(gapmer, HuddingeDistance1) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const sb("ACGT");
	gapmer_type const lhs_(sb.data(), sb.size(), 0, 0);
	gapmer_type const rhs_(sb.data(), sb.size(), 0, 0);
	ASSERT_EQ(0, lhs_.huddinge_distance(rhs_));
	ASSERT_EQ(0, rhs_.huddinge_distance(lhs_));
}

TEST(gapmer, HuddingeDistance2) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const lhs("ACGT");
	string_buffer_type const rhs("A.GT");
	gapmer_type const lhs_(lhs.data(), lhs.size(), 0, 0);
	gapmer_type const rhs_(rhs.data(), rhs.size() - 1, 1, 1);
	ASSERT_EQ(1, lhs_.huddinge_distance(rhs_));
	ASSERT_EQ(1, rhs_.huddinge_distance(lhs_));
}

TEST(gapmer, HuddingeDistance3) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const lhs("ACGT");
	string_buffer_type const rhs("AC.T");
	gapmer_type const lhs_(lhs.data(), lhs.size(), 0, 0);
	gapmer_type const rhs_(rhs.data(), rhs.size() - 1, 2, 1);
	ASSERT_EQ(1, lhs_.huddinge_distance(rhs_));
	ASSERT_EQ(1, rhs_.huddinge_distance(lhs_));
}

TEST(gapmer, HuddingeDistance4) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const lhs("AAAAA");
	string_buffer_type const rhs("AAAAAC");
	gapmer_type const lhs_(lhs.data(), lhs.size(), 0, 0);
	gapmer_type const rhs_(rhs.data(), rhs.size(), 0, 0);
	ASSERT_EQ(1, lhs_.huddinge_distance(rhs_));
	ASSERT_EQ(1, rhs_.huddinge_distance(lhs_));
}

TEST(gapmer, HuddingeDistance5) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const lhs("CCCCC");
	string_buffer_type const rhs("CCCCCG");
	gapmer_type const lhs_(lhs.data(), lhs.size(), 0, 0);
	gapmer_type const rhs_(rhs.data(), rhs.size(), 0, 0);
	ASSERT_EQ(1, lhs_.huddinge_distance(rhs_));
	ASSERT_EQ(1, rhs_.huddinge_distance(lhs_));
}

TEST(gapmer, HuddingeDistance6) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const lhs("AAAAA");
	string_buffer_type const rhs("AAAAA....A");
	gapmer_type const lhs_(lhs.data(), lhs.size(), 0, 0);
	gapmer_type const rhs_(rhs.data(), rhs.size() - 4, 5, 4);
	ASSERT_EQ(1, lhs_.huddinge_distance(rhs_));
	ASSERT_EQ(1, rhs_.huddinge_distance(lhs_));
}

TEST(gapmer, HuddingeDistance7) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const lhs("GGGGG");
	string_buffer_type const rhs("GGGGG....G");
	gapmer_type const lhs_(lhs.data(), lhs.size(), 0, 0);
	gapmer_type const rhs_(rhs.data(), rhs.size() - 4, 5, 4);
	ASSERT_EQ(1, lhs_.huddinge_distance(rhs_));
	ASSERT_EQ(1, rhs_.huddinge_distance(lhs_));
}

TEST(gapmer, HuddingeDistance8) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const lhs("AAAAAAAAAAAAAAAAAAAAAAA");
	string_buffer_type const rhs("AAAAAAAAAAAAAAAAAAAAAAAC");
	gapmer_type const lhs_(lhs.data(), lhs.size(), 0, 0);
	gapmer_type const rhs_(rhs.data(), rhs.size(), 0, 0);
	ASSERT_EQ(1, lhs_.huddinge_distance(rhs_));
	ASSERT_EQ(1, rhs_.huddinge_distance(lhs_));
}

TEST(gapmer, HuddingeDistance9) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const lhs("TGTTTTGTGCATAATCG...TCGC");
	string_buffer_type const rhs("CCATGTTTTGTGCATAAT...CGTCGC");
	gapmer_type const lhs_(lhs.data(), lhs.size() - 3, 17, 3);
	gapmer_type const rhs_(rhs.data(), rhs.size() - 3, 18, 3);
	ASSERT_EQ(5, lhs_.huddinge_distance(rhs_));
	ASSERT_EQ(5, rhs_.huddinge_distance(lhs_));
}

TEST(gapmer, HuddingeNeighbourhood1) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const ss("AAAAAAAAAAAAAAAAAAAAAAA");
	gapmer_type gg(ss.data(), 23, 0, 0);
	gg.all_gap_neighbours<false, false, false>([&](gapmer_type const gg_){
		EXPECT_EQ(1, gg.huddinge_distance(gg_)) << gg.to_string() << " (length: " << gg.length() << "), " << gg_.to_string() << " (length: " << gg_.length() << ")\n" << gg.bits() << '\n' << gg_.bits();
	});
}

TEST(gapmer, HuddingeNeighbourhood2) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const ss("CCCCCCCCCCCCCCCCCCCCCCC");
	gapmer_type gg(ss.data(), 23, 0, 0);
	gg.all_gap_neighbours<false, false, false>([&](gapmer_type const gg_){
		EXPECT_EQ(1, gg.huddinge_distance(gg_)) << gg.to_string() << " (length: " << gg.length() << "), " << gg_.to_string() << " (length: " << gg_.length() << ")\n" << gg.bits() << '\n' << gg_.bits();
	});
}

TEST(gapmer, HuddingeNeighbourhood3) {
	typedef sf::string_buffer<uint64_t> string_buffer_type;
	typedef sf::gapmer<> gapmer_type;
	string_buffer_type const ss("TGCTCAT..TAGGAAGCCGCCTTGA");
	gapmer_type gg(ss.data(), 23, 7, 2);
	gg.all_gap_neighbours<false, false, false>([&](gapmer_type const gg_){
		EXPECT_EQ(1, gg.huddinge_distance(gg_)) << gg.to_string() << " (length: " << gg.length() << "), " << gg_.to_string() << " (length: " << gg_.length() << ")\n" << gg.bits() << '\n' << gg_.bits();
	});
}

}  // namespace sf
