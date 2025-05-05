#pragma once

#include <string>

namespace sf {
TEST(gapmer, StrConstructor1) {
  std::string s = "TTTTTTTT";
  gapmer g(s.c_str(), 8);
  ASSERT_TRUE(g.is_valid());
  auto sm = g.to_string();
  ASSERT_EQ(sm, s);
}

TEST(gapmer, StrConstructor2) {
  std::string s = "TTT..TTT";
  gapmer g(s.c_str(), 6, 3, 2);
  ASSERT_TRUE(g.is_valid());
  auto sm = g.to_string();
  ASSERT_EQ(sm, s);
}

TEST(gapmer, StrConstructor3) {
  std::string s = "CATTATTAC";
  gapmer g(s.c_str(), 9, 0, 0);
  ASSERT_TRUE(g.is_valid());
  auto sm = g.to_string();
  ASSERT_EQ(sm, s);
}

TEST(gapmer, StrConstructor4) {
  std::string s = "ATATATATAT...CATTATTAC";
  gapmer g(s.c_str(), 19, 10, 3);
  ASSERT_TRUE(g.is_valid());
  auto sm = g.to_string();
  ASSERT_EQ(sm, s);
}

TEST(gapmer, StrConstructor5) {
  std::string s = "AATGATT..........TCTGTGG";
  gapmer<true, 10> g(s.c_str(), 14, 7, 10);
  ASSERT_TRUE(g.is_valid()) << g.to_string() << " " << g.bits();
  auto sm = g.to_string();
  ASSERT_EQ(sm, s);
}

TEST(gapmer, NextComp1) {
  std::string s = "CATTATAC";
  gapmer g(s.c_str(), 5);
  g = g.next(s[5]);
  gapmer o(s.c_str() + 1, 5);
  ASSERT_TRUE(g.is_valid());
  ASSERT_TRUE(g == o) << g.to_string() << " <-> " << o.to_string();
}

TEST(gapmer, NextComp2) {
  std::string s = "CATTATAC";
  gapmer g(s.c_str(), 5, 2, 1);
  g = g.next(s[2], s[6]);
  gapmer o(s.c_str() + 1, 5, 2, 1);
  ASSERT_TRUE(g.is_valid());
  ASSERT_TRUE(g == o) << g.to_string() << " <-> " << o.to_string();
}

TEST(gapmer, Neighbour1) {
  std::string s = "TTTTTTTT";
  gapmer g(s.c_str(), 8);
  ASSERT_FALSE(g.is_neighbour(g));
}

TEST(gapmer, Neighbour2) {
  std::string s = "TTTTTTTT";
  gapmer g(s.c_str(), 8);
  gapmer h(s.c_str(), 7);
  ASSERT_TRUE(g.is_neighbour(h));
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour3) {
  std::string s = "TTTTTTTT";
  gapmer g(s.c_str(), 8);
  gapmer h(s.c_str(), 6);
  ASSERT_FALSE(g.is_neighbour(h));
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour4) {
  std::string s = "TTTTTTTT";
  std::string t = "AAAAAAAA";
  gapmer g(s.c_str(), 8);
  gapmer h(t.c_str(), 6);
  ASSERT_FALSE(g.is_neighbour(h));
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour5) {
  std::string s = "TTTTTTTT";
  std::string t = "TTTTTTAA";
  gapmer g(s.c_str(), 8);
  gapmer h(t.c_str(), 6);
  ASSERT_FALSE(g.is_neighbour(h));
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour6) {
  std::string s = "CATTATTAC";
  std::string t = "ATTATTAC";
  gapmer g(s.c_str(), 9);
  gapmer h(t.c_str(), 8);
  ASSERT_TRUE(g.is_neighbour(h));
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour7) {
  std::string s = "CATTATTAC";
  std::string t = "CATTATTA";
  gapmer g(s.c_str(), 9);
  gapmer h(t.c_str(), 8);
  ASSERT_TRUE(g.is_neighbour(h));
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour8) {
  std::string s = "CATTATTAC";
  std::string t = "CCATTATTA";
  gapmer g(s.c_str(), 9);
  gapmer h(t.c_str(), 9);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour9) {
  std::string s = "CATTATTAC";
  gapmer g(s.c_str(), 9);
  gapmer h(s.c_str(), 8, 4, 1);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour10) {
  std::string s = "AA.AAT";
  std::string t = "AAAA";
  gapmer g(s.c_str(), 5, 2, 1);
  gapmer h(t.c_str(), 4);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour11) {
  std::string s = "AC.ATG";
  std::string t = "C.ATGG";
  gapmer g(s.c_str(), 5, 2, 1);
  gapmer h(t.c_str(), 5, 1, 1);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour12) {
  std::string s = "ACATG";
  std::string t = "ACTGC.T";
  gapmer g(s.c_str(), 5);
  gapmer h(t.c_str(), 6, 5, 1);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour13) {
  std::string s = "ACATG";
  std::string t = "C...ACATG";
  gapmer g(s.c_str(), 5);
  gapmer h(t.c_str(), 6, 1, 3);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour14) {
  std::string s = "ACATG";
  std::string t = "T.CATG";
  gapmer g(s.c_str(), 5);
  gapmer h(t.c_str(), 5, 1, 1);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour15) {
  std::string s = "ACATG";
  std::string t = "G.....GTGTA";
  gapmer g(s.c_str(), 5);
  gapmer h(t.c_str(), 6, 1, 5);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour16) {
  std::string s = "ACATG";
  std::string t = "C.....GCCA";
  gapmer g(s.c_str(), 5);
  gapmer h(t.c_str(), 6, 1, 5);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour17) {
  std::string s = "ACATG";
  std::string t = "T....GATG";
  gapmer g(s.c_str(), 5);
  gapmer h(t.c_str(), 5, 1, 4);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour18) {
  std::string s = "AA.....AAA";
  std::string t = "AA......AAG";
  gapmer g(s.c_str(), 5, 2, 5);
  gapmer h(t.c_str(), 5, 2, 6);
  ASSERT_TRUE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_TRUE(h.is_neighbour(g));
}

TEST(gapmer, Neighbour19) {
  std::string s = "AAAAAA";
  std::string t = "AGA.AAAA";
  gapmer g(s.c_str(), 6);
  gapmer h(t.c_str(), 7, 3, 1);
  ASSERT_FALSE(g.is_neighbour(h)) << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour(g));
}

TEST(gapmer, RCNeighbour1) {
  std::string s = "CATTA...TTACC";
  gapmer g(s.c_str(), 10, 3, 5);
  gapmer o = g.reverse_complement();
  auto callback = [&](const decltype(o)& n) {
    ASSERT_TRUE(g.is_neighbour<true>(n));
  };
  o.huddinge_neighbours(callback);
}

TEST(gapmer, RCNeighbour2) {
  std::string s = "CAAT.TTAC";
  std::string t = "AAAA.AAAA";
  gapmer g(s.c_str(), 8, 1, 4);
  gapmer h(t.c_str(), 8, 1, 4);
  ASSERT_FALSE(g.is_neighbour<true>(h))
      << g.to_string() << " <-> " << h.to_string();
  ASSERT_FALSE(h.is_neighbour<true>(g));
}

TEST(gapmer, Huddinge1) {
  std::string s = "ACGTGT";
  typedef gapmer<false, 5> G;
  G g(s.c_str(), 6);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge2) {
  std::string s = "A.GTGT";
  typedef gapmer<false, 5> G;
  G g(s.c_str(), 5, 1, 1);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge3) {
  std::string s = "ACGT.T";
  typedef gapmer<false, 5> G;
  G g(s.c_str(), 5, 4, 1);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge4) {
  std::string s = "AC.TCT";
  typedef gapmer<false, 5> G;
  G g(s.c_str(), 5, 4, 1);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge5) {
  std::string s = "ACGTGT";
  typedef gapmer<true, 5> G;
  G g(s.c_str(), 6);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge6) {
  std::string s = "ACG.GT";
  typedef gapmer<true, 5> G;
  G g(s.c_str(), 5, 3, 1);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge7) {
  std::string s = "ACG.TGT";
  typedef gapmer<true, 5> G;
  G g(s.c_str(), 5, 3, 1);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge8) {
  std::string s = "AC..TCT";
  typedef gapmer<true, 5> G;
  G g(s.c_str(), 5, 2, 2);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge9) {
  std::string s = "GAC..TCT";
  typedef gapmer<true, 5> G;
  G g(s.c_str(), 6, 3, 2);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge10) {
  std::string s = "GAC..TC";
  typedef gapmer<true, 5> G;
  G g(s.c_str(), 5, 3, 2);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge11) {
  std::string s = "GACTC";
  typedef gapmer<true, 5> G;
  G g(s.c_str(), 5, 0, 0);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge12) {
  std::string s = "AC.TCT";
  typedef gapmer<true, 5> G;
  G g(s.c_str(), 5, 2, 1);
  bool res = sf::compare_generation<G, true, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge13) {
  std::string s = "A...CTCT";
  typedef gapmer<false, 5> G;
  G g(s.c_str(), 5, 1, 3);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge14) {
  std::string s = "AG.TCT";
  typedef gapmer<false, 5> G;
  G g(s.c_str(), 5, 2, 1);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge15) {
  std::string s = "AGTC...T";
  typedef gapmer<false, 5> G;
  G g(s.c_str(), 5, 4, 3);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, Huddinge16) {
  std::string s = "AG...TCT";
  typedef gapmer<false, 5> G;
  G g(s.c_str(), 5, 2, 3);
  bool res = sf::compare_generation<G, false, 5>(g);
  ASSERT_TRUE(res);
}

TEST(gapmer, IsCanonical1) {
  std::string s = "AAAA...TTT";
  // rc = AAA...TTTT so technically smaller but ignore gap
  gapmer<true, 5> g(s.c_str(), 7, 4, 3);
  ASSERT_TRUE(g.is_canonical());
}

TEST(gapmer, IsCanonical2) {
  std::string s = "ACGCCGT";
  gapmer<true, 5> g(s.c_str(), 7);
  ASSERT_TRUE(g.is_canonical());
}

TEST(gapmer, IsCanonical3) {
  std::string s = "ACGGCGT";
  gapmer<true, 5> g(s.c_str(), 7);
  ASSERT_FALSE(g.is_canonical());
}

TEST(gapmer, IsCanonical4) {
  std::string s = "AAATTT";
  gapmer<true, 5> g(s.c_str(), 6);
  ASSERT_TRUE(g.is_canonical());
}

TEST(gapmer, IsCanonical5) {
  std::string s = "GGGTTT";
  gapmer<true, 5> g(s.c_str(), 6);
  ASSERT_FALSE(g.is_canonical());
}

TEST(gapmer, RevComp1) {
  std::string s = "AAAATTTT";
  gapmer<true, 5> g(s.c_str(), 8);
  g = g.reverse_complement();
  ASSERT_EQ(g.to_string(), s);
}

TEST(gapmer, RevComp2) {
  std::string s = "AAA.TTTT";
  std::string rc = "AAAA.TTT";
  gapmer<true, 5> g(s.c_str(), 7, 3, 1);
  g = g.reverse_complement();
  ASSERT_EQ(g.to_string(), rc);
}

TEST(gapmer, RevComp3) {
  std::string s = "ACGTA";
  std::string rc = "TACGT";
  gapmer<true, 5> g(s.c_str(), 5);
  g = g.reverse_complement();
  ASSERT_EQ(g.to_string(), rc);
}

TEST(gapmer, RevComp4) {
  std::string s = "ACGT";
  gapmer<true, 5> g(s.c_str(), 4);
  g = g.reverse_complement();
  ASSERT_EQ(g.to_string(), s);
}

TEST(gapmer, Align1) {
  std::string s = "AAAAA";
  gapmer<true, 5> g(s.c_str(), 5);
  bool res = g.aligns_to(g);
  ASSERT_TRUE(res);
  gapmer<true, 5> o = g.reverse_complement();
  res = g.aligns_to(o);
  ASSERT_TRUE(res);
}

TEST(gapmer, Align2) {
  std::string s = "AAAAA";
  std::string os = "GGGGAAAAAT";
  gapmer<true, 5> g(s.c_str(), 5);
  gapmer<true, 5> o(os.c_str(), 10);
  bool res = g.aligns_to(o);
  ASSERT_TRUE(res);
}

TEST(gapmer, Align3) {
  std::string s = "AAAAA";
  std::string os = "GGGG...AAAAA";
  gapmer<true, 5> g(s.c_str(), 5);
  gapmer<true, 5> o(os.c_str(), 9, 4, 3);
  bool res = g.aligns_to(o);
  ASSERT_TRUE(res);
}

TEST(gapmer, Align4) {
  std::string s = "AAA...AA";
  std::string os = "GAAA...AATTT";
  gapmer<true, 5> g(s.c_str(), 5, 3, 3);
  gapmer<true, 5> o(os.c_str(), 9, 4, 3);
  bool res = g.aligns_to(o);
  ASSERT_TRUE(res);
}

}  // namespace sf
