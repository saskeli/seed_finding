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

}  // namespace sf
