#include "gtest/gtest.h"

extern "C"{
#include "dbase.h"
}

namespace {

TEST(FindLevelByName, First) {
  char filename[] = "reference_data/li.lev.b";
  int nele = 3;
  int i;

  i = FindLevelByName(filename, nele, "1*2 2*1", "2s1", "2s+1(1)1");
  EXPECT_EQ(i, 0);

  i = FindLevelByName(filename, nele, "1*2 2*1", "2p1", "2p+1(3)3");
  EXPECT_EQ(i, 1);

  i = FindLevelByName(filename, nele, "1*2 2*1", "2p1", "2p-1(1)1");
  EXPECT_EQ(i, 2);
}
}
