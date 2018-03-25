#include "gtest/gtest.h"

extern "C"{
#include "dbase.h"
}

namespace {

TEST(FindLevelByName, First) {
  // This test is named "Negative", and belongs to the "FactorialTest"
  // test case.
  char filename[] = "c/examples/li.lev.b";
  int nele = 3;
  char nc[] = "1*2 2*1";
  char cnr[] = "2s1";
  char cr[] = "2s+1(1)1";
  int i;

  i = FindLevelByName(filename, nele, "1*2 2*1", "2s1", "2s+1(1)1");
  EXPECT_EQ(i, 0);

  i = FindLevelByName(filename, nele, "1*2 2*1", "2p1", "2p+1(3)3");
  EXPECT_EQ(i, 1);

  i = FindLevelByName(filename, nele, "1*2 2*1", "2p1", "2p-1(1)1");
  EXPECT_EQ(i, 2);
}
}
