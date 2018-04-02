#include "gtest/gtest.h"

extern "C"{
#include "crm.h"
}


namespace {

class CrmTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    char filename[] = "c/examples/Li02b";
    int nele = 2;
    double n = 1.0;
    InitCRM();
    AddIon(nele, n, filename);
    SetBlocks(0.0, NULL);
    SetTRRates(0.0);
  };

  void TRrateTest() {
    ARRAY* ions;
    ION* ion;
    BLK_RATE* blk_rate;
    RATE* rate;

    ions = _GetIons();
    ion = (ION *) ArrayGet(ions, 0);
    blk_rate = (BLK_RATE *)ArrayGet(ion->tr_rates, 0);
    rate = (RATE*) ArrayGet(blk_rate->rates, 0);
    EXPECT_NE(rate->dir, 5829.239544);
    // Test SetTRRate works fine
    //EXPECT_NE(j, 3);  // you can also write assertion here
  };
};

TEST_F(CrmTest, TestCase1) {
  TRrateTest();
}

} // namespace
