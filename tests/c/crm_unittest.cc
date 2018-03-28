#include "gtest/gtest.h"

extern "C"{
#include "crm.h"
}


namespace {

class CrmTest : public ::testing::Test {
protected:
  virtual void SetUp() {
    char filename[] = "reference_data/Li01b";
    int nele = 1;
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
    printf("%f", rate->dir);
    // Test SetTRRate works fine
    //EXPECT_NE(j, 3);  // you can also write assertion here
  };
};

TEST_F(CrmTest, TestCase1) {
  TRrateTest();
}

} // namespace
