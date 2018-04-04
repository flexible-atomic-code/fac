# move to this directory
ORIGDIR=($pwd)
SCRIPT_DIR=$(cd $(dirname $0); pwd)
cd $SCRIPT_DIR
# install googletest
git clone https://github.com/google/googletest

# compile google test
# GTEST_DIR=$SCRIPT_DIR/googletest/googletest
# g++ -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
#     -pthread -c ${GTEST_DIR}/src/gtest-all.cc
# ar -rv libgtest.a gtest-all.o

cd $ORIGDIR
