# move to this directory
ORIGDIR=($pwd)
SCRIPT_DIR=$(cd $(dirname $0); pwd)
cd $SCRIPT_DIR

GTEST_DIR=$SCRIPT_DIR/googletest/googletest


cd $ORIGDIR
