# move to this directory
ORIGDIR=($pwd)
SCRIPT_DIR=$(cd $(dirname $0); pwd)
cd $SCRIPT_DIR

# We need to prepare output directory before running tests
TEST_OUTPUT_DIR=python/data
if [ ! -e $TEST_OUTPUT_DIR ]; then
        mkdir $TEST_OUTPUT_DIR
fi

# run pytest
pytest -v python --boxed

cd $ORIGDIR
