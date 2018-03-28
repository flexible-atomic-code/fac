import os

from pfac.atom import atomic_data
from pfac.spm import *
from pfac import const
from . import read_ascii


THIS_DIR = os.path.abspath(os.path.join(__file__, os.pardir))
# generate an output directory for this test
dir0 = THIS_DIR + '/data_{}/'.format(__name__)
if not os.path.exists(dir0):
    os.mkdir(dir0)
# directory that contains reference data. see run_pytest.sh
ref_dir = THIS_DIR + '../reference_data/'

# test with hydrogen like Lithium
asym = 'Be'

neles = list(range(3, 4))
# This generates several files below dir0
atomic_data(neles, asym, iprint=1, dir=dir0)


def assert_header_identical(actual, expected):
    assert actual['NBLOCKS'] == expected['NBLOCKS']


# validate en file
header_actual, blocks_actual = read_ascii.en(dir0 + 'Be03a.en')
header_expected, blocks_expected = read_ascii.en(dir0 + 'Be03a.en')
assert_header_identical(header_actual, header_expected)

for bl_actual, bl_expected in zip(blocks_actual, blocks_expected):
    for key in ['NELE', 'TYPE', 'IBLK', 'ICOMP', 'FBLK', 'FCOMP']:
        assert bl_actual[key] == bl_expected[key]
    for key in ['upper', 'lower']:
        assert (bl_actual[key] == bl_expected[key])
    assert np.allclose(bl_actual['Delta E'], bl_expected['Delta E'])
    assert np.allclose(bl_actual['emissivity'], bl_expected['emissivity'])

# validate tr file
