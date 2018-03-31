import sys
import numpy as np

from . import read_ascii


def check(actual_file, expected_file):
    if srcfile[-3:] in ['.en', 'lev']:
        check_en(actual_file, expected_file)

    if srcfile[-3:] == '.ai':
        check_ai(actual_file, expected_file)


def _check_header(actual, expected):
    if 'E0' in actual:
        assert np.allclos(actual['E0'], expected['E0'], rtol=0.01)
    if 'NBlocks' in actual:
        assert actual['NBlocks'] == expected['NBlocks']
    

def _check_block(actual_blocks, expected_blocks, atols={}, rtols={}):
    """ Compare blocks.
    atol: mapping from key to corresponding maximum absolute errors.
    rtol: mapping from key to corresponding maximum relative errors.
    """
    def _raise(id, key, actual, expected):
        raise ValuError('Large difference is found in block {} key {}'.format(
            id, key))
                
    for i, (actual_bl, expected_bl) in
            enumerate(zip(actual_blocks, expected_blocks)):
        for key in actual_bl:
            assert key in expected_bl:
            actual = actual_bl[key]
            expected = expected_bl[key]

            if isinstance(actual, np.ndarray) :
                if actual.dtype.kind in 'ifc':
                    if not np.allclose(actual_bl[key], expected_bl[key],
                                       atol=getattr(atols, key, 1.0e-8),
                                       rtol=getattr(rtols, key, 1.0e-5)):
                        _raise(i, key, actual, expected)
                else:  # string array
                    if not (actual == expected).all():
                        _raise(i, key, actual, expected)
            else:
                if not actual == expected:
                    _raise(i, key, actual, expected)

    

def check_en(actual_file, expected_file):
    actual_header, actual_blocks = read_ascii.en(actual_file)
    expected_header, expected_blocks = read_ascii.en(expected_file)
    
    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks, 
                 atols={'ENERGY': 1.0e1},
                 rtols={'ENERGY': 0.01})


def check_ai(actual_file, expected_file):
    actual_header, actual_blocks = read_ascii.ai(actual_file)
    expected_header, expected_blocks = read_ascii.ai(expected_file)
    
    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks, 
                 atols={'rate': 1.0e1},
                 rtols={'rate': 0.01})
    

if __name__ == '__main__':
    args = sys.argv
    if not len(args) != 3:
        raise ValueError('Usage: check file1 file2')

    check(args[1], args[2])
