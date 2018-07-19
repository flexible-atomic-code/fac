import sys
import numpy as np
from collections import OrderedDict
sys.path = sys.path + ['../python','../../python']
from rfac import *

def _check_header(actual, expected):
    if 'E0' in actual:
        assert np.allclose(actual['E0'], expected['E0'], rtol=0.01)
    if 'NBlocks' in actual:
        assert actual['NBlocks'] == expected['NBlocks']


def _check_block(actual_blocks, expected_blocks, atols={}, rtols={}):
    """ Compare each values in two blocks.
    atol: mapping from key to corresponding maximum absolute errors.
    rtol: mapping from key to corresponding maximum relative errors.
    """
    def _raise(id, key, actual, expected):
        message = 'Large difference is found in block {} key {}.'.format(
            id, key)
        if isinstance(actual, np.ndarray):
            ind_a = np.unravel_index(np.argmax(np.abs(actual - expected)),
                                     actual.shape)
            ind_r = np.unravel_index(
                np.argmax(0.5 * np.abs(actual - expected) /
                    (np.abs(actual) + np.abs(expected))), actual.shape)
            message += ' The largest abs difference is in {}'.format(ind_a)
            message += ' with {} vs {} '.format(actual[ind_a], expected[ind_a])
            message += ' The largest rel difference is in {}'.format(ind_r)
            message += ' with {} vs {} '.format(actual[ind_r], expected[ind_r])
        else:
            message += ' {} vs {}'.format(actual, expected)
        raise ValueError(message)

    for i, (actual_bl, expected_bl) in  enumerate(
            zip(actual_blocks, expected_blocks)):
        for key in actual_bl:
            assert key in expected_bl
            actual = actual_bl[key]
            expected = expected_bl[key]

            if isinstance(actual, np.ndarray) :
                if actual.dtype.kind in 'ifc':
                    if not np.allclose(actual_bl[key], expected_bl[key],
                                       atol=atols.get(key, 1.0e-8),
                                       rtol=rtols.get(key, 1.0e-5)):
                        _raise(i, key, actual, expected)
                else:  # string array
                    pass

            elif isinstance(actual, list):  # list of np.ndarray
                for ac, ex in zip(actual_bl[key], expected_bl[key]):
                    if not np.allclose(ac, ex,
                                       atol=atols.get(key, 1.0e-8),
                                       rtol=rtols.get(key, 1.0e-5)):
                        _raise(i, key, actual, expected)

            elif isinstance(actual, float):
                if not np.allclose(actual_bl[key], expected_bl[key],
                                   atol=atols.get(key, 1.0e-8),
                                   rtol=rtols.get(key, 1.0e-5)):
                   _raise(i, key, actual, expected)
            else:
                if not actual == expected:
                    pass


def _sort_array(arrays, ref_arrays, *keys):
    """ Sort arrays so that arrays[key] == ref_arrays[key] """
    size = len(ref_arrays[keys[0]])

    new_indexes = np.zeros(size, dtype=int)
    for j in range(size):
        for i in range(size):
            if all(ref_arrays[k][i] == arrays[k][j] for k in keys):
                new_indexes[i] = j
                break

    new_arrays = OrderedDict()
    for k, v in arrays.items():
        if isinstance(v, np.ndarray) and len(v) == size:
            new_arrays[k] = v[new_indexes]
        elif isinstance(v, list) and len(v) == size:
            new_arrays[k] = []
            for ind in new_indexes:
                new_arrays[k].append(v[ind])
        else:
            new_arrays[k] = v

    return new_arrays


def check_en(actual_file, expected_file):
    actual_header, actual_blocks = read_lev(actual_file)
    expected_header, expected_blocks = read_lev(expected_file)

    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks,
                 atols={'ENERGY': 1.0e1},
                 rtols={'ENERGY': 0.01})


def check_ai(actual_file, expected_file):
    actual_header, actual_blocks = read_ai(actual_file)
    expected_header, expected_blocks = read_ai(expected_file)

    # we need to sort arrays, because with openmp, the order may differ
    actual_blocks = [_sort_array(ac, ex, 'bound_index', 'free_index')
                     for ac, ex in zip (actual_blocks, expected_blocks)]

    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks,
                 atols={'rate': 1.0e1, 'Delta E': 0.5},
                 rtols={'rate': 0.01, 'DC strength': 0.01,
                        'AI rate': 0.01, 'Delta E': 0.01,
                        'EGRID': 1.0e-4})


def check_tr(actual_file, expected_file):
    actual_header, actual_blocks = read_tr(actual_file)
    expected_header, expected_blocks = read_tr(expected_file)

    # we need to sort arrays, because with openmp, the order may differ
    actual_blocks = [_sort_array(ac, ex, 'lower_index', 'upper_index')
                     for ac, ex in zip (actual_blocks, expected_blocks)]

    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks,
                 atols={'Delta E': 1.0e0, 'gf': 1.0e-4, 'rate': 1.0e4,
                        'multipole': 1.0e-5},
                 rtols={'Delta E': 0.01, 'gf': 0.02, 'rate': 0.05,
                        'multipole': 0.01})


def check_ce(actual_file, expected_file):
    actual_header, actual_blocks = read_ce(actual_file)
    expected_header, expected_blocks = read_ce(expected_file)

    # we need to sort arrays, because with openmp, the order may differ
    actual_blocks = [_sort_array(ac, ex, 'lower_index', 'upper_index')
                     for ac, ex in zip (actual_blocks, expected_blocks)]

    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks,
                 atols={'Delta E': 1.0e0, 'bethe': 1.0e0, 'born': 1.0,
                        'collision strength': 0.001, 'crosssection': 0.01,
                        'TEGRID': 1.0, 'EGRID': 1.0, 'USR': 1.0,
                        'TE0': 1.0},
                 rtols={'Delta E': 0.01, 'bethe': 0.01, 'born': 0.01,
                        'collision strength': 0.02, 'crosssection': 0.1,
                        'TEGRID': 0.01, 'EGRID': 0.01, 'USR': 0.01,
                        'TE0': 0.01})


def check_ci(actual_file, expected_file):
    actual_header, actual_blocks = read_ci(actual_file)
    expected_header, expected_blocks = read_ci(expected_file)

    # we need to sort arrays, because with openmp, the order may differ
    actual_blocks = [_sort_array(ac, ex, 'bound_index', 'free_index')
                     for ac, ex in zip (actual_blocks, expected_blocks)]

    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks,
                 atols={'Delta E': 1.0e0, 'bethe': 1.0e0, 'born': 1.0,
                        'collision strength': 0.001, 'crosssection': 0.01,
                        'TEGRID': 1.0, 'EGRID': 1.0, 'USR': 1.0,
                        'TE0': 1.0},
                 rtols={'Delta E': 0.01, 'bethe': 0.01, 'born': 0.01,
                        'collision strength': 0.02, 'crosssection': 0.1,
                        'TEGRID': 0.01, 'EGRID': 0.01, 'USR': 0.01,
                        'TE0': 0.01, 'parameters': 0.01})


def check_rr(actual_file, expected_file):
    actual_header, actual_blocks = read_rr(actual_file)
    expected_header, expected_blocks = read_rr(expected_file)

    # we need to sort arrays, because with openmp, the order may differ
    actual_blocks = [_sort_array(ac, ex, 'bound_index', 'free_index')
                     for ac, ex in zip (actual_blocks, expected_blocks)]

    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks,
                 atols={'Delta E': 1.0e0, 'bethe': 1.0e0, 'born': 1.0,
                        'TEGRID': 1.0, 'EGRID': 1.0, 'USR': 1.0,
                        'TE0': 1.0},
                 rtols={'Delta E': 0.01, 'bethe': 0.01, 'born': 0.01,
                        'RR crosssection': 0.01, 'PI crosssection': 0.01,
                        'TEGRID': 0.01, 'EGRID': 0.01, 'USR': 0.01,
                        'TE0': 0.01, 'parameters': 0.01, 'gf': 0.01})


def check_sp(actual_file, expected_file):
    actual_header, actual_blocks = read_sp(actual_file)
    expected_header, expected_blocks = read_sp(expected_file)

    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks,
                 atols={'Delta E': 1.0e0, 'emissivity': 1.0e-10},
                 rtols={'Delta E': 0.01, 'emissivity': 0.01})


def check_rt(actual_file, expected_file):
    actual_header, actual_blocks = read_rt(actual_file)
    expected_header, expected_blocks = read_rt(expected_file)

    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks)


def check(actual_file, expected_file):
    if actual_file[-3:] in ['.en', 'lev']:
        check_en(actual_file, expected_file)

    elif actual_file[-3:] == '.ai':
        check_ai(actual_file, expected_file)

    elif actual_file[-3:] == '.tr':
        check_tr(actual_file, expected_file)

    elif actual_file[-3:] in ['.ce', 'ceM']:
        check_ce(actual_file, expected_file)

    elif actual_file[-3:] == '.ci':
        check_ci(actual_file, expected_file)

    elif actual_file[-3:] == '.rr':
        check_rr(actual_file, expected_file)

    elif actual_file[-3:] == '.sp':
        check_sp(actual_file, expected_file)

    elif actual_file[-3:] == '.rt':
        check_rt(actual_file, expected_file)

    else:
        raise TypeError('Unknown file extension: {}'.format(actual_file))


if __name__ == '__main__':
    args = sys.argv
    if len(args) != 3:
        raise ValueError('Usage: python check.py file1 file2')

    print(args[1])
    check(args[1], args[2])
