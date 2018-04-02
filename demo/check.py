import sys
import numpy as np
from collections import OrderedDict


def _read_value(lines, cls):
    """
    Pop one line from lines,
    split line by '=' and cast the value, then rturn the value and the rest of
    lines.
    """
    vals = lines[0].split('=')
    if len(vals) > 1:
        val = cls(vals[1].strip())
    else:
        val = ''
    return val, lines[1:]


def _get_header(lines):
    """ Read fac header """
    header = {}
    header['FAC'] = lines[0][4:-1]
    lines = lines[1:]
    header['Endian'], lines = _read_value(lines, int)
    header['TSess'], lines = _read_value(lines, int)
    header['Type'], lines = _read_value(lines, int)
    header['Verbose'], lines = _read_value(lines, int)
    key = lines[0].split('\t')[0]
    header[key], lines = _read_value(lines, float)
    header['NBlocks'], lines = _read_value(lines, int)
    return header, lines


def read_lev(filename):
    """ read *a.lev / *a.en file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    ind, e0 = lines[0].split('=')[-1].split(',')
    header['E0_index'] = int(ind)
    header['E0'] = float(e0)
    lines = lines[2:]

    def read_blocks(lines):
        block = {}
        block['NELE'], lines = _read_value(lines, int)
        nlev, lines = _read_value(lines, int)
        # read the values
        block['ILEV'] = np.zeros(nlev, dtype=int)
        block['IBASE'] = np.zeros(nlev, dtype=int)
        block['ENERGY'] = np.zeros(nlev, dtype=float)
        block['P'] = np.zeros(nlev, dtype=int)
        block['VNL'] = np.chararray(nlev, itemsize=3)
        block['2J'] = np.zeros(nlev, dtype=int)
        block['ncomplex'] = np.chararray(nlev, itemsize=20)
        block['sname'] = np.chararray(nlev, itemsize=20)
        block['name'] = np.chararray(nlev, itemsize=30)
        lines = lines[1:]
        
        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks
            block['ILEV'][i] = int(line[:6])
            block['IBASE'][i] = int(line[7:13])
            block['ENERGY'][i] = float(line[14:29])
            block['P'][i] = int(line[30:31])
            block['VNL'][i] = line[34:37].strip()
            block['2J'][i] = int(line[41:42])
            block['ncomplex'][i] = line[43:63].strip()
            block['sname'][i] = line[64:84].strip()
            block['name'][i] = line[85:].strip()

        return (block, )
            
    return header, read_blocks(lines)


def read_tr(filename):
    """ read *a.tr file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    header['Nblocks'] = _read_value(lines, int)
    lines = lines[1:]

    def read_blocks(lines):
        block = {}
        block['NELE'], lines = _read_value(lines, int)
        ntrans, lines = _read_value(lines, int)
        block['MULTIP'], lines = _read_value(lines, int)
        block['GAUGE'], lines = _read_value(lines, int)
        block['MODE'], lines = _read_value(lines, int)
        # read the values
        block['lower_index'] = np.zeros(ntrans, dtype=int)
        block['lower_2J'] = np.zeros(ntrans, dtype=int)
        block['upper_index'] = np.zeros(ntrans, dtype=int)
        block['upper_2J'] = np.zeros(ntrans, dtype=int)
        block['Delta E'] = np.zeros(ntrans, dtype=float)
        block['gf'] = np.zeros(ntrans, dtype=float)
        block['rate'] = np.zeros(ntrans, dtype=float)
        block['multipole'] = np.zeros(ntrans, dtype=float)
        
        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks
            block['upper_index'][i] = int(line[:7])
            block['upper_2J'][i] = int(line[8:10])
            block['lower_index'][i] = int(line[14:17])
            block['lower_2J'][i] = int(line[18:20])
            block['Delta E'][i] = float(line[21:34])
            block['gf'][i] = float(line[35:48])
            block['rate'][i] = float(line[49:62])
            block['multipole'][i] = float(line[63:76])
  
        return (block, )

    return header, read_blocks(lines)


def read_ai(filename):
    """ read *a.ai file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    header['Nblocks'] = _read_value(lines, int)
    lines = lines[1:]

    def read_blocks(lines):
        block = {}
        block['NELE'], lines = _read_value(lines, int)
        ntrans, lines = _read_value(lines, int)
        block['EMIN'], lines = _read_value(lines, float)
        negrid, lines = _read_value(lines, int)
        block['EGRID'] = np.zeros(negrid, dtype=float)
        for i in range(negrid):
            block['EGRID'][i] = float(lines.pop(0))
        # read the values
        block['bound_index'] = np.zeros(ntrans, dtype=int)
        block['bound_2J'] = np.zeros(ntrans, dtype=int)
        block['free_index'] = np.zeros(ntrans, dtype=int)
        block['free_2J'] = np.zeros(ntrans, dtype=int)
        block['Delta E'] = np.zeros(ntrans, dtype=float)
        block['AI rate'] = np.zeros(ntrans, dtype=float)
        block['DC strength'] = np.zeros(ntrans, dtype=float)
        
        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks
            block['bound_index'][i] = int(line[:7])
            block['bound_2J'][i] = int(line[8:10])
            block['free_index'][i] = int(line[14:17])
            block['free_2J'][i] = int(line[18:20])
            block['Delta E'][i] = float(line[21:32])
            block['AI rate'][i] = float(line[33:44])
            block['DC strength'][i] = float(line[45:56])
  
        return (block, )

    return header, read_blocks(lines)
    

def read_ce(filename):
    """ read *a.ce file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    lines = lines[1:]

    def read_blocks(lines):
        block = {}
        block['NELE'], lines = _read_value(lines, int)
        ntrans, lines = _read_value(lines, int)
        block['QKMODE'], lines = _read_value(lines, int)
        nparams, lines = _read_value(lines, int)
        block['MSUB'], lines = _read_value(lines, int)
        block['PWTYPE'], lines = _read_value(lines, int)
        ntegrid, lines = _read_value(lines, int)
        block['TEGRID'] = np.zeros(ntegrid, dtype=float)
        for i in range(ntegrid):
            block['TEGRID'][i] = float(lines.pop(0))
    
        block['TE0'], lines = _read_value(lines, float)
        block['ETYPE'], lines = _read_value(lines, int)
        negrid, lines = _read_value(lines, int)
        block['EGRID'] = np.zeros(negrid, dtype=float)
        for i in range(negrid):
            block['EGRID'][i] = float(lines.pop(0))
        block['UTYPE'], lines = _read_value(lines, int)
        nusr, lines = _read_value(lines, int)
        block['USR'] = np.zeros(nusr, dtype=float)
        for i in range(nusr):
            block['USR'][i] = float(lines.pop(0))

        block['lower_index'] = np.zeros(ntrans, dtype=int)
        block['lower_2J'] = np.zeros(ntrans, dtype=int)
        block['upper_index'] = np.zeros(ntrans, dtype=int)
        block['upper_2J'] = np.zeros(ntrans, dtype=int)
        block['Delta E'] = np.zeros(ntrans, dtype=float)
        block['bethe'] = np.zeros(ntrans, dtype=float)
        block['born'] = np.zeros((ntrans, 2), dtype=float)
        if block['MSUB']:
            block['collision strength'] = [None] * ntrans
            block['crosssection'] = [None] * ntrans
        else:
            block['collision strength'] = np.zeros((ntrans, nusr), dtype=float)
            block['crosssection'] = np.zeros((ntrans, nusr), dtype=float)
        
        if block['MSUB']:
            block['ratio collision strength'] = np.zeros(ntrans, dtype=float)
            
        nsub = np.zeros(ntrans, dtype=int)
        if block['QKMODE'] == 2:
            block['params'] == np.zeros((ntrans, 4), dtype=float)

        for tr in range(ntrans):
            line = lines[0]
            lines = lines[1:]
            block['lower_index'][tr] = int(line[:7].strip())
            block['lower_2J'][tr] = int(line[8:10].strip())
            block['upper_index'][tr] = int(line[15:17].strip())
            block['upper_2J'][tr] = int(line[18:20].strip())
            block['Delta E'][tr] = float(line[21:31].strip())
            nsub = int(line[32:])
            if block['MSUB']:
                block['collision strength'][tr] = np.zeros(
                    (nusr, nsub), dtype=float)
                block['crosssection'][tr] = np.zeros(
                    (nusr, nsub), dtype=float)

            line = lines[0]
            lines = lines[1:]
            block['bethe'][tr] = float(line[:12].strip())
            block['born'][tr, 0] = float(line[12:23].strip())
            block['born'][tr, 1] = float(line[24:36].strip())

            for sub in range(nsub):
                if block['MSUB']:
                    line = lines[0]
                    lines = lines[1:]
                    block['ratio collision strength'][tr] = float(line)                    
                    for i in range(nusr):
                        line = lines[0]
                        lines = lines[1:]
                        block['collision strength'][tr][i, sub] = float(
                            line[12:23])
                        block['crosssection'][tr][i, sub] = float(line[24:])
                    if sub < nsub - 1:
                        line = lines[0]
                        lines = lines[1:]  # skip separator -----
                    
                else:
                    for i in range(nusr):
                        line = lines[0]
                        lines = lines[1:]
                        block['collision strength'][tr, i] = float(line[12:23])
                        block['crosssection'][tr, i] = float(line[24:])                    
                    
        if len(lines) < 3:
            return (block, )

        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks

        raise ValueError('Bad file format.')

    return header, read_blocks(lines)


def read_ci(filename):
    """ read *a.ci file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    header['Nblocks'] = _read_value(lines, int)
    lines = lines[1:]

    def read_blocks(lines):
        block = {}
        block['NELE'], lines = _read_value(lines, int)
        ntrans, lines = _read_value(lines, int)
        block['QKMODE'], lines = _read_value(lines, int)
        nparams, lines = _read_value(lines, int)
        block['PWTYPE'], lines = _read_value(lines, int)
        ntegrid, lines = _read_value(lines, int)
        block['TEGRID'] = np.zeros(ntegrid, dtype=float)
        for i in range(ntegrid):
            block['TEGRID'][i] = float(lines.pop(0))
        block['ETYPE'], lines = _read_value(lines, int)
        negrid, lines = _read_value(lines, int)
        block['EGRID'] = np.zeros(negrid, dtype=float)
        for i in range(negrid):
            block['EGRID'][i] = float(lines.pop(0))
        block['UTYPE'], lines = _read_value(lines, int)
        nusr, lines = _read_value(lines, int)
        block['USR'] = np.zeros(nusr, dtype=float)
        for i in range(nusr):
            block['USR'][i] = float(lines.pop(0))

        # read the values
        block['bound_index'] = np.zeros(ntrans, dtype=int)
        block['bound_2J'] = np.zeros(ntrans, dtype=int)
        block['free_index'] = np.zeros(ntrans, dtype=int)
        block['free_2J'] = np.zeros(ntrans, dtype=int)
        block['Delta E'] = np.zeros(ntrans, dtype=float)
        block['Delta L'] = np.zeros(ntrans, dtype=int)
        block['parameters'] = np.zeros((ntrans, nparams), dtype=float)
        block['collision strength'] = np.zeros((ntrans, nusr), dtype=float)
        block['crosssection'] = np.zeros((ntrans, nusr), dtype=float)
        
        for tr in range(ntrans):
            line = lines[0]
            lines = lines[1:]
            block['bound_index'][tr] = int(line[:7])
            block['bound_2J'][tr] = int(line[8:10])
            block['free_index'][tr] = int(line[14:17])
            block['free_2J'][tr] = int(line[18:20])
            block['Delta E'][tr] = float(line[21:32])
            block['Delta L'][tr] = int(line[33:])
            block['parameters'][tr] = [float(l) for l in lines[0].split()]
            lines = lines[1:]
            for i in range(nusr):
                line = lines[0]
                lines = lines[1:]
                block['collision strength'][tr, i] = float(line[12:23])
                block['crosssection'][tr, i] = float(line[24:])
              
            if len(lines) < 3:
                return (block, )

        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks

        raise ValueError('Bad file format.')

    return header, read_blocks(lines)


def read_rr(filename):
    """ read *a.rr file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    header['Nblocks'] = _read_value(lines, int)
    lines = lines[1:]

    def read_blocks(lines):
        block = {}
        block['NELE'], lines = _read_value(lines, int)
        ntrans, lines = _read_value(lines, int)
        block['QKMODE'], lines = _read_value(lines, int)
        block['MULTIP'], lines = _read_value(lines, int)
        nparams, lines = _read_value(lines, int)
        ntegrid, lines = _read_value(lines, int)
        block['TEGRID'] = np.zeros(ntegrid, dtype=float)
        for i in range(ntegrid):
            block['TEGRID'][i] = float(lines.pop(0))
        block['ETYPE'], lines = _read_value(lines, int)
        negrid, lines = _read_value(lines, int)
        block['EGRID'] = np.zeros(negrid, dtype=float)
        for i in range(negrid):
            block['EGRID'][i] = float(lines.pop(0))
        block['UTYPE'], lines = _read_value(lines, int)
        nusr, lines = _read_value(lines, int)
        block['USR'] = np.zeros(nusr, dtype=float)
        for i in range(nusr):
            block['USR'][i] = float(lines.pop(0))

        # read the values
        block['bound_index'] = np.zeros(ntrans, dtype=int)
        block['bound_2J'] = np.zeros(ntrans, dtype=int)
        block['free_index'] = np.zeros(ntrans, dtype=int)
        block['free_2J'] = np.zeros(ntrans, dtype=int)
        block['Delta E'] = np.zeros(ntrans, dtype=float)
        block['Delta L'] = np.zeros(ntrans, dtype=int)
        block['parameters'] = np.zeros((ntrans, nparams), dtype=float)
        block['RR crosssection'] = np.zeros((ntrans, nusr), dtype=float)
        block['PI crosssection'] = np.zeros((ntrans, nusr), dtype=float)
        block['gf'] = np.zeros((ntrans, nusr), dtype=float)
        
        for tr in range(ntrans):
            line = lines[0]
            lines = lines[1:]
            block['bound_index'][tr] = int(line[:7])
            block['bound_2J'][tr] = int(line[8:10])
            block['free_index'][tr] = int(line[14:17])
            block['free_2J'][tr] = int(line[18:20])
            block['Delta E'][tr] = float(line[21:32])
            block['Delta L'][tr] = int(line[33:])
            block['parameters'][tr] = [float(l) for l in lines[0].split()]
            lines = lines[1:]
            for i in range(nusr):
                line = lines[0]
                lines = lines[1:]
                block['RR crosssection'][tr, i] = float(line[12:23])
                block['PI crosssection'][tr, i] = float(line[24:35])
                block['gf'][tr, i] = float(line[36:])
              
            if len(lines) < 3:
                return (block, )

        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks

        raise ValueError('Bad file format.')

    return header, read_blocks(lines)


def read_sp(filename):
    """ read *a.sp file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    header['Nblocks'] = _read_value(lines, int)
    lines = lines[1:]

    def read_blocks(lines):
        block = {}
        block['NELE'], lines = _read_value(lines, int)
        ntrans, lines = _read_value(lines, int)
        block['TYPE'], lines = _read_value(lines, str)
        block['IBLK'], lines = _read_value(lines, int)
        block['ICOMP'], lines = _read_value(lines, str)
        block['FBLK'], lines = _read_value(lines, int)
        block['FCOMP'], lines = _read_value(lines, str)

        # read the values
        block['block'] = np.zeros(ntrans, dtype=int)
        block['level'] = np.zeros(ntrans, dtype=int)
        block['abs. energy'] = np.zeros(ntrans, dtype=float)
        block['population'] = np.zeros(ntrans, dtype=float)
        block['Delta E'] = np.zeros(ntrans, dtype=float)
        block['emissivity'] = np.zeros(ntrans, dtype=float)
        
        for tr in range(ntrans):
            line = lines[0]
            lines = lines[1:]

            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks

            block['block'][tr] = int(line[:7])
            block['level'][tr] = int(line[8:14])
            block['abs. energy'][tr] = float(line[15:28])
            block['population'][tr] = float(line[29:40])
            block['Delta E'][tr] = float(line[15:28])
            block['emissivity'][tr] = float(line[29:40])
              
        return (block, )

    return header, read_blocks(lines)


def read_rt(filename):
    """ read *a.rt file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    header['Nblocks'] = _read_value(lines, int)
    lines = lines[1:]

    def read_blocks(lines):
        block = {}
        ntrans, lines = _read_value(lines, int)
        block['EDEN'], lines = _read_value(lines, float)
        block['EDIST'], lines = _read_value(lines, int)
        npedis, lines = _read_value(lines, int)
        block['EDIS'] = np.zeros(npedis, float)
        for i in range(npedis):
            line = lines[0]
            lines = lines[1:]
            block['EDIS'][i] = float(line)
        block['PDEN'], lines = _read_value(lines, float)
        block['PDIST'], lines = _read_value(lines, int)
        nppdis, lines = _read_value(lines, int)
        block['PPDIS'] = np.zeros(nppdis, float)
        for i in range(npedis):
            line = lines[0]
            lines = lines[1:]
            block['PPDIS'][i] = float(line)
        lines = lines[1:] # skip header

        # read the values
        block['block'] = np.zeros(ntrans, dtype=int)
        block['level'] = np.zeros(ntrans, dtype=int)
        block['NB'] = np.zeros(ntrans, dtype=int)
        block['TR'] = np.zeros(ntrans, dtype=int)
        block['CE'] = np.zeros(ntrans, dtype=int)
        block['RR'] = np.zeros(ntrans, dtype=int)
        block['AI'] = np.zeros(ntrans, dtype=int)
        block['CI'] = np.zeros(ntrans, dtype=int)
        block['ncomplex'] = np.chararray(ntrans, itemsize=20)
        
        for tr in range(ntrans):
            line = lines[0]
            lines = lines[1:]

            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks

            block['block'][tr] = int(line[:7])
            block['level'][tr] = int(line[8:12])
            block['NB'][tr] = float(line[13:25])
            block['TR'][tr] = float(line[26:37])
            block['CE'][tr] = float(line[38:49])
            block['RR'][tr] = float(line[50:61])
            block['AI'][tr] = float(line[62:73])
            block['CI'][tr] = float(line[74:85])
            block['ncomplex'][tr] = line[86:].strip()
              
        return (block, )

    return header, read_blocks(lines)


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
                    if not (actual == expected).all():
                        _raise(i, key, actual, expected)

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
                    _raise(i, key, actual, expected)


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
    
    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks, 
                 atols={'rate': 1.0e1},
                 rtols={'rate': 0.01})
    

def check_tr(actual_file, expected_file):
    actual_header, actual_blocks = read_tr(actual_file)
    expected_header, expected_blocks = read_tr(expected_file)
    
    _check_header(actual_header, expected_header)
    _check_block(actual_blocks, expected_blocks, 
                 atols={'Delta E': 1.0e0, 'gf': 1.0e-4, 'rate': 1.0e4, 
                        'multipole': 1.0e-5},
                 rtols={'Delta E': 0.01, 'gf': 0.02, 'rate': 0.05, 
                        'multipole': 0.01})
    

def check_ce(actual_file, expected_file):
    actual_header, actual_blocks = read_ce(actual_file)
    expected_header, expected_blocks = read_ce(expected_file)
    
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
