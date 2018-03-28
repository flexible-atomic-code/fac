# Simple functions to read ascii files.
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


def lev(filename):
    """ read *.lev file """
    

def en(filename):
    """ read *a.en file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
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
        # convert to array
        block = {key: np.full(ntrans, val) for key, val in block.items()}

        # read the values
        block['upper'] = np.zeros(ntrans, dtype=int)
        block['lower'] = np.zeros(ntrans, dtype=int)
        block['Delta E'] = np.zeros(ntrans, dtype=float)
        block['emissivity'] = np.zeros(ntrans, dtype=float)

        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:].strip())
                return (block, ) + blocks
            block['upper'][i] = int(line[:7])
            block['lower'][i] = int(line[8:14])
            block['Delta E'][i] = float(line[15:27])
            block['emissivity'][i] = float(line[28:40])

        return (block, )

    return header, read_blocks(lines)


def ce(filename):
    """ read *a.en file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(f)
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
            block['TEGRID'][i] = float(lines[0].strip())
            lines = lines[1:]
        block['TE0'], lines = _read_value(lines, float)
        block['ETYPE'], lines = _read_value(lines, int)
        negrid, lines = _read_value(lines, int)
        block['EGRID'] = np.zeros(negrid, dtyp=float)
        for i in range(ntegrid):
            block['EGRID'][i] = float(lines[0].strip())
            lines = lines[1:]
        block['UTYPE'], lines = _read_value(lines, int)
        nusr, lines = _read_value(lines, int)
        block['USR'] = np.zeros(nusr, dtyp=float)
        for i in range(ntegrid):
            block['USR'][i] = float(lines[0].strip())
            lines = lines[1:]

        block['upper'] = np.zeros(ntrans, dtype=int)
        block['lower'] = np.zeros(ntrans, dtype=int)
        nsub = np.zeros(ntrans, dtype=int)
        block['bethe'] = np.zeros(ntrans, dtype=float)
        block['born'] = np.zeros((ntrans, 2), dtype=float)
        # TODO this part looks different from the description in the manual
        block['strength'] = np.zeros((ntrans, nusr, 3), dtype=float)

        for tr in range(ntrans):
            line = lines[0]
            lines = lines[1:]
            block['upper'][tr] = int(line[:6].strip())
            block['lower'][tr] = int(line[7:11].strip())
            nsub = int(line[11:16].strip())
            block['bethe'][tr] = float(line[17:19].strip())
            block['born'][tr, 0] = float(line[20:31].strip())
            block['born'][tr, 1] = float(line[32:].strip())
            for i in range(nusr):
                line = lines[0]
                lines = lines[1:]
                block['strength'][tr, i, 0] = float(line[:11])
                block['strength'][tr, i, 1] = float(line[12:23])
                block['strength'][tr, i, 2] = float(line[24:])

        if len(lines) < 3:
            return (block, )

        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks

        raise ValueError('Bad file format.')

    return header, read_blocks(lines)
