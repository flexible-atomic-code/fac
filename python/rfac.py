#
#   FAC - Flexible Atomic Code
#   Copyright (C) 2001-2015 Ming Feng Gu
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>.
#

# the read_*** functions were written by Keisuke Fujii of Kyoto Univ.

import numpy as np
from collections import OrderedDict
from distutils.version import LooseVersion


def _wrap_get_length(line0):
    """ Returns get_length functions for lev file, depending on the version """
    
    slice_lcomplex = slice(43, 75)
    slice_lsname = slice(76, 124)
    slice_lname = slice(125, None)
    if line0 != None and len(line0) > 65:
        n = line0[64:65]
        s = line0[65:66]
        if n.isdigit() and (s.isalpha() or s == '['):
            slice_lcomplex = slice(43, 64)
            slice_lsname = slice(64, 85)
            slice_lname = slice(85, None)

    def get_lcomplex(line):
        return line[slice_lcomplex].strip()

    def get_lsname(line):
        return line[slice_lsname].strip()

    def get_lname(line):
        return line[slice_lname].strip()

    return get_lcomplex, get_lsname, get_lname


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
    if len(lines) > 3:
        line0 = lines[3]
    else:
        line0 = None
    get_lcomplex, get_lsname, get_lname = _wrap_get_length(line0)

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
        block['ncomplex'] = np.chararray(nlev, itemsize=32)
        block['sname'] = np.chararray(nlev, itemsize=48)
        block['name'] = np.chararray(nlev, itemsize=128)
        lines = lines[1:]

        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks
            block['ILEV'][i] = int(line[:6])
            block['IBASE'][i] = int(line[7:13])
            block['ENERGY'][i] = float(line[14:30])
            block['P'][i] = int(line[30:31])
            block['VNL'][i] = line[32:37].strip()
            block['2J'][i] = int(line[38:42])
            block['ncomplex'][i] = get_lcomplex(line)
            block['sname'][i] = get_lsname(line)
            block['name'][i] = get_lname(line)

        return (block, )

    return header, read_blocks(lines)


def read_en(filename):
    return read_lev(filename)


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
            block['lower_index'][i] = int(line[11:17])
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
            block['free_index'][i] = int(line[11:17])
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
            block['upper_index'][tr] = int(line[11:17].strip())
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
            block['free_index'][tr] = int(line[11:17])
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
            block['free_index'][tr] = int(line[11:17])
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


class FLEV:
    def sort(self):
        i = np.argsort(self.e)
        self.e = self.e[i]
        self.p = self.p[i]
        self.j = self.j[i]
        self.c = self.c[i]
        self.v = self.v[i]
        self.s = self.s[i]
        self.n = self.n[i]
        self.ig = self.ig[i]

    def add(self, a):
        r = FLEV(None)
        r.concat(self, a)
        return r

    def concat(self, g, c):
        ng = len(g.e)
        nc = len(c.e)
        n = ng+nc
        self.e0 = min(g.e0, c.e0)
        self.e = np.zeros(n, dtype=float)
        self.p = np.zeros(n, dtype=int)
        self.j = np.zeros(n, dtype=int)
        self.c = np.chararray(n, itemsize=32)
        self.v = np.chararray(n, itemsize=3)
        self.s = np.chararray(n, itemsize=48)
        self.n = np.chararray(n, itemsize=128)
        self.ig = np.zeros(n, dtype=int)
        self.e[:ng] = g.e
        self.e[ng:] = c.e
        self.p[:ng] = g.p
        self.p[ng:] = c.p
        self.j[:ng] = g.j
        self.j[ng:] = c.j
        self.c[:ng] = g.c
        self.c[ng:] = c.c
        self.v[:ng] = g.v
        self.v[ng:] = c.v
        self.s[:ng] = g.s
        self.s[ng:] = c.s
        self.n[:ng] = g.n
        self.n[ng:] = c.n
        self.ig[:ng] = g.ig
        self.ig[ng:] = c.ig

    def combine(self, g, c):
        wg = np.where(g.ig == 1)
        wc = np.where(c.ig == 0)
        ng = len(wg[0])
        nc = len(wc[0])
        n = ng+nc
        self.eg = g.e0
        self.ec = c.e0
        self.e0 = g.e0
        self.e = np.zeros(n, dtype=float)
        self.p = np.zeros(n, dtype=int)
        self.j = np.zeros(n, dtype=int)
        self.c = np.chararray(n, itemsize=32)
        self.v = np.chararray(n, itemsize=3)
        self.s = np.chararray(n, itemsize=48)
        self.n = np.chararray(n, itemsize=128)
        self.ig = np.zeros(n, dtype=int)
        self.e[:ng] = g.e[wg]
        self.e[ng:] = c.e[wc]
        self.e = self.e - self.e0
        self.p[:ng] = g.p[wg]
        self.p[ng:] = c.p[wc]
        self.j[:ng] = g.j[wg]
        self.j[ng:] = c.j[wc]
        self.c[:ng] = g.c[wg]
        self.c[ng:] = c.c[wc]
        self.v[:ng] = g.v[wg]
        self.v[ng:] = c.v[wc]
        self.s[:ng] = g.s[wg]
        self.s[ng:] = c.s[wc]
        self.n[:ng] = g.n[wg]
        self.n[ng:] = c.n[wc]
        self.ig[:ng] = g.ig[wg]
        self.ig[ng:] = c.ig[wc]

    def __init__(self, f, ig=0):
        if f == None:
            return
        (hlev,blev) = read_lev(f)
        b0 = blev[0]
        self.nele = b0['NELE']
        ks = {'s':0, 'p':1, 'd':2, 'f':3, 'g':4, 'h':5}
        if self.nele <= 12:
            gc = [2, 2, 6, self.nele-10, 0, 0]
        else:
            gc = [2, 2, 6, 2, self.nele-12, 0]
        for m in range(4,6):
            for k in range(m):
                gc.append(0)
        if ig > 0:
            n = f.split('/')
            n = n[-2]
            cc = gc+[]
            if n != 'g':
                n1 = int(n[3:4])
                k1 = int(n[5:6])
                n2 = int(n[8:9])
                k2 = int(n[10:11])
                i1 = n1*(n1-1)/2 + k1
                i2 = n2*(n2-1)/2 + k2
                cc[i1] -= 1
                cc[i2] += 1
            self.cc = cc
        self.e0 = hlev['E0']
        self.e = b0['ENERGY']
        self.p = b0['P']
        self.j = b0['2J']
        self.c = b0['ncomplex']
        self.v = b0['VNL']
        self.s = b0['sname']
        self.n = b0['name']
        self.e = self.e + self.e0
        self.ig = np.zeros(len(self.e), dtype=int)
        if ig > 0:
            nq = cc+[]
            for i in range(len(self.e)):
                x = self.s[i].split('.')
                for s in range(len(nq)):
                    nq[s] = 0
                for y in x:
                    n = int(y[0:1])
                    k = ks[y[1:2]]
                    s = n*(n-1)/2 + k
                    nq[s] = int(y[2:])
                self.ig[i] = 1
                for s in range(len(nq)):
                    if nq[s] != cc[s]:
                        self.ig[i] = 0
                        break

    def match(self, m):
        self.idx = np.arange(len(self.e))
        self.em = np.zeros(len(self.e), dtype=float)
        self.em[:] = -1.0
        self.im = np.zeros(len(self.e), dtype=int)
        self.cm = np.chararray(len(self.e),itemsize=32)
        self.cm[:]=''
        if m == None:
            return
        for p in [0, 1]:
            for j in range(min(m.j),max(m.j)+1,2):
                w0 = np.where(np.logical_and(self.p == p, self.j == j))
                n0 = len(w0[0])
                if (n0 == 0):
                    continue
                ew0 = self.e[w0]
                w1 = np.where(np.logical_and(m.p == p, m.j == j))
                ew1 = m.e[w1]
                n1 = len(w1[0])
                for i in range(n1):
                    ade = abs(ew0 - ew1[i])
                    wm0 = np.argmin(ade)
                    if ade[wm0] > 3:
                        continue
                    wm = w0[0][wm0]
                    if self.em[wm] < 0 or np.fabs(self.em[wm]-self.e[wm])>ade[wm0]:
                        self.im[wm] = w1[0][i]
                        self.em[wm] = ew1[i]
                        self.cm[wm] = m.c[w1[0][i]]

    def write(self, fn):
        f = open(fn, 'w')
        for i in range(len(self.e)):
            em = self.em[i]
            if em >= 0:
                de = self.e[i]-self.em[i]
            else:
                em = 0.0
                de = 0.0
            s = '%4d %4d %11.5E %11.5E %10.3E %d %2d %-32s %-84s %-48s\n'%(i,self.im[i],self.e[i],em,de,self.p[i],self.j[i],self.s[i],self.n[i],self.cm[i])
            f.write(s)
        f.close()

class MLEV:
    def __init__(self, f):
        if f == None:
            return
        self.e = np.transpose(np.loadtxt(f, usecols=3, dtype='float', skiprows=1, delimiter=' ; '))
        j = np.transpose(np.loadtxt(f, usecols=1, dtype='float', skiprows=1, delimiter=' ; ', converters={1:lambda x: eval('2*'+x)}))
        self.j = np.int32(j)
        pc = np.transpose(np.loadtxt(f, usecols=2, dtype='string', skiprows=1, delimiter=' ; '))
        self.p = np.int32(pc == 'o')
        self.c = np.transpose(np.loadtxt(f, usecols=0, dtype='string', skiprows=1, delimiter=' ; '))

def aflev(d0, d1, a, n):
    if (d0 != None and len(d0) > 0):
        r0 = cflev(d0, a, 0, n)
    else:
        r0 = None
    if (d1 != None and len(d1) > 0):
        r1 = cflev(d1, a, 12, n)
    else:
        r1 = None
    if r0 == None:
        return r1
    if r1 == None:
        return r0
    r = FLEV(None)
    r.concat(r0, r1)
    r.sort()
    return r

def cflev(d, a, nj, n):
    if nj == 0:
        return pjflev(d, a, 'a', 'a', n)
    f0 = None
    for p in [0, 1]:
        for j in range(nj):
            j2 = j*2
            if f0 == None:
                f0 = pjflev(d, a, p, j2, n)
            else:
                f1 = pjflev(d, a, p, j2, n)
                if f1 != None:
                    fc = FLEV(None)
                    fc.concat(f0, f1)
                    f0 = fc
    if f0 != None:
        f0.sort()
    return f0

def pjflev(d, a, p, j, n):
    if (type(p) == type(0)):
        p = '%d'%p
        j = '%d'%j

    f = '%s/p%sj%s/g/%s%02di00a.en'%(d,p, j, a, n)
    try:
        f0 = FLEV(f)
    except:
        f0 = None
    for n1 in [2, 3]:
        for k1 in range(0, n1):
            for n2 in [3, 4, 5]:
                for k2 in range(0, n2):
                    f = '%s/p%sj%s/g_n%dk%d_n%dk%d/%s%02di00a.en'%(d,p,j,n1,k1,n2,k2,a,n)
                    if f0 == None:
                        try:
                            f0 = FLEV(f)
                        except:
                            pass
                    else:
                        try:
                            f1 = FLEV(f)
                            fc = FLEV(None)
                            fc.concat(f0, f1)
                            f0 = fc
                        except:
                            pass
    if f0 != None:
        f0.sort()
    return f0

def mflev(a, n, df0, df1, dm, fw):
    #dm = '/Users/yul20/atomic/juan/mrmp_results'
    #df0 = '/Users/yul20/atomic/juan/mbpt2l'
    #df1 = '/Users/yul20/atomic/juan/mbpt2m'
    if n == 11:
        iso = 'na'
    elif n == 12:
        iso = 'mg'
    elif n == 13:
        iso = 'al'
    elif n == 17:
        iso = 'cl'
    elif n == 18:
        iso = 'ar'
    else:
        iso = None
    if iso != None and dm != None and dm != 'None':
        fm = '%s/%s26_mrmp2_levels.txt'%(dm, iso)
    else:
        fm = None
    if (df0 == 'None'):
        df0 = None
    if (df1 == 'None'):
        df1 = None
    r0 = aflev(df0, df1, a, n)
    r0.e = r0.e - r0.e0
    if fm != None:
        r1 = MLEV(fm)
    else:
        r1 = None
    r0.match(r1)
    #fw = 'mf%d.txt'%n
    r0.write(fw)
    return r0
