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
import struct
from sys import version_info
from pfac import fac
from pfac import const
from pfac import util
import os
from multiprocessing import Pool, cpu_count

def e2v(e, m=0):
    if (m == 0):
        k = e/const.Hartree_eV;
        k = 2.0*k*(1.0 + 0.5*const.Alpha**2*k);
        k = const.Alpha**2*k;
        k = np.sqrt(k/(1.0+k));
        k /= const.Alpha;
        k *= const.RBohr*1e-8*const.Rate_AU*1E-10;
    else:
        m = const.Mp_keV*m/const.Me_keV;
        k = e/const.Hartree_eV;
        k = 2.0*m*k*(1.0+0.5*const.Alpha**2*k/m);
        k = const.Alpha**2*k;
        k = np.sqrt(k/(m*m+k));
        k /= const.Alpha;
        k *= const.RBohr*1e-8*const.Rate_AU*1E-10;  

    return k

def k2lj(k):
    l = abs(k)
    if k < 0:
        l = l-1
        j = 1
    else:
        j = -1
    return l,j

def lj2k(l, j):
    k = l
    if j > 0:
        k = -k-1
    return k

def cfgnr(nq):
    s = ''
    if len(nq) == 0:
        return s
    for n,k,q in nq:
        if type(q) == int:
            s += ' %d%s%d'%(n,fac.SPECSYMBOL[k],q)
        else:
            s += ' %d%s%g'%(n,fac.SPECSYMBOL[k],q)
    return s[1:]

def nlq(s):
    i = 0
    n = 0
    l = -1
    q = 0
    while True:
        if s[i] == '[':
            n = int(s[:i])
            i += 1
            k = i
            while (s[i] != ']'):
                i += 1
                if i >= len(s):
                    break
            try:
                l = int(s[k:i])
            except:
                l = -1
            if i < len(s)-1:
                q = int(s[i+1:])
            else:
                q = 1
            break
        elif not s[i].isdigit():
            n = int(s[:i])
            try:
                l = fac.SPECSYMBOL.index(s[i])
            except:
                l = -1
            if i < len(s)-1:
                q = int(s[i+1:])
            else:
                q = 1
            break
        i += 1
        if i >= len(s):
            break
    return n,l,q

def nlqs(s):
    r = np.array([nlq(x) for x in s.split(' ')])
    i = np.argsort(r[:,0]*100 + r[:,1])
    return r[i,:]

def voigt_fwhm(gw, lw):
    return 0.5346*lw + np.sqrt(0.2166*lw**2 + gw**2)

def doppler_fwhm(ti, z=1, m=0.0):
    if m <= 0:
        m = fac.ATOMICMASS[z]
    return np.sqrt(8*np.log(2)*ti/(m*const.AMU*const.Me_eV))

def voigt(alpha, x):
    v = x/1.41421
    a=np.zeros(8)
    b=np.zeros(8)
    c=np.zeros(8)
    
    a[1]=122.607931777104326
    a[2]=214.382388694706425
    a[3]=181.928533092181549
    a[4]=93.155580458134410
    a[5]=30.180142196210589
    a[6]=5.912626209773153
    a[7]=0.564189583562615
    
    b[1]=122.607931773875350
    b[2]=352.730625110963558
    b[3]=457.334478783897737
    b[4]=348.703917719495792
    b[5]=170.354001821091472
    b[6]=53.992906912940207
    b[7]=10.479857114260399
    
    c[1]=0.5641641
    c[2]=0.8718681
    c[3]=1.474395
    c[4]=-19.57862
    c[5]=802.4513
    c[6]=-4850.316
    c[7]=8031.468
    
    n = len(v)
    H = np.zeros(n)
    vb = 2.5
    if (alpha <= .001):
        w = np.where(abs(v) >= vb)[0]
        if (len(w) > 0):
            v2   = v[w]* v[w]
            v3   = 1.0
            fac1 = c[1]
            fac2 = c[1] * (v2 - 1.0)
            for i in range(1,8):
                v3     = v3 * v2
                fac1 = fac1 + c[i] / v3
                fac2 = fac2 + c[i] / v3 * (v2 - i)
                
            H[w] = np.exp(-v2)*(1. + fac2*alpha**2 * (1. - 2.*v2)) + fac1 * (alpha/v2);
        w = np.where(abs(v) < vb)
    else:
        w = np.arange(0,n)
        
    if (len(w) > 0):
        p1 = alpha
        vw = v[w]
        o1 = -vw
        p2 = (p1 * alpha + o1 * vw)
        o2 = (o1 * alpha - p1 * vw)
        p3 = (p2 * alpha + o2 * vw)
        o3 = (o2 * alpha - p2 * vw)
        p4 = (p3 * alpha + o3 * vw)
        o4 = (o3 * alpha - p3 * vw)
        p5 = (p4 * alpha + o4 * vw)
        o5 = (o4 * alpha - p4 * vw)
        p6 = (p5 * alpha + o5 * vw)
        o6 = (o5 * alpha - p5 * vw)
        p7 = (p6 * alpha + o6 * vw)
        o7 = (o6 * alpha - p6 * vw)

        q1 = a[1] + p1 * a[2] + p2 * a[3] + p3 * a[4] + p4 * a[5] + p5 * a[6] + p6 * a[7];
        r1 =        o1 * a[2] + o2 * a[3] + o3 * a[4] + o4 * a[5] + o5 * a[6] + o6 * a[7];
        q2 = b[1] + p1 * b[2] + p2 * b[3] + p3 * b[4] +  p4 * b[5] + p5 * b[6] + p6 * b[7] + p7;
        r2 =        o1 * b[2] + o2 * b[3] + o3 * b[4] + o4 * b[5] + o5 * b[6] + o6 * b[7] + o7;

        H[w] = (q1 * q2 + r1 * r2) / (q2 * q2 + r2 * r2);
        
    return H;

def convd(xd, yd, s0, gw=None, lw=None, x0=None, x1=None):
    dx = s0/3.0
    xd = np.array(xd)
    yd = np.array(yd)
    if x0 is None:
        x0 = xd.min()-20.0*s0
    if x1 is None:
        x1 = xd.max()+20.0*s0
    x = np.arange(x0, x1, dx)
    y = np.zeros(len(x))
    p0 = 1.0/np.sqrt(2*np.pi)
    if lw is None:
        for i in range(len(xd)):
            s = s0
            if not gw is None:
                s = s0 + gw[i]
            p = p0/s
            t = (x-xd[i])/s
            w = np.nonzero((t > -20.0) & (t < 20.0))[0]
            if (len(w) > 0):
                y += p*yd[i]*np.exp(-0.5*t*t)
    else:
        for i in range(len(xd)):
            s = s0
            if not gw is None:
                s = s0 + gw[i]
            p = p0/s
            t = (x-xd[i])/s
            a = lw[i]/(1.414*s)
            y += p*yd[i]*voigt(a, t)
            
    return (x, y)

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
    a = lines[0][4:-1].split('[')
    header['FAC'] = a[0]
    if len(a) > 1:
        a = a[1][:-1].split('.')
        header['nthreads'] = int(a[0])
        header['uta'] = int(a[1])
        header['utaci'] = int(a[2])
    else:
        header['nthreads'] = 0
        header['uta'] = 0
        header['utaci'] = 0        
    lines = lines[1:]
    header['Endian'], lines = _read_value(lines, int)
    header['TSess'], lines = _read_value(lines, int)
    header['Type'], lines = _read_value(lines, int)
    header['Verbose'], lines = _read_value(lines, int)
    key = lines[0].split('\t')[0]
    header[key], lines = _read_value(lines, float)
    a = key.split(' ')
    header['asym'] = a[0]
    header['Z'] = int(header[key])
    header['NBlocks'], lines = _read_value(lines, int)
    return header, lines


def read_lev(filename):
    """ read *a.lev / *a.en file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    if header['NBlocks'] == 0:
        return header,()
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
        block['VNL'] = np.zeros(nlev, dtype=int)
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
            block['VNL'][i] = int(line[32:37])
            block['2J'][i] = int(line[38:42])
            block['ncomplex'][i] = get_lcomplex(line)
            block['sname'][i] = get_lsname(line)
            block['name'][i] = get_lname(line)

        return (block, )

    return header, read_blocks(lines)


def read_en(filename):
    return read_lev(filename)

def read_enf(filename):
    """ read en file with B&E """
    with open(filename, 'r') as f:
        lines = f.readlines()
    header, lines = _get_header(lines)
    if header['NBlocks'] == 0:
        return header,()
    lines = lines[1:]

    def read_blocks(lines):
        block = {}
        block['NELE'], lines = _read_value(lines, int)
        nlev, lines = _read_value(lines, int)
        block['EFIELD'], lines = _read_value(lines, float)
        block['BFIELD'], lines = _read_value(lines, float)
        block['FANGLE'], lines = _read_value(lines, float)
        lines = lines[1:]
        block['ilev'] = np.zeros(nlev, dtype=int)
        block['energy'] = np.zeros(nlev, dtype=float)
        block['pbasis'] = np.zeros(nlev, dtype=int)
        block['mbasis'] = np.zeros(nlev, dtype=int)
        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks
            block['ilev'][i] = int(line[:6])
            block['energy'][i] = float(line[6:29])
            block['pbasis'][i] = int(line[29:36])
            block['mbasis'][i] = int(line[36:41])
        return (block, )

    return header, read_blocks(lines)

def read_trf(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    header, lines = _get_header(lines)
    if header['NBlocks'] == 0:
        return header,()
    lines = lines[1:]

    def read_blocks(lines):
        block = {}
        block['NELE'], lines = _read_value(lines, int)
        ntrans, lines = _read_value(lines, int)
        block['MULTIP'], lines = _read_value(lines, int)
        block['GAUGE'], lines = _read_value(lines, int)
        block['MODE'], lines = _read_value(lines, int)
        block['EFIELD'], lines = _read_value(lines, float)
        block['BFIELD'], lines = _read_value(lines, float)
        block['FANGLE'], lines = _read_value(lines, float)

        block['upper_index'] = np.zeros(ntrans, dtype=int)
        block['lower_index'] = np.zeros(ntrans, dtype=int)
        block['upper_pbasis'] = np.zeros(ntrans, dtype=int)
        block['lower_pbasis'] = np.zeros(ntrans, dtype=int)
        block['upper_mbasis'] = np.zeros(ntrans, dtype=int)
        block['lower_mbasis'] = np.zeros(ntrans, dtype=int)
        block['energy'] = np.zeros(ntrans, dtype=float)
        block['rate'] = np.zeros(ntrans, dtype=float)
        nm = 2*abs(block['MULTIP'])+1
        block['mrate'] = np.zeros((ntrans,nm), dtype=float)
        j = 0
        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks
            im = i%nm
            block['mrate'][j,im] = float(line[66:80])
            if im != nm-1:
                continue
            block['upper_index'][j] = int(line[:6])
            block['upper_pbasis'][j] = int(line[6:13])
            block['upper_mbasis'][j] = int(line[13:17])
            block['lower_index'][j] = int(line[17:24])
            block['lower_pbasis'][j] = int(line[24:31])
            block['lower_mbasis'][j] = int(line[31:35])
            block['energy'][j] = float(line[38:52])
            block['rate'][j] = float(line[94:108])
            j += 1
            
        return (block, )

    return header, read_blocks(lines)

def read_tr(filename):
    """ read *a.tr file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    if header['NBlocks'] == 0:
        return header,()
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
        if header['uta'] > 0:
            block['sdev'] = np.zeros(ntrans, dtype=float)
            block['rci'] = np.zeros(ntrans, dtype=float)

        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks
            if header['uta'] == 0:
                block['upper_index'][i] = int(line[:6])
                block['upper_2J'][i] = int(line[7:9])
                block['lower_index'][i] = int(line[10:16])
                block['lower_2J'][i] = int(line[17:19])
                block['Delta E'][i] = float(line[20:33])
                block['gf'][i] = float(line[34:47])
                block['rate'][i] = float(line[48:61])
                block['multipole'][i] = float(line[62:75])
            else:
                a = line.split()
                block['upper_index'][i] = int(a[0])
                block['upper_2J'][i] = int(a[1])
                block['lower_index'][i] = int(a[2])
                block['lower_2J'][i] = int(a[3])
                block['Delta E'][i] = float(a[4])
                block['sdev'][i] = float(a[5])
                block['gf'][i] = float(a[6])
                block['rate'][i] = float(a[7])
                block['multipole'][i] = float(a[8])
                block['rci'][i] = float(a[9])
                
        return (block, )

    return header, read_blocks(lines)

def read_ai(filename):
    """ read *a.ai file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    if header['NBlocks'] == 0:
        return header,()
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
            block['bound_index'][i] = int(line[:6])
            block['bound_2J'][i] = int(line[7:9])
            block['free_index'][i] = int(line[10:16])
            block['free_2J'][i] = int(line[17:19])
            block['Delta E'][i] = float(line[20:31])
            block['AI rate'][i] = float(line[32:43])
            block['DC strength'][i] = float(line[44:55])

        return (block, )

    return header, read_blocks(lines)


def read_ce(filename):
    """ read *a.ce file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    if header['NBlocks'] == 0:
        return header,()
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
            block['lower_index'][tr] = int(line[:6].strip())
            block['lower_2J'][tr] = int(line[7:9].strip())
            block['upper_index'][tr] = int(line[10:16].strip())
            block['upper_2J'][tr] = int(line[17:19].strip())
            block['Delta E'][tr] = float(line[20:30].strip())
            nsub = int(line[31:])
            if block['MSUB']:
                block['collision strength'][tr] = np.zeros(
                    (nusr, nsub), dtype=float)
                block['crosssection'][tr] = np.zeros(
                    (nusr, nsub), dtype=float)

            line = lines[0]
            lines = lines[1:]
            block['bethe'][tr] = float(line[:11].strip())
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
    if header['NBlocks'] == 0:
        return header,()
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
            block['bound_index'][tr] = int(line[:6])
            block['bound_2J'][tr] = int(line[7:9])
            block['free_index'][tr] = int(line[10:16])
            block['free_2J'][tr] = int(line[17:19])
            block['Delta E'][tr] = float(line[20:31])
            block['Delta L'][tr] = int(line[32:])
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
    if header['NBlocks'] == 0:
        return header,()
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
            block['bound_index'][tr] = int(line[:6])
            block['bound_2J'][tr] = int(line[7:9])
            block['free_index'][tr] = int(line[10:16])
            block['free_2J'][tr] = int(line[17:19])
            block['Delta E'][tr] = float(line[20:31])
            block['Delta L'][tr] = int(line[32:])
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
    if header['NBlocks'] == 0:
        return header,()
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

            block['block'][tr] = int(line[:6])
            block['level'][tr] = int(line[7:13])
            block['abs. energy'][tr] = float(line[14:27])
            block['population'][tr] = float(line[28:39])
            block['Delta E'][tr] = float(line[14:27])
            block['emissivity'][tr] = float(line[28:39])

        for i, line in enumerate(lines):
            if line.strip() == '':  # if empty
                blocks = read_blocks(lines[i+1:])
                return (block, ) + blocks
            
        return (block, )

    return header, read_blocks(lines)

def read_wfun(fn, npi=0, rmax=None):
    r = np.loadtxt(fn, unpack=1)
    d = np.loadtxt(fn, max_rows=24, comments='@',
                   usecols=0, delimiter='=', dtype=str)
    for i in range(24):
        if d[i][1:].strip() == 'ilast':
            break
    h = np.loadtxt(fn, max_rows=24, comments='@', delimiter='=', usecols=1)
    i = int(h[i])+1
    h = np.loadtxt(fn, comments='@', max_rows=3, usecols=3, unpack=1)
    n = int(h[0])
    k = int(h[1])
    e = float(h[2])
    r0 = r[1][:i]
    p0 = r[4][:i]
    q0 = r[5][:i]
    pa = ()
    if n == 0:
        pa = (r[1][i:],r[3][i:],r[2][i:])
        if npi == 0:
            re = r[1][i:]
            te = r[3][i:]
            ae0 = r[2][i:]
            ae1 = r[4][i:]
            ae2 = r[5][i:]
        else:
            a0 = list(r[2][i:])
            t0 = list(r[3][i:])
            a1 = list(r[4][i:])
            a2 = list(r[5][i:])
            r1 = list(r[1][i:])
            ts = []
            nr = len(r1)
            tnr = 0
            for i in range(1,nr):
                dt = t0[i]-t0[i-1]
                nt = int(2+dt/(np.pi/npi))
                ts.append(np.linspace(t0[i-1],t0[i],nt)[:-1])
                tnr += nt-1
                if not rmax is None:
                    if (r1[i] >= rmax):
                        break
            re = np.zeros(tnr)
            te = np.zeros(tnr)
            ae0 = np.zeros(tnr)
            ae1 = np.zeros(tnr)
            ae2 = np.zeros(tnr)
            i0 = 0
            for x in ts:
                x = list(x)
                i1 = i0 + len(x)
                te[i0:i1] = x
                re[i0:i1] = util.UVIP3P(t0, r1, x)
                ae0[i0:i1] = util.UVIP3P(t0, a0, x)
                ae1[i0:i1] = util.UVIP3P(t0, a1, x)
                ae2[i0:i1] = util.UVIP3P(t0, a2, x)
                i0 = i1
        r0 = np.append(r0, re)
        p0 = np.append(p0, ae0*np.sin(te))
        q0 = np.append(q0, ae1*np.cos(te)+ae2*np.sin(te))
    return n,k,e,r0,p0,q0,pa

def read_pot(fn, cfg=None, header=None):
    if cfg is None and header is None:
        return np.loadtxt(fn, unpack=1)
    nw = 0
    with open(fn, 'r') as f:
        rs = f.readlines(20000)
        for i in range(len(rs)):
            if len(rs[i]) > 6 and rs[i][:6] == '# Mean':
                nw = i
                break
    if not header is None:
        rs = np.loadtxt(fn, unpack=1, max_rows=nw, comments='@',
                        usecols=1, dtype=str)
        rd = np.loadtxt(fn, unpack=1, max_rows=nw, comments='@', usecols=3)
        r = {}
        for i in range(len(rs)):
            r[rs[i]] = rd[i]
        if len(header) == 0:
            return r
        return r[header]
    nc = np.int32(np.loadtxt(fn, unpack=1, skiprows=nw, max_rows=1,
                             comments='@', usecols=3))
    d = np.loadtxt(fn, unpack=1, comments='@', skiprows=nw+1,
                   max_rows=nc, usecols=range(1,10))
    if len(cfg) == 0:
        return d
    if cfg == 'bnd':
        ns = np.int32(d[1])
        ks = np.int32(d[2])
        nlq = []
        for n in range(min(ns),max(ns)+1):
            for k in range(n):
                w = np.where((ns==n)&((ks==k)|(ks==-(k+1))))[0]
                if len(w) == 0:
                    continue
                if np.mean(d[-1][w]) >= 0:
                    continue
                nlq.append((n,k,round(np.sum(d[3][w]),2)))
        return cfgnr(nlq)
    if (cfg[-1] == '+'):
        j = 1
    elif (cfg[-1] == '-'):
        j = -1
    else:
        j = 0
    if j == 0:
        k = fac.SPECSYMBOL.index(cfg[-1:])
        n = int(cfg[:-1])
    else:
        k = fac.SPECSYMBOL.index(cfg[-2:-1])
        n = int(cfg[:-2])
    dn = np.int32(d[1])
    dk = np.int32(d[2])
        
    w1 = np.where((dn == n) & (dk == k))[0]
    w2 = np.where((dn == n) & (dk == -(k+1)))[0]

    if j == 0:
        fb = np.sum(d[3][w1])+np.sum(d[3][w2])
        eb = -(np.sum(d[3][w1]*d[8][w1]) +
               np.sum(d[3][w2]*d[8][w2]))/fb * 27.21
    elif j == 1:
        fb = np.sum(d[3][w2])
        eb = np.sum(d[3][w2]*d[8][w2])/fb * 27.21
    else:
        fb = np.sum(d[3][w1])
        eb = np.sum(d[3][w1]*d[8][w1])/fb * 27.21
    return fb,eb

def valence_shells(fn, nr=0):
    z = read_pot(fn, header='Z')
    r = read_pot(fn, cfg='')
    i = r[-1].argmax()
    r[-1] = r[-1]
    z1 = 1+(z-r[4][-1])
    de1 = ((z1+1)**2-z1**2)/(2*r[1][-1]**2)
    z1 = 1+(z-r[4])
    de0 = ((z1+1)**2-(z1**2))/(2*r[1]**2)
    de = de0-de1
    ei = -r[-1][i]
    ex = r[-1][i]-r[-1]+de
    w = np.where(ex < ei)[0]
    x = []
    nm = int(r[1].max())
    f = np.zeros((nm,nm),dtype=np.int32)
    for i in w:
        if nr == 0:
            x.append((int(r[1][i]),int(r[2][i])))
        else:
            k,j = k2lj(int(r[2][i]))
            if f[int(r[1][i])-1,k] == 0:
                x.append((int(r[1][i]),k))
                f[int(r[1][i])-1,k] = 1
    return x
    
def read_rra(fn):
    """ read the output of the Asymmetry function """
    r = {}
    with open(fn, 'r') as f:
        lines = f.readlines()
    if len(lines) < 3:
        return r

    a = ''
    for i in range(len(lines)):
        x = lines[i]
        if len(x) < 2:
            continue
        x = x.split()
        if (x[0] == '#'):
            if a == '':
                a = x[1]
            else:
                if len(x) == 8 or len(x) == 9:
                    s = x[1]
                    nt = int(x[6])
                    nq = int(x[-1])
                else:
                    s = x[1]+'_'+x[2]
                    nt = int(x[4])
                    nq = int(x[-1])
                if len(x) == 6 or len(x) == 8:
                    if (nq > 1000):
                        nq = nq%1000
                if (1 == nq%2):
                    nq = nq+1
                nq = nq+1
                sk = a+':'+s
                a = ''
                p = 0
                d = np.zeros((nq*2+9,nt))
                r[sk] = d
        else:
            if p < nt:
                for j in range(len(x)):
                    d[j,p] = float(x[j])                
            else:
                k = p%nt
                q = int(p/nt)+8
                d[q,k] = float(x[1])
                d[q+nq,k] = float(x[2])                
            p = p+1
    return r

def interp_rra(d, ea, aa=None):    
    if aa is None:
        x = np.log(ea)
        x0 = np.log(d[0])
        if (len(x0) < 1):
            r = np.zeros((10,1))
            r[:9] = d.copy()
            r[9] = (1-r[6])/(1+r[6])
        else:
            r = np.zeros((10,len(x)))
            r[0] = ea.copy()
            r[1] = r[0]+d[1][0]-d[0][0]
            for i in range(2,9):
                r[i] = np.exp(np.interp(x, x0, np.log(d[i])))
            r[9] = (1-r[6])/(1+r[6])
        return r
    x0 = np.log(d[0])
    nq = int((d.shape[0]-8)/2)
    if len(x0) == 1:
        b = d[9:nq+9,0]
        bp = d[nq+9:,0]
    else:
        x = np.log(ea)
        b = np.zeros(nq)
        bp = np.zeros(nq)
        for q in range(nq):
            b[q] = np.interp(x, x0, d[9+q])
            bp[q] = np.interp(x, x0, d[9+nq+q])
    na =len(aa)
    r = np.zeros((5,na))
    r[0] = aa.copy()
    r[1] = np.cos(r[0]*np.pi/180)
    for i in range(na):
        x = r[1][i]
        x = max(-0.999999,x)
        x = min( 0.999999,x)
        b0 = 0.0
        b1 = 0.0
        for q in range(nq):
            if abs(x) > 0.99999:
                y0 = 1.0
            else:
                y0 = fac.Legendre(abs(x), q, 0)[0]
            if (q%2 == 1 and x < 0):
                y0 = -y0
            b0 = b0 + b[q]*y0
            if q >= 2:
                if abs(x) > 0.99999:
                    y1 = 0.0
                else:
                    y1 = fac.Legendre(abs(x), q, 2)[0]
                if (q%2 == 1 and x < 0):
                    y1 = -y1
                b1 = b1 + bp[q]*y1/(q*(q-1))
        r[2,i] = b0
        r[3,i] = b1
        r[4,i] = -b1/b0
    return r
        
def read_rt(filename):
    """ read *a.rt file. """
    with open(filename, 'r') as f:
        lines = f.readlines()

    # header
    header, lines = _get_header(lines)
    if header['NBlocks'] == 0:
        return header,()
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

            block['block'][tr] = int(line[:6])
            block['level'][tr] = int(line[7:11])
            block['NB'][tr] = float(line[12:24])
            block['TR'][tr] = float(line[25:36])
            block['CE'][tr] = float(line[37:48])
            block['RR'][tr] = float(line[49:60])
            block['AI'][tr] = float(line[61:72])
            block['CI'][tr] = float(line[73:84])
            block['ncomplex'][tr] = line[85:].strip()

        return (block, )

    return header, read_blocks(lines)


MAX_SYMMETRIES = 256

def read_ham(filename):
    """ Read hamiltonian """
    header = OrderedDict()
    hamiltonian = []
    with open(filename, 'rb') as f:
        header['ng0'] = struct.unpack('i', f.read(4))[0]
        header['ng'] = struct.unpack('i', f.read(4))[0]
        header['kg'] = [struct.unpack('i', f.read(4))[0] 
                        for i in range(header['ng'])]
        header['ngp'] = struct.unpack('i', f.read(4))[0]
        header['kgp'] = [struct.unpack('i', f.read(4))[0]
                         for i in range(header['ngp'])]
            
        for s in range(MAX_SYMMETRIES):
            h = OrderedDict()
            hamiltonian.append(h)  # store the reference first
            
            h['s'] = struct.unpack('i', f.read(4))[0]
            h['dim'] = struct.unpack('i', f.read(4))[0]
            if h['dim'] <= 0:
                continue
            
            h['orig_dim'] = struct.unpack('i', f.read(4))[0]
            h['n_basis'] = struct.unpack('i', f.read(4))[0]
            h['basis'] = [struct.unpack('i', f.read(4))[0]
                             for i in range(h['n_basis'])]
            # number of elements
            n = struct.unpack('i', f.read(4))[0]
            
            h['i'] = np.zeros(n, int)
            h['j'] = np.zeros(n, int)
            h['value'] = np.zeros(n, float)
            for i in range(n):
                h['i'][i] = struct.unpack('i', f.read(4))[0]
                h['j'][i] = struct.unpack('i', f.read(4))[0]
                h['value'][i] = struct.unpack('d', f.read(8))[0]

    return header, hamiltonian
    
def comments(fn, hint=2048):
    with open(fn, 'r') as f:
        r = f.readlines(hint)
        r = [x[:2] for x in r]
        for i in range(len(r)):
            if r[i][0] == ' ' or r[i][0].isdigit():
                break
        return r[:i]

def valid_lines(fn):
    with open(fn) as f:
        for i,line in enumerate(f):
            line = line.lstrip()
            if len(line) > 20 and not line[0].isalpha():
                yield line

def mix_header_lines(fn):
    mh = '============Mixing Coefficients==================='
    n = len(mh)
    nh = 0
    with open(fn) as f:
        for i,line in enumerate(f):
            nh += 1
            if line[:n] == mh:
                break
    return nh

def read_bst(fn):
    nh = mix_header_lines(fn)
    h0 = np.loadtxt(fn, skiprows=1, max_rows=nh-3, unpack=1,
                    usecols=range(7), dtype=int)
    h1 = np.loadtxt(fn, skiprows=1, max_rows=nh-3, unpack=1,
                    usecols=(7,8,9), dtype=str)
    return h0,h1

def read_mix(fn):
    nh = mix_header_lines(fn)
    r = np.loadtxt(fn, unpack=1, skiprows=nh)
    return r
    
def load_fac(fn):
    r = np.loadtxt(valid_lines(fn), unpack=1, ndmin=2)
    return r

def load_atbase(fn):
    ks = {'s':0, 'p':1, 'd':2, 'f':3, 'g':4, 'h':5}
    with open(fn) as f: 
        d = f.readlines()
        nx = len(d)
        ilev = np.zeros(nx, dtype=np.int32)
        e = np.zeros(nx)
        j = ilev.copy()
        wj = j.copy()
        k = j.copy()
        p = j.copy()
        s = np.chararray(nx, itemsize=128)
        tm = s.copy()
        cn = s.copy()
        iup = j.copy()
        ilo = j.copy()
        te = e.copy()
        tf = e.copy()
        tt = j.copy()
        ti = j.copy()
        t = 0
        q = -1
        for i in range(nx):
            x = d[i]
            if len(x) < 6:
                continue
            if x[2:6] == 'Tran':
                q = 0                
            if x[0] == '#':
                continue
            if len(x) < 90:
                continue
            if q >= 0:
                iup[q] = int(x[:6])
                ilo[q] = int(x[13:21])
                ti[q] = int(x[30:39])
                tt[q] = int(x[48:57])
                te[q] = float(x[57:74])
                tf[q] = float(x[74:])
                q += 1
                continue
            ilev[t] = int(x[:6])
            k[t] = int(x[6:15])
            e[t] = float(x[15:33])
            j[t] = int(2*float(x[53:63]))
            wj[t] = int(float(x[33:53]))
            tm[t] = x[72:76]
            a = x[82:]
            cn[t] = a
            a = a.replace('(','').replace(')','.')
            a = a.strip()
            if (len(a) > 0):
                a = a[:-1]
                b = a.split('.')
                kt = 0
                a = ''
                for c in b:
                    if c[0] == '{':
                        c = c[4:]                    
                    kt += ks[c[1:2]]*int(c[2:])
                    a += '.' + c
                p[t] = kt%2
                s[t] = a[1:]
                t += 1
        return ilev[:t],k[:t],e[:t],j[:t],p[:t],s[:t],iup[:q],ilo[:q],ti[:q],tt[:q],te[:q],tf[:q],tm[:t],wj[:t],cn[:t]

def remove_closed(s):
    b = s.split('.')
    a = []
    for x in b:
        if x.rfind(')') < 0:
            a.append(x)
    r = ''
    na = len(a)
    if na == 0:
        return r
    for i in range(na):
        n,l,q = nlq(a[i])
        if (i == na-1 and r == '') or (q > 0 and q < 2*(2*l+1)):
            r = r + '.%d%s%d'%(n,fac.SPECSYMBOL[l],q)
    return r[1:]

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
        self.ib = self.ib[i]
        self.ibk = self.ibk[i]

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
        self.ib = np.zeros(n, dtype=int)
        self.ibk = np.zeros(n, dtype=int)
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
        self.ib[:ng] = g.ib
        self.ib[ng:] = c.ib
        self.ibk[:ng] = g.ibk
        self.ibk[ng:] = c.ibk

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
        self.ib = np.zeros(n, dtype=int)
        self.ibk = np.zeros(n, dtype=int)
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
        self.ib[:ng] = g.ib[wg]
        self.ib[ng:] = c.ib[wc]
        self.ibk[:ng] = g.ibk[wg]
        self.ibk[ng:] = c.ibk[wc]
        
    def read_atbase(self, f, zi=18, ki=2):
        if type(f) == type(''):
            r = load_atbase(f)
        else:
            r = f
        w = np.where(r[1] == ki)[0]
        wi = np.where(r[1] == ki+1)[0]
        w = np.append(w,wi[0])
        nw = len(w)
        self.z = zi
        self.asym = fac.ATOMICSYMBOL[zi]
        self.ilev = r[0][w]
        self.ib = self.ilev.copy()
        self.ib[:] = -1
        self.nele = np.int32(zi-r[1][w]+1)
        self.e = r[2][w]
        self.e0 = self.e[0]
        self.j = r[3][w]
        self.p = r[4][w]
        self.s = np.array([x.decode() for x in r[5][w]])
        self.n = np.array([x.decode() for x in r[12][w]])
        self.wj = r[13][w]
        w = np.where(self.nele == self.nele[0]-1)[0]
        if len(w) > 0:
            self.ei = self.e[w[0]]-self.e0
        else:
            self.ei = 0.0
        
    def __init__(self, f, ig=0, zi=0, ki=0):
        if f is None:
            return
        if zi > 0:
            self.read_atbase(f, zi=zi, ki=ki)
            return
        (hlev,blev) = read_lev(f)
        b0 = blev[0]
        b0['NELE'] = np.repeat(b0['NELE'],len(b0['ILEV']))
        b0['ibk'] = np.repeat(0, len(b0['ILEV']))
        for i in range(1, len(blev)):
            b = blev[i]
            b['NELE'] = np.repeat(b['NELE'],len(b['ILEV']))
            b['ibk'] = np.repeat(i,len(b['ILEV']))
            for kn in b0.keys():
                b0[kn] = np.append(b0[kn], b[kn])

        self.nele = b0['NELE']
        if ki > 0:
            w = np.where(self.nele == ki)[0]
            for kn in b0.keys():
                b0[kn] = b0[kn][w]
            self.nele = self.nele[w]
        ks = {'s':0, 'p':1, 'd':2, 'f':3, 'g':4, 'h':5}
        if self.nele[0] <= 12:
            gc = [2, 2, 6, self.nele[0]-10, 0, 0]
        else:
            gc = [2, 2, 6, 2, self.nele[0]-12, 0]
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
        self.z = hlev['Z']
        self.asym = hlev['asym']
        self.e0 = hlev['E0']
        self.e = b0['ENERGY']
        self.ib = b0['IBASE']
        self.ibk = b0['ibk']
        self.p = b0['P']
        self.j = b0['2J']
        self.wj = self.j+1
        self.c = np.array([x.decode() for x in b0['ncomplex']])
        self.v = b0['VNL']
        self.s = np.array([x.decode() for x in b0['sname']])
        self.n = np.array([x.decode() for x in b0['name']])
        self.e = self.e + self.e0
        self.e0 = min(self.e)
        self.ig = np.zeros(len(self.e), dtype=int)
        self.ilev = None
        w = np.where(self.nele == self.nele[0]-1)[0]
        if len(w) > 0:
            self.ei = self.e[w[0]]-self.e0
        else:
            self.ei = 0.0
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
        self.im[:] = -1
        self.cm = np.chararray(len(self.e),itemsize=64)
        self.cm[:]=b'.'
        if m == None:
            return
        w = np.where(self.nele == self.nele[0]-1)[0]
        if len(w) > 0 and m.ei > 0:
            w0 = w[0]
            ei = self.e[w0]-self.e0            
            self.im[w] = -1
            self.em[w] = (self.e[w]-self.e0)+(m.ei-ei)
            self.cm[w] = b'.'
        cs = np.array([remove_closed(m.s[i]) for i in range(len(m.s))])
        cs0 = np.array([remove_closed(self.s[i]) for i in range(len(self.s))])
        uc = np.unique(cs)
        imd = np.zeros(len(m.s),dtype=np.int32)
        for c in uc:
            ns = len(c)
            for p in [0, 1]:
                jmin = max(min(self.j),min(m.j))
                jmax = min(max(self.j),max(m.j))
                for j in range(jmin,jmax+1,1):
                    #print([p,j,c])
                    w0 = np.where((self.p == p) &
                                  (self.j == j) &
                                  (cs0 == c))[0]
                    n0 = len(w0)
                    ew0 = self.e[w0]-self.e0
                    w1 = np.where((m.p == p) &
                                  ((m.j == j)|((m.j < 0)&(imd==0))) &
                                  (cs == c))[0]
                    ew1 = m.e[w1] - m.e0
                    n1 = len(w1)
                    if (n1 == 0):
                        continue
                    if (n0 == 0):
                        continue
                    i0 = 0
                    i1 = 0
                    while (i0 < n0 and i1 < n1):
                        wi0 = w0[i0]
                        wi1 = w1[i1]
                        dex = min(50.0, 5*self.nele[0])
                        if (m.j[wi1] < 0):
                            dex *= 0.05
                        dex = max(0.25, dex)
                        if abs(ew0[i0]-ew1[i1]) < dex:
                            if m.ilev is None:
                                self.im[wi0] = wi1
                            else:
                                self.im[wi0] = m.ilev[wi1]
                            self.em[wi0] = ew1[i1]
                            self.cm[wi0] = m.s[wi1]
                            imd[wi1] = 1
                            i0 += 1
                            i1 += 1
                        elif ew0[i0] < ew1[i1]:
                            i0 += 1
                        else:
                            i1 += 1
        self.cm = np.array([x.decode() for x in self.cm])
        
    def write(self, fn):
        f = open(fn, 'w')
        for i in range(len(self.e)):
            if self.ilev is None:
                ik = i
            else:
                ik = self.ilev[i]
            em = self.em[i]
            if i == 0 or em >= 0:
                de = self.e[i]-self.e0-self.em[i]
            else:
                em = 0.0
                de = 0.0
            s = '%4d %4d %15.8E %15.8E %10.3E %4d %4d %4d %-32s %-84s %-48s\n'%(ik,self.im[i],self.e[i]-self.e0,em,de,self.wj[i],self.j[i],self.p[i],self.s[i],self.n[i],self.cm[i])
            f.write(s)
        f.close()

def strnum(s):
    s = s.replace('"','').replace('[','').replace(' ','')
    if len(s) == 0:
        return 0.0        
    for i in range(len(s)):
        if not (s[i].isdigit() or s[i]=='.'):
            break
    return float(s[:i])

def valid_nistlev(fn):
    with open(fn) as f:
        for i,line in enumerate(f):
            line = line.replace('\,',' or ')
            yield line
            
class MLEV:
    def __init__(self, f, md=1):
        if f == None:
            return
        self.ilev = None
        if md == 0:
            self.e = np.transpose(np.loadtxt(f, usecols=3, dtype='float', skiprows=1, delimiter=' ; '))
            j = np.transpose(np.loadtxt(f, usecols=1, dtype='float', skiprows=1, delimiter=' ; ', converters={1:lambda x: eval('2*'+x)}))
            self.j = np.int32(j)
            pc = np.transpose(np.loadtxt(f, usecols=2, dtype='string', skiprows=1, delimiter=' ; '))
            self.p = np.int32(pc == 'o')
            self.c = np.transpose(np.loadtxt(f, usecols=0, dtype='string', skiprows=1, delimiter=' ; '))
            self.s = self.c
            self.ei = 0.0
            self.e0 = self.e[0]
        else:
            r = np.loadtxt(valid_nistlev(f), unpack=1, delimiter=',', dtype=str)
            r[1] = np.array([str(x).strip() for x in r[1]], dtype='<U128')
            w0 = np.where(r[1] == 'Limit')
            ri = r[:,w0[0]]
            w0 = np.where(r[1] != 'Limit')
            r = r[:,w0[0]]
            self.c = np.array([str(x).replace('?','').strip() for x in r[0]],dtype='<U128')
            self.t = np.array([str(x).replace('?','').strip() for x in r[1]])
            self.j = np.array([int(2*eval(x.replace('?','').split('or')[0])) for x in r[2]])
            self.wj = self.j+1
            self.p = np.array([int(len(x)>0 and x[-1]=='*') for x in self.t])        
            self.e = np.array([strnum(x) for x in r[4]])*const.Ryd_eV
            if len(ri[0]) > 0:
                self.ei = strnum(ri[4,0])*const.Ryd_eV
            else:
                self.ei = 0.0
            self.e0 = self.e[0]
            self.nele = np.zeros(len(self.c),dtype=np.int32)
            fs = f.split('/')[-1].split('-')
            self.z = fac.ATOMICSYMBOL.index(fs[0])
            self.nele[:] = 1+self.z-int(fs[1].split('.')[0])
            for i in range(len(self.c)):
                a = self.c[i].split(" ")
                a = a[0].split(".")
                tc = ''
                for b in a:
                    b = b.split('<')[0]
                    if b[0].isdigit():                        
                        if (not b[-1].isdigit()):
                            b += '1'
                        tc += '.'+b
                self.c[i] = tc[1:]
            self.s = self.c
            
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

def NISTCorr(ff, fn, fo):
    r0 = FLEV(ff)
    r1 = MLEV(fn, md=1)
    r0.match(r1)
    r0.write(fo)

def read_rp(f):
    d = {}
    if not os.path.exists(f):
        return d
    r = np.loadtxt(f, unpack=1, usecols=6, ndmin=1)
    s = np.loadtxt(f, unpack=1, usecols=5, dtype=str, ndmin=1)
    if (len(r) == 0):
        return d
    sn = [x[:-1] for x in s]
    for i in range(len(s)):
        d[s[i]] = r[i]
        if sn[i] in d:
            d[sn[i]] = 0.5*(d[sn[i]]+r[i])
        else:
            d[sn[i]] = r[i]
    return d

def read_rps(fs, ds):
    nd = len(ds)
    d = []
    for f in fs:
        d.append(read_rp(f))

    r = {}
    r['ds'] = np.array(ds)
    r['dx'] = r['ds']**(1/3.)
    for s in d[0].keys():
        r[s] = np.zeros(nd)
        for i in range(nd):
            if s in d[i]:
                r[s][i] = d[i][s]*const.Hartree_eV
            else:
                r[s][i] = 0.0
    return r

def read_rps_zk(z, k, i, odir):
    p = '%s/%s/k%02d/%s%02d'%(odir,fac.ATOMICSYMBOL[z],k,fac.ATOMICSYMBOL[z],k)
    t = np.loadtxt(p+'a.tdg', unpack=1, ndmin=2)
    nd = int(t[3][i])
    ds = np.exp(np.linspace(np.log(t[1][i]), np.log(t[2][i]), nd))
    fs = ['%sd%02dt%02db.rp'%(p,j,i) for j in range(nd)]    
    ds = np.append([0.0], ds)
    fs = [p+'b.rp'] + fs
    r = read_rps(fs, ds)
    r['z'] = z
    r['k'] = k
    r['zp'] = z-k+1
    r['it'] = i
    r['ts'] = t[0]
    if len(t) > 4:
        r['ts0'] = t[5]
        r['tzs'] = t[4][0]
    else:
        r['ts0'] = r['ts']
        r['tzs'] = 0.0
    return r

def interp_ipd_at(s, r0, d, eth=1e-3, ood=0):
    xd = r0['ds'][1:]
    yd = r0[s][1:] - r0[s][0]
    e0 = r0[s][0]
    w = np.where(r0[s][1:] < e0*eth)[0]
    n = len(w)
    if n == 0:
        return jsp(d, t0, zp, 1.0)
    xd = xd[w]
    yd = yd[w]
    x = np.log(d)
    xi = np.log(xd)
    yi = np.log(yd)
    i0 = r0['it']
    t0 = r0['ts'][i0]
    zp = r0['zp']
    dd = 1.025
    dx = np.log(dd)
    if x < xi[0]:
        if n == 1 or ood//10 == 0:
            y0 = jsp(d, t0, zp, 1.0)
            y1 = jsp(xd[0], t0, zp, 1.0)
            r = yd[0]*(y0/y1)
            if n > 1:
                y2 = jsp(d*dd, t0, zp, 1.0)
                k0 = (yi[1]-yi[0])/(xi[1]-xi[0])
                k1 = np.log(y2/y0)/dx
                r = r*(d/xd[0])**(k0-k1)
        else:
            k0 = (yi[1]-yi[0])/(xi[1]-xi[0])
            r = np.exp(yi[0] + k0*(x-xi[0]))
        return r
    if x > xi[-1]:
        if n == 1 or ood%10 == 0:
            y0 = jsp(d, t0, zp, 1.0)
            y1 = jsp(xd[-1], t0, zp, 1.0)
            r = yd[-1]*(y0/y1)
            if n > 1:
                y2 = jsp(d*dd, t0, zp, 1.0)
                k0 = (yi[-1]-yi[-2])/(xi[-1]-xi[-2])
                k1 = np.log(y2/y0)/dx
                r = r*(d/xd[-1])**(k0-k1)
        else:
            k0 = (yi[-1]-yi[-2])/(xi[-1]-xi[-2])
            r = np.exp(yi[-1] + k0*(x-xi[-1]))
        return r
    return np.exp(np.interp(x, xi, yi))

def interp_ipd(z, k, ds, ts, ss, odir=None, rs=None,
               eth=1e-3, rr=False, oot=1, ood=0):
    if (odir is None) and (rs is None):
        print('must have odir or rs')
        return
    if type(ss) == type(''):        
        ss = [ss]
    if type(ds) == type(0) or type(ds) == type(0.0):
        ds = [ds]
    if type(ts) == type(0) or type(ts) == type(0.0):
        ts = [ts]
    ns = len(ss)
    nd = len(ds)
    nti = len(ts)    
    ys = np.zeros((ns,nti,nd))
    if not (odir is None):
        p = '%s/%s/k%02d/%s%02d'%(odir,fac.ATOMICSYMBOL[z],k,
                                  fac.ATOMICSYMBOL[z],k)
        r = np.loadtxt(p+'a.tdg', unpack=1, ndmin=2)
        xt = r[0]
    else:
        xt = rs[k][0]['ts']
    nt = len(xt)
    for m in range(nti):
        t = ts[m]
        w = np.where(xt > t)[0]
        if len(w) == 0:
            i1 = nt-1
            if oot//10 > 0:
                i0 = i1-1
            else:
                i0 = i1
        else:
            i1 = w[0]
            if i1 == 0:
                if oot%10 > 0:
                    i0 = 0
                    i1 = 1
                else:
                    i0 = 0
            else:
                i0 = i1-1
        if rs is None:
            r0 = read_rps_zk(z, k, i0, odir)
            if i1 != i0:
                r1 = read_rps_zk(z, k, i1, odir)
            else:
                r1 = r0
        else:
            if odir is None:
                rp0 = rs[k][i0]
                rp1 = rs[k][i1]
            else:
                rp0 = rs[0]
                rp1 = rs[1]
            r0 = None
            r1 = None
            if i0 == rp0['it']:
                r0 = rp0
            elif i0 == rp1['it']:
                r0 = rp1
            if i1 == rp1['it']:
                r1 = rp1
            elif i1 == rp0['it']:
                r1 = rp0
            if r0 is None:        
                r0 = read_rps_zk(z, k, i0, odir)
            if i1 == i0:
                r1 = r0
            elif r1 is None:
                r1 = read_rps_zk(z, k, i1, odir)
        rsr = (r0,r1)
        zp = z-k+1
        for i in range(ns):
            s = ss[i]
            for j in range(nd):
                d = ds[j]
                y0 = interp_ipd_at(s, r0, d, eth=eth, ood=ood)
                if i1 != i0:
                    y1 = interp_ipd_at(s, r1, d, eth=eth, ood=ood)
                    dt = np.log(xt[i1]/xt[i0])
                    dt0 = np.log(t/xt[i0])
                    y = np.exp(np.log(y0)+dt0*np.log(y1/y0)/dt)
                elif i0 == 0:
                    j0 = jsp(d, t, zp, 1.0)
                    j1 = jsp(d, xt[0], zp, 1.0)
                    y = y0*(j0/j1)
                else:
                    j0 = jsp(d, t, zp, 1.0)
                    j1 = jsp(d, xt[-1], zp, 1.0)
                    y = y0*(j0/j1)
                ys[i,m,j] = y

    ys = np.squeeze(ys)
    if (len(ys.shape) == 0):
        ys = float(ys)
    if rr:
        return ys,rsr
    return ys

def read_ipd(fn):
    with open(fn, 'rb') as f:
        d = f.read()
    i0 = 0
    i1 = 28
    zk,ndm,nt,t0,dt = struct.unpack('=iiidd', d[i0:i1])
    tzs = (dt//100)/1e4
    dt = (dt - 100*(dt//100))
    t1 = t0 + (nt-1)*dt
    ts0 = np.exp(np.linspace(t0, t1, nt))
    z = zk%256
    kk = zk//256
    k1 = kk%256
    k0 = kk//256
    if k0 == 0:
        k0 = 1
    if k1 == 0:
        k1 = z
    rs = []
    for k in range(z+1):
        rs.append([])
    for k in range(k0, k1+1):
        rk = []
        for p in range(nt):
            i0 = i1
            i1 = i1+4        
            ns = struct.unpack('i', d[i0:i1])[0]
            if ns == 0:
                continue
            r = {}
            r['z'] = z
            r['k'] = k
            r['zp'] = z-k+1
            r['it'] = p
            r['tzs'] = tzs
            r['ts0'] = ts0
            r['ts'] = ts0
            if tzs > 0:
                r['ts'] = ts0*(z-k+1)**tzs
            for j in range(ns):
                i0 = i1
                i1 = i1 + 3*4 + 3*8
                n,m,nd,d0,dd,e0 = struct.unpack('=iiiddd', d[i0:i1])
                if nd == 0:
                    continue
                ipd = np.zeros(nd)
                for t in range(nd):
                    i0 = i1
                    i1 = i1 + 8
                    ipd[t] = struct.unpack('d', d[i0:i1])[0]
                if j == 0:
                    d1 = d0 + (nd-1)*dd
                    ds = np.exp(np.linspace(d0, d1, nd))
                    ds = np.append([0.0], ds)
                    r['ds'] = ds
                    r['dx'] = ds**(1/3.)
                es = np.append([-e0], ipd-e0)
                ss = '%d%s'%(n,fac.SPECSYMBOL[m])
                r[ss] = es
            rk.append(r)
        rs[k] = rk
    return rs

def tab_ipd(z, odir, wdir='', k0=0, k1=0, md=0, ifill=True, eth=1e-3):
    if k0 == 0:
        k0 = 1
    if k1 == 0:
        k1 = z
    a = fac.ATOMICSYMBOL[z]
    if wdir == '':
        fn = '%s/%s/%s%02d%02db.es'%(odir,a,a,k0,k1)
    else:
        fn = '%s/%s%02d%02db.es'%(wdir,a,k0,k1)
    if md == 0:
        f = open(fn, 'wb')
    else:
        f = open(fn, 'w')

    print('header: z=%d k=%d'%(z,k0))
    r = read_rps_zk(z, k0, 0, odir)
    ts = r['ts0']
    tzs = r['tzs']
    nt = len(ts)
    t0 = np.log(ts[0])
    t1 = np.log(ts[-1])
    dt = np.log(ts[1]/ts[0])
    if (tzs > 0):
        dt = dt + 1e6*tzs
    ndm = len(r['ds'])-1
    if md == 0:
        zk = z
        if not ifill:
            if k1 < z:
                zk = 256*k1 + zk
            if k0 > 1:
                zk = 256*256*k0 + zk
        f.write(struct.pack('iii', zk, ndm, nt))
        f.write(struct.pack('dd', t0, dt))
    else:
        s = '%3d %3d %3d %15.8E %15.8E %15.8E\n'%(z, ndm, nt, t0, t1, dt)
        f.write(s)
    
    for k in range(1, z+1):
        if k < k0 or k > k1:
            if ifill:
                for p in range(nt):
                    if md == 0:
                        f.write(struct.pack('i', 0))
                    else:
                        s = '%3d %15.8E %3d\n'%(p, ts[p], 0)
                        f.write(s)            
            continue        
        print('z=%d k=%d'%(z, k))
        for p in range(nt):
            r = read_rps_zk(z, k, p, odir)
            ds = r['ds'][1:]
            xd = np.log(ds)
            nn = []
            mm = []
            kk = []
            for ks in r.keys():
                try:
                    m = fac.SPECSYMBOL.index(ks[-1])
                    n = int(ks[:-1])
                    nn.append(n)
                    mm.append(m)
                    kk.append(ks)
                except:
                    m = -1
                    n = 0
            if md == 0:
                f.write(struct.pack('i', len(nn)))
            else:
                s = '%3d %15.8E %3d\n'%(p, ts[p], len(nn))
                f.write(s)
            nm = np.array(nn)*100+np.array(mm)
            ik = np.argsort(nm)
            for j in range(len(ik)):
                i = ik[j]
                ks = kk[i]
                n = nn[i]
                m = mm[i]
                y = r[ks][1:]
                e0 = r[ks][0]
                w = np.where(y < e0*eth)[0]
                if len(w) == 0:
                    i0 = 0
                    i1 = 0
                    nd = 0
                    d0 = 0
                    dd = 0
                else:
                    i0 = w[0]
                    i1 = w[-1]+1
                    nd = i1-i0
                    if (nd > len(w)):
                        y[i0:i1] = np.interp(xd[i0:i1], xd[w], y[w])
                    d0 = xd[i0]
                    d1 = xd[i1-1]
                    if nd <= 1:
                        dd = 0.0
                    else:
                        dd = xd[i0+1]-xd[i0]
                if md == 0:
                    f.write(struct.pack('iii', n, m, nd))
                    f.write(struct.pack('ddd', d0, dd, -e0))
                else:
                    s = '%3d %3d %2d %2d %3s %3d %15.8E %3d %15.8E %15.8E %15.8E %15.8E\n'%(z, k, n, m, ks, p, ts[p], nd, d0, d1, dd, -e0)
                    f.write(s)
                if nd > 0:
                    for i in range(i0,i1):
                        a = y[i]-e0
                        if md == 0:
                            f.write(struct.pack('d', a))
                        else:
                            if i < i1-1:
                                f.write('%15.8E '%a)
                            else:
                                f.write('%15.8E\n'%a)
    f.close()
                    
def CorrCLow(r, d, eden, md=0):
    de = np.zeros(len(r.e))
    for j in range(len(r.e)):
        k = r.nele[j]
        c = r.c[j]
        s = r.s[j]
        n = r.n[j]
        nc, ns, nn = fac.FillClosedShell(k, c, s, n)
        a = nn.split('.')
        for b in a:
            i = b.rfind('(')
            b = b[:i]
            i = b.rfind('-')
            if i < 0:
                i = b.rfind('+')
            q = int(b[i+1])
            if md == 0:
                s = b[:i+1]
            else:
                s = b[:i]
            if s in d:
                e0 = d[s][0]
                xi = np.log(eden)
                y = d[s]
                w = np.where((d['ds'] > 0) & (d[s] < 0.0))[0]
                if len(w) == 0:
                    de[j] = 1e31
                    break
                x = np.log(d['ds'][w])
                y = (d[s][w]-d[s][0])
                yi = np.interp(xi, x, y)
                if yi >= -d[s][0]:
                    de[j] = 1e31
                    break
                de[j] = de[j] + q*yi
            else:
                de[j] = 1e31
                break
    return de

def jsp(d, t, zp, zs):
    zs1 = zs+1.0
    ta = t/const.Hartree_eV
    da = d*const.RBohr**3
    di2 = 4*np.pi*zs1*da/ta
    di = np.sqrt(di2)
    k = zp*di/ta
    j = ((3*zs1*k+1)**(2/3)-1.)/(2*zs1)
    
    return j*t

def ipd_more(d, t, zp, zs):
    di = d/zs
    ra = 0.62*(1/(di*1e24))**(1./3)
    rd = 743.4*np.sqrt(t/(d*1e24)/(1+zs))
    r = rd/ra
    j = zp*(((1+r**3))**(2./3)-r**2)*2.16e-7/ra
    
    return j

def ipd_ek(d, t, zp, zs):
    j = 2.3e-7*(zp+1.)*((1+1./zs)*d*1e24)**(1/3.)

    return j

def jsp_zd(de, t, zp, zs):
    ta = t/const.Hartree_eV
    ea = de/const.Hartree_eV
    a = 1.5*zp/ea
    d0 = zp/(4*np.pi/3*a**3)
    d1 = zp/ea
    d1 = ta/(4*np.pi*(zs+1.)*d1*d1)
    u = 1/const.RBohr**3
    d0 = d0*u
    d1 = d1*u
    dmin = min(d0, d1)
    while (jsp(dmin, t, zp, zs) > de):
        dmin = dmin*2.0
    dmax = max(d0, d1)
    while (jsp(dmax, t, zp, zs) < de):
        dmax = dmax*2.0

    d = 0.5*(dmin+dmax)
    while (dmax-dmin > 1e-5*d):
        r = jsp(d, t, zp, zs)
        if r > de:
            dmax = d
        elif r < de:
            dmin = d
        else:
            break
        d = 0.5*(dmin+dmax)
    return d


def jj2lsj(pref, dd='.', ex='jj2lsj'):
    i = 0;
    fb = '%s_%03d'%(pref,i)
    ofn = '%s.LS'%(pref)
    r = os.system('rm -rf %s/%s'%(dd,ofn))
    while (True):
        if (not os.path.exists('%s/%s.c'%(dd,fb))):
            break
        dt = datetime.datetime.now()
        print('%s %s %s @ %s'%(dd,fb,ofn,dt))
        with open('%s/jj2lsj.in'%(dd),'w') as f:
            f.write('%s\n'%fb)
            f.write('Y\n')
            f.write('Y\n')
            f.write('Y\n')
        r = os.system('cd %s; %s < jj2lsj.in > jj2lsj.out 2>&1'%(dd,ex))
        if (os.path.exists('%s/%s.lsj.lbl'%(dd,fb))):
            if (not os.path.exists('%s/%s'%(dd,ofn))):
                r = os.system('cd %s; cat %s.lsj.lbl > %s'%(dd,fb,ofn))
            else:
                r = os.system('cd %s; cat %s.lsj.lbl >> %s'%(dd,fb,ofn))
        i = i + 1        
        fb = '%s_%03d'%(pref,i)
        ofn = '%s.LS'%(pref)
