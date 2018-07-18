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

from fac import *
from types import *
import time

# reinitialize FAC
def reinit(m = 0):
    Reinit(m)

    return

# the index of n-th shell in a complex
def get_index(n, complex):
    for i in range(len(complex)):
        if (complex[i][0] == n):
            return i
        elif (complex[i][0] > n):
            return -i-1

    return -i-2

def get_terms(n, nq):
    if (nq == 1):
        t = ['%d*1'%n]
    elif (nq == 2):
        if (n == 2):
            t = ['2s2', '2s1 2p1', '2p2']
        else:
            t = ['%d*2'%n]
    else:
        if (n == 2):
            t = ['2s2 2p%d'%(nq-2)]
            if (nq < 8):
                t.append('2s1 2p%d'%(nq-1))
                if (nq < 7):
                    t.append('2p%d'%nq)
            if (nq > 8):
                raise ValueError('nq <= 8 for n = 2')
        elif (n == 3):
            if (nq < 8):
                t = ['3s2 3p%d'%(nq-2)]
                t.append('3s1 3p%d'%(nq-1))
                t.append('3s1 3p%d 3d1'%(nq-2))
                if (nq == 3):
                    t.append('3s2 3d1')
                else:
                    t.append('3s2 3p%d 3d1'%(nq-3))
            else:
                if (nq != 8):
                    t = ['3s2 3p6 3d%d'%(nq-8)]
                else:
                    t = ['3s2 2p6']
                t.append('3s2 3p5 3d%d'%(nq-7))
                t.append('3s1 3p6 3d%d'%(nq-7))
        else:
            raise ValueError('only n <= 3 supports nq > 2')

    return t


# a configuration complex
class COMPLEX:
    def __init__(self, bname, nele = 0):
        self.complex = []
        self.terms = []
        self.nrec = 0
        self.nrec_ext = 0
        self.bname = bname
        self.name = []

        if (nele > 0):
            self.set_ground(nele)

        return

    def set_terms(self):
        t = []
        for c in self.complex:
            t.append(get_terms(c[0], c[1]))
        if (len(t) == 0):
            self.terms = ['']
            return
        a = t[0]
        for b in t[1:]:
            c = []
            for i in a:
                for j in b:
                    c.append(i+' '+j)
            a = c
        self.terms = a

        return

    def set_name(self, name = []):
        if (len(name) > 0):
            self.name = name
            return

        a = []
        for i in range(len(self.terms)):
            a.append('%s%d'%(self.bname, i))
        self.name = a

        return

    def set_ground(self, nelectrons):
        n = 1
        nele = nelectrons
        g = []
        while (nele > 0):
            nqm = 2*n*n
            if (nele < nqm):
                g.append((n, nele))
            else:
                g.append((n, nqm))
            nele = nele - nqm
            n = n+1

        self.complex = g
        self.set_terms()
        self.set_name()

        return

    def set_excited(self, i0, i1, cbase):
        base = cbase.complex
        if (i0 >= i1):
            return -1
        i = get_index(i0, base)
        if (i < 0):
            return -1
        j = get_index(i1, base)
        ex = base[:]
        if (j >= 0):
            nq1 = ex[j][1]
            if (nq1 == 2*i1*i1):
                return 0
            else:
                ex[j] = (i1, nq1+1)
        else:
            j = -j - 1
            ex.insert(j, (i1, 1))

        nq0 = ex[i][1]
        if (nq0 == 1):
            ex.remove(ex[i])
        else:
            ex[i] = (i0, nq0-1)

        self.complex = ex
        self.set_terms()
        self.set_name()

        return 0

    def set_ionized(self, i0, cbase):
        base = cbase.complex
        i = get_index(i0, base)
        if (i < 0):
            return -1
        ion = base[:]
        nq0 = ion[i][1]
        if (nq0 == 1):
            ion.remove(ion[i])
        else:
            ion[i] = (i0, nq0-1)

        self.complex = ion
        self.set_terms()
        self.set_name()

        return 0

    def set_recombined(self, i0, bname):
        self.set_name(bname)
        self.nrec = i0
        self.complex = (bname, i0)

        return 0


# a group of complexes
class CGROUP:
    def __init__(self, type):
        self.cgroup = []
        self.type = type

        return

    def add_complex(self, complex):
        self.cgroup.append(complex)

        return


# a generic atomic ion
class ATOM:
    def __init__(self, nele, asym=0, dir=''):
        if (nele <= 0):
            raise ValueError('NELE must > 0')

        self.process = {'ce': 1,
                        'tr': 1,
                        'rr': 1,
                        'ci': 1,
                        'ai': 1}

        self.ecorrections = 1
        self.nele = nele
        self.nele_max = [0, 2, 10, 28]
        self.nele_sim = [6,7,8,9]

        self.grd_complex = COMPLEX('grd.', nele)
        self.exc_complex = []
        self.ion_complex = CGROUP('ion')

        if (self.nele <= self.nele_max[1]):
            self.n_shells = 1
            self.nterms = [-1,-1,-1]
            self.nexc_max = [9, 9, 9]
            self.nfrozen = [9, 9, 9]
            self.nexc_rec = [10, 8, 6]
            self.nrec_max = [25, 16, 10]
            self.rec_pw_max = [10, 9, 6]
            self.nrec_ext = 45
            self.n_decay = [10, 3, -1]
        elif (self.nele <= self.nele_max[2]):
            self.n_shells = 2
            self.nterms = [-1,-1,-1,-1]
            if (self.nele > 5):
                self.nterms = [3, 2, 2, 2]
            self.nexc_max = [7, 7, 7, 7]
            self.nfrozen = [10, 7, 7, 7]
            self.nexc_rec = [10, 4, 4, 4]
            self.nrec_max = [25, 16, 12, 8]
            self.rec_pw_max = [10, 9, 9, 5]
            self.nrec_ext = 45
            self.n_decay = [10, 4, 4, -1]
        elif (self.nele <= self.nele_max[3]):
            self.n_shells = 3
            self.nterms = [-1,-1,-1,-1,-1]
            self.nexc_max = [7, 7, 7, 7, 7]
            self.nfrozen = [10, 7, 7, 10, 7]
            self.nexc_rec = [4, 4, 4, 4, 0]
            self.nrec_max = [25, 10, 10, 10, 8]
            self.rec_pw_max = [10, 9, 6, 5, 5]
            self.nrec_ext = 45
            self.n_decay = [10, 4, 4, 4, 4]
        else:
            raise ValueError('ion with NELE >= %d not '
                             ' supported'%(self.nele_max[3]))

        if (type(asym) == StringType):
            self.set_atom(asym, dir=dir)

        return


    def set_atom(self, asym, dir=''):
        self.asym = asym
        self.bfiles = {'en': dir+'%s%02db.en'%(asym, self.nele),
                       'tr': dir+'%s%02db.tr'%(asym, self.nele),
                       'ce': dir+'%s%02db.ce'%(asym, self.nele),
                       'rr': dir+'%s%02db.rr'%(asym, self.nele),
                       'ai': dir+'%s%02db.ai'%(asym, self.nele),
                       'ci': dir+'%s%02db.ci'%(asym, self.nele)}

        self.afiles = {'en': dir+'%s%02da.en'%(asym, self.nele),
                       'tr': dir+'%s%02da.tr'%(asym, self.nele),
                       'ce': dir+'%s%02da.ce'%(asym, self.nele),
                       'rr': dir+'%s%02da.rr'%(asym, self.nele),
                       'ai': dir+'%s%02da.ai'%(asym, self.nele),
                       'ci': dir+'%s%02da.ci'%(asym, self.nele)}

        return


    def add_excited(self, n0, n1, ibase):
        if (ibase == -1):
            base = self.grd_complex
        else:
            base = self.exc_complex[0].cgroup[ibase]

        n1.sort()

        cg = CGROUP('exc')
        k = len(self.exc_complex)
        bn = 'exc.%d'%(k)
        for i1 in n1:
            ex = COMPLEX('%s.%d.'%(bn,i1))
            if (ex.set_excited(n0, i1, base) == 0):
                if (self.nele in self.nele_sim and
                    k > self.n_shells):
                    ex.terms = ex.terms[:1]
                    ex.name = ex.name[:1]
                elif (len(ex.terms) == 3 and self.nele > 3 and n0 > 1):
                    if (k > 0):
                        ex.terms = ex.terms[:2]
                        ex.name = ex.name[:2]
                cg.add_complex(ex)

        self.exc_complex.append(cg)

        return


    def add_ionized(self, n0, ibase):
        if (ibase == -1):
            base = self.grd_complex
        else:
            base = self.exc_complex[0].cgroup[ibase]

        k = len(self.ion_complex.cgroup)
        ion = COMPLEX('%s.%d.'%('ion', k))
        ion.set_ionized(n0, base)
        if (len(ion.terms) == 3 and self.nele > 4 and n0 > 1):
            if (ibase == -1 and n0 != self.n_shells):
                ion.terms = ion.terms[:2]
                ion.name = ion.name[:2]
            elif (ibase > 0 and self.nele > 6):
                ion.terms = ion.terms[:2]
                ion.name = ion.name[:2]
        self.ion_complex.add_complex(ion)

        return


    def add_recombined(self, n0, ibase):
        base = self.ion_complex.cgroup[ibase]
        n0.sort()
        cg = self.exc_complex[ibase]
        if (ibase > self.n_shells and self.nele in self.nele_sim):
            bname = base.name[:1]
        else:
            bname = base.name
            if (len(bname) > 2 and ibase != 0):
                    bname = bname[:2]
        for i0 in n0:
            rec = COMPLEX('rec')
            if (rec.set_recombined(i0, bname) == 0):
                cg.add_complex(rec)

        return


    def run_tr_ce(self, a, b, c, d, tr = [-1], ce = 1):
        if (type(a) == StringType or type(a) == IntType):
            low = [a]
        else:
            low = a
        if (type(b) == StringType or type(b) == IntType):
            up = [b]
        else:
            up = b

        if (self.process['tr'] != 0):
            for m in tr:
                s = 'TR: %s -> %s %d'%(str(a), str(b), m)
                Print(s)
                TransitionTable(self.bfiles['tr'], low, up, m)

        if (ce != 0 and self.process['ce'] != 0):
            if (type(c) == StringType or type(c) == IntType):
                low = [c]
            else:
                low = c
            if (type(d) == StringType or type(d) == IntType):
                up = [d]
            else:
                up = d
            s = 'CE: %s -> %s'%(str(c), str(d))
            Print(s)
            CETable(self.bfiles['ce'], low, up)

        Reinit(radial = 1, excitation = 1)

        return


    def run_ci_rr(self, a, b, c, d, ci = 1, rr = 1):
        if (rr != 0 and self.process['rr'] != 0):
            if (type(c) == StringType or type(c) == IntType):
                low = [c]
            else:
                low = c
            if (type(d) == StringType or type(d) == IntType):
                up = [d]
            else:
                up = d
            s = 'RR: %s -> %s'%(c, d)
            Print(s)

            RRTable(self.bfiles['rr'], low, up)

            SetPEGrid(0)
            SetUsrPEGrid(0)
            SetRRTEGrid(0)

        if (ci != 0 and self.process['ci'] != 0):
            if (type(a) == StringType or type(a) == IntType):
                low = [a]
            else:
                low = a
            if (type(b) == StringType or type(b) == IntType):
                up = [b]
            else:
                up = b
            s = 'CI: %s -> %s'%(a, b)
            Print(s)
            CITable(self.bfiles['ci'], low, up)


        if (ci != 0 or rr != 0):
            Reinit(radial = 1, recombination = 1, ionization = 1)

        return


    def run_ai(self, a, b, k):
        if (self.process['ai'] == 0):
            return

        s = 'AI: %s -> %s %d'%(a, b, k)
        Print(s)
        if (type(a) == StringType or type(a) == IntType):
            low = [a]
        else:
            low = a
        if (type(b) == StringType or type(b) == IntType):
            up = [b]
        else:
            up = b

        AITable(self.bfiles['ai'], low, up)

        Reinit(radial = 1, recombination = 1)

        return


    def energy_corrections(self):
        e = []
        ec = []
        nmin = self.nexc_max[0] + 1
        if (self.nele == 2):
            iec = range(0, -7, -1)
            if (self.asym == 'C'):
                etable = [-1.25,-0.098,-0.537,-0.537,-0.537,-0.537, -1.059]
                e3 = -0.339
            elif (self.asym == 'N'):
                etable = [-1.299,-0.126,-0.568,-0.568,-0.568,-0.568,-1.118]
                e3 = -0.3345
            elif (self.asym == 'O'):
                etable = [-1.245,-0.1558,-0.517,-0.517,-0.517,-0.517,-1.083]
                e3 = -0.285
            elif (self.asym == 'F'):
                etable = [-1.329,-0.136,-0.581,-0.581,-0.581,-0.581,-1.188]
                e3 = -0.369
            elif (self.asym == 'Ne'):
                etable = [-1.308,-0.104,-0.552,-0.552,-0.552,-0.552,-1.123]
                e3 = -0.36
            elif (self.asym == 'Na'):
                etable = [-1.326,-0.119,-0.554,-0.554,-0.554,-0.554,-1.131]
                e3 = -0.365
            elif (self.asym == 'Mg'):
                etable = [-1.317,-0.089,-0.538,-0.538,-0.538,-0.538,-1.133]
                e3 = -0.372
            elif (self.asym == 'Al'):
                etable = [-1.330,-0.095,-0.547,-0.547,-0.547,-0.547,-1.145]
                e3 = -0.390
            elif (self.asym == 'Si'):
                etable = [-1.310,-0.077,-0.517,-0.517,-0.517,-0.517,-1.128]
                e3 = -0.35
            elif (self.asym == 'S'):
                etable = [-1.375,-0.089,-0.539,-0.539,-0.539,-0.539,-1.146]
                e3 = -0.403
            elif (self.asym == 'Ar'):
                etable = [-1.378,-0.135,-0.512,-0.512,-0.512,-0.512,-1.139]
                e3 = -0.393
            elif (self.asym == 'Ca'):
                etable = [-1.382,-0.095,-0.505,-0.505,-0.505,-0.505,-1.136]
                e3 = -0.365
            elif (self.asym == 'Fe'):
                etable = [-1.497,-0.041,-0.491,-0.491,-0.491,-0.491,-1.135]
                e3 = -0.495
            else:
                etable = [-1.25,-0.12,-0.54,-0.54,-0.54,-0.54,-1.05]
                e3 = -0.335

            for i in range(-16,-26,-1):
                iec.append(i)
                etable.append(e3)
            ec.append((iec, etable))
            nmin = 1
        elif (self.nele >= 4 and self.nele <= 10):
            ie0 = [46, 125, 236, 272, 226, 113, 37]
            etable = [[]]*7
            if (self.asym == 'C'):
                etable = [[0.0, 7.995, 8.00835],
                          [0.0, 6.49269, 6.495627, 6.502615, 12.69004,
                           17.03862, 17.04218, 17.04808, 18.08634, 22.62958],
                          [0.0, 0.007863, 5.33173, 5.33446, 5.33797,
                           9.290152, 9.290464, 11.9637, 13.71565, 13.72078,
                           17.60912, 18.65486, 18.65549, 20.9198, 20.92213]]
            elif (self.asym == 'N'):
                etable = [[0.0, 9.97617, 10.00824],
                          [0.0, 8.33288, 8.34071, 8.35856, 16.20398,
                           21.76362, 21.77264, 21.78811,
                           23.41845, 29.1845],
                          [0.0, 0.02162, 7.0903, 7.0977, 7.10776,
                           12.52537, 12.52620, 16.24247, 18.08629,
                           18.09994, 23.1599, 25.17805, 25.17982,
                           28.5665, 28.56704],
                          [0.0, 0.00604, 0.01622, 1.89897, 4.0529,
                           5.80055, 11.43596, 11.43758, 11.43777,
                           13.54126, 13.54114, 13.54199, 17.87703,
                           19.23327, 20.67631]]
            elif (self.asym == 'O'):
                etable = [[0.0, 11.94898, 12.015],
                          [0.0, 10.15958, 10.17645, 10.21448,
                           19.68841, 26.46599, 26.48529, 26.51862,
                           28.7298, 35.69634],
                          [0.0, 0.04785, 8.85741, 8.87356, 8.89655,
                           15.73810, 15.73982, 20.37884, 22.37678,
                           22.40695, 28.707, 31.63531, 31.63891,
                           35.83336, 35.83436],
                          [0.0, 0.014032, 0.03796, 2.513566, 5.354351,
                           7.479323, 14.88123, 14.88473, 14.88533,
                           17.653, 17.6531, 17.65455, 23.19175,
                           24.43577, 26.09395, 35.18173, 35.20872, 35.22044],
                          [0.0, 3.324086, 3.326568, 5.017396, 5.017642,
                           14.85792, 14.87816, 14.88838, 20.57995, 20.58095,
                           24.26501, 26.35828, 26.37916]]
            elif (self.asym == 'Ne'):
                etable = [[0.0,15.88882,16.0933],
                          [0.0,13.7939,13.8503,13.9735,26.6507,
                           35.8726,35.936,36.0454,39.6402,49.3705],
                          [0.0,0.1624,12.4308,12.4857,12.5654,22.1922,
                           22.1958,28.6221,30.9083,31.0099,39.8707,44.5674,
                           44.5674,50.3356,50.3445],
                          [0.0,0.05113,0.13763,3.75567,7.92428,10.9552,21.8006,
                           21.8094,21.8122,25.8082,25.8085,25.8127,
                           33.5457,34.6369,37.6679,51.1659,51.2633,51.3051],
                          [0.0,5.11244,5.11801,7.74091,7.74174,22.7957,
                           22.8722,22.9122,31.50198,31.50466,37.14914,
                           39.67863,39.76677],
                          [0.0,0.07971,0.11413,3.20385,6.9122,25.329,
                           25.4018,25.4421,35.8908],
                          [37,38,39], [0.0,0.09675,26.91048]]
            elif (self.asym == 'Mg'):
                etable = [[0.0,19.8393,20.332],
                          [0.0,17.4203,17.5600,17.8650,33.6849,
                           45.3604,45.5219,45.791,50.22,61.9466],
                          [0.0,0.4094,16.104,16.245,16.453,28.7983,
                           28.8024,36.9823,39.5164,39.7646,51.281,57.7450,
                           57.7541,65.0486,65.072],
                          [0.0,0.1373,0.3625,5.0769,10.5576,14.64,28.8701,
                           28.8830,28.8913,34.0838,34.0829,34.0891,
                           43.9401,44.8968,49.2407,67.2386,67.4962,67.6041,
                           71.449,81.636],
                          [0.0,6.8633,6.86535,10.40476,10.41820,30.7416,
                           30.9445,31.051,42.3717,42.3769,49.8196,
                           52.716,52.9580,80.8212,81.1443],
                          [0.0,0.22108,0.31266,4.4543,9.5814,35.11386,
                           35.31422,35.42379,49.2815,82.197],
                          [0.0,0.2762,38.6251]]
            elif (self.asym == 'Al'):
                etable = [[0.0,21.8236,22.5413],
                          [0.0,19.2359,19.4405,19.8907,37.256,
                           50.1608,50.4017,50.795,55.7597,68.6604],
                          [0.0,0.60,18.011,18.224,18.530,32.202,
                           32.206,41.250,43.9,44.250,57.148,64.481,
                           64.489,72.553,72.597],
                          [0.0,0.21,0.54,5.792,11.93,16.594,32.506,
                           32.517,32.524,38.324,38.325,38.326,
                           49.224,50.11,55.119,75.39,75.776,75.941,
                           80.256,91.561],
                          [0.0,7.7258,7.7328,11.7293,11.7623,34.74,
                           35.046,35.207,47.84,47.846,56.172,
                           59.226,59.594,91.221,91.71],
                          [0.0,0.3387,0.4747,5.1041,10.937,40.0472,
                           40.353,40.5199,55.966,93.686],
                          [0.0,0.4268,44.4875]]
            elif (self.asym == 'Si'):
                etable = [[0.0,23.8127,24.8264],
                          [0.0,21.0528,21.3431,21.9846,40.875,55.008,
                           55.3582,55.9126,61.3971,75.4764],
                          [0.0,0.86672,19.962,20.270,20.712,35.688,35.692,
                           45.585,48.358,48.853,63.148,71.344,71.346,
                           80.188,80.266],
                          [0.0,0.31554,0.7952,6.56198,13.3654,18.693,
                           36.2322,36.2401,36.2581,42.6652,42.6517,42.6599,
                           54.6030,55.4138,61.0939,83.6601,84.22,84.4431,
                           89.2069,101.6285],
                          [0.0,8.57575,8.60705,13.0615,13.128,38.766,39.210,
                           39.444,53.357,53.361,62.566,65.765,66.05,
                           101.6,102.4],
                          [0.0,0.49,0.69,5.77392,12.3167,45.027,45.4757,
                           45.7206,62.692,105.2697],
                          [0.0,0.63, 50.3992]]
            elif (self.asym == 'S'):
                etable = [[0.0,27.8178,29.6854],
                          [0.0,24.6953,25.3376,26.4312,48.3022,64.8731,
                           65.5624,66.5617,73.0824,89.4949],
                          [0.0,1.62857,24.0383,24.6326,25.4695,42.98,
                           43.0232,54.501,57.6223,58.449,75.6397,
                           85.538,85.608,95.966,96.1876],
                          [0.0,0.64571,1.53593,8.32508,16.4811,23.0922,
                           44.0238,44.057,44.0595,51.6997,51.6949,51.7534,
                           65.7336,66.358,73.458,100.6383,101.7329,
                           102.1190,107.6756,122.3397],
                          [0.0,10.22155,10.36445,15.7429,15.9697,
                           46.9228,47.7788,48.2154,64.5614,64.5789,
                           75.4796,78.9653,80.051,122.8993,124.2783],
                          [0.0,0.99,1.32,7.22753,15.21,55.1714,56.0403,
                           56.523,76.3833,128.8468],
                          [0.0,1.2504,62.4439]]
            elif (self.asym == 'Ar'):
                etable = [[0.0,31.8672,35.0383],
                          [0.0,28.352,29.2433,31.3287,56.0672,75.0060,
                           76.2657,77.8972,85.5021,104.2226],
                          [0.0,2.809,28.549,29.626,31.048,50.85,50.983,
                           63.780,67.602,68.771,89.13,100.4,100.7,
                           112.6,113.0],
                          [0.0,1.2224,2.709,10.5426,20.1026,26.7,52.41,52.477,
                           52.50,61.426,61.472,61.625,77.596,77.940,86.626,
                           118.4,120.3,120.9,127.24,144.7],
                          [0.0,11.7563,12.1688,18.4957,19.099,55.289,56.798,
                           57.534,76.110,76.182,88.694,92.462,94.452,
                           146.0,148.0],
                          [0.0,1.7919,2.270,8.9062,18.412,65.667,
                           67.197,68.064,90.555,153.17],
                          [0.0,2.2398,74.8973]]
            elif (self.asym == 'Ca'):
                etable = [[0.0,35.9625,41.0286],
                          [0.0,32.023,33.408,36.817,64.300,85.435,
                           87.617,90.068,98.955,119.91],
                          [0.0,4.527,33.226,35.02,37.29,59.440,59.743,
                           73.421,78.576,80.051,103.50,116.5,117.1,
                           130.51,131.675],
                          [0.0,2.177,4.4539,13.46,24.508,34.20,61.580,61.69,
                           62.020,72.125,72.255,72.613,90.369,90.465,
                           100.97,137.318,140.579,141.338,148.176,167.489],
                          [0.0,13.126,14.074,21.37,22.733,63.95,66.439,
                           67.582,88.116,88.33,102.29,106.40,109.80,
                           167.114,171.111],
                          [0.0,3.032,3.5817,10.9364,22.1452,76.6555,
                           79.1315,80.6028,105.4,178.577],
                          [0.0,3.7246,87.908]]
            elif (self.asym == 'Fe'):
                etable = [[0.0,48.6,64.561],
                          [0.0,43.168,47.006,58.493,93.340,118.5,
                           127.35,132.87,149.30,176.4],
                          [0.0,14.663,50.157,57.05,63.636,91.316,
                           94.18,105.81,121.28,123.02,155.68,
                           173.134,176.91,194.609,201.811],
                          [0.0,9.156,14.5499,30.321,46.10,60.374,96.308,
                           96.379,99.674,113.61,114.67,116.83,135.83,
                           139.70,156.3,204.11,215.19,215.79,225.31,253.94],
                          [0.0,17.1867,21.8373,32.2694,40.089,93.3266,
                           101.769,104.486,129.262,131.220,148.193,154.042,
                           166.144,242.330,255.680],
                          [0.0,9.3298,11.0893,20.935,40.3122,114.4238,
                           122.0922,127.7063,157.1624,264.6047],
                          [0.0,12.7182,132.0063]]
            elif (self.asym == 'Ni'):
                etable = [[0.0,52.963,74.970],
                          [0.0,46.889,51.914,68.12,105.1,129.97,
                           143.11,149.74,171.04,199.80],
                          [0.0,20.328,56.79,67.0,75.0,104.6,109.6,
                           118.4,141.72,140.05,176.66,196.1,201.64,
                           220.82,232.21],
                          [0.0,13.609,20.084,40.250,57.51,72.765,
                           110.8,100.0,116.7,132.03,133.71,137.0,
                           155.04,161.74,180.97,231.19,245.16,247.89,
                           257.96,291.13],
                          [0.0,19.531,25.959,37.50,49.60,105.1,116.0,
                           120.0,145.83,149.20,166.70,173.4,190.50,
                           271.58, 290.23],
                          [0.0,11.50,15.905,26.315,50.43,129.35,
                           139.6,147.92,178.09,297.99],
                          [0.0,17.8486,149.06]]

            i = self.nele - 4
            if (len(etable[i])):
                iec = range(ie0[i], ie0[i]+len(etable[i]))
                ec.append((iec, etable[i]))

        return (ec, nmin)


    def run_en(self):
        SetAtom(self.asym)

        g = self.grd_complex.name
        c = self.exc_complex[0].cgroup[0].name

        if (self.nele > 1):
            ConfigEnergy(0)
        OptimizeRadial(g)
        if (self.nele > 1):
            ConfigEnergy(1)

        if (self.ecorrections):
            (ec, nmin) = self.energy_corrections()
            for a in ec:
                CorrectEnergy(a[0], a[1], nmin)

        # ground and the first excited complex are interacting
        Print('Structure: ground complex')
        Structure(self.bfiles['en'], g+c)
        PrepAngular(g+c)

        Print('Structure: ionized complexes')
        c = self.ion_complex.cgroup
        for a in c:
            if (len(a.name) == 0):
                continue
            s = '    %s'%str(a.name)
            Print(s)
            Structure(self.bfiles['en'], a.name)
        for i in range(len(c)):
            PrepAngular(c[i].name)
            for j in range(i+1, len(c)):
                PrepAngular(c[i].name, c[j].name)
        PrepAngular(c[0].name, g)

        Print('Structure: excited complexes')
        for i in range(len(self.exc_complex)):
            c = self.exc_complex[i].cgroup
            SetRecPWOptions(self.rec_pw_max[i])
            SetRecSpectator(self.nexc_max[i]+1, self.nfrozen[i])
            for j in range(len(c)):
                a = c[j]
                if (len(a.name) == 0):
                    continue
                if (i == 0 and j == 0):
                    x = self.ion_complex.cgroup[0]
                    continue
                if (a.nrec > 0):
                    s = '    (%s, %d)'%(str(a.name), a.nrec)
                    Print(s)
                    RecStates(self.bfiles['en'], a.name, a.nrec)
                else:
                    s = '    %s'%str(a.name)
                    Print(s)
                    Structure(self.bfiles['en'], a.name)
                    if (i == 0):
                        x = self.ion_complex.cgroup[0]
                        PrepAngular(x.name, a.name)
                        PrepAngular(g, a.name)
        return


    def run_tr_ce_all(self):
        g = self.grd_complex.name
        self.run_tr_ce(g, g, g, g, tr=[-1,1,-2,2], ce=1)
        for i in range(self.n_shells):
            c = self.exc_complex[i].cgroup
            for j in range(len(c)):
                if (len(c[j].name) == 0):
                    continue
                a1 = g
                b = c[j].name
                nb = len(b)
                if (c[j].nrec > 0 and
                    self.nterms[i] > 0 and nb > self.nterms[i]):
                    nb = self.nterms[i]
                    b = b[:nb]
                if (c[j].nrec == 0):
                    a2 = g
                else:
                    b = (b, c[j].nrec)
                    a2 = g[:1]
                if (i == 0 and (self.nele <= 2 or self.nele == 10 or j == 0)):
                    tr = [-1, 1, -2, 2]
                else:
                    tr = [-1]
                nmax = c[j].nrec
                if (nmax == 0):
                    nmax = c[j].complex[-1][0]
                if (nmax <= self.nexc_rec[i]):
                    ce = 1
                else:
                    ce = 0
                if (ce != 0 or len(tr) != 0):
                    self.run_tr_ce(a1, b, a2, b, tr=tr, ce=ce)

        if (self.nele < 11):
            ce0 = 1
        else:
            ce0 = 0
        for i in range(len(self.exc_complex)):
            c = self.exc_complex[i].cgroup
            for j in range(len(c)):
                n2 = c[j].nrec
                b = c[j].name
                nb = len(b)
                if (nb == 0):
                    continue
                if (self.nterms[i] > 0 and nb > self.nterms[i]):
                    nb = self.nterms[i]
                    b = b[:nb]
                if (n2 == 0):
                    n2 = c[j].complex[-1][0]
                else:
                    b = (b, n2)
                if (self.nele <= 2):
                    if (i == 0 and n2 <= self.nexc_max[i]):
                        ce = ce0
                    else:
                        ce = 0
                else:
                    if (i == 0 and j == 0):
                        ce = ce0
                    else:
                        ce = 0
                for m in range(j+1):
                    a = c[m].name
                    na = len(a)
                    if (na == 0):
                        continue
                    if (self.nterms[i] > 0 and na > self.nterms[i]):
                        na = self.nterms[i]
                        a = a[:na]
                    nmax = c[m].nrec
                    if (nmax == 0):
                        nmax = c[m].complex[-1][0]
                    tr = []
                    if (i == 0 and nmax <= self.n_decay[0]):
                        if (j == 0 or self.nele == 10):
                            tr = [-1, 1, -2, 2]
                        else:
                            tr = [-1]
                    if (c[m].nrec > 0):
                        a = (a, c[m].nrec)
                    if (m == 0):
                        ce1 = ce
                    else:
                        ce1 = 0
                    if (ce1 != 0 or tr != []):
                        self.run_tr_ce(a, b, a[:1], b, tr=tr, ce=ce1)

        for i in range(len(self.exc_complex)):
            c = self.exc_complex[i].cgroup
            d = self.exc_complex[0].cgroup
            for j in range(len(c)):
                b = c[j].name
                nb = len(b)
                if (nb == 0):
                    continue
                if (self.nterms[i] > 0 and nb > self.nterms[i]):
                    nb = self.nterms[i]
                    b = b[:nb]
                n2 = c[j].nrec
                if (n2 == 0):
                    n2 = c[j].complex[-1][0]
                else:
                    b = (b, n2)
                if (i == 0 and n2 <= self.n_decay[i]):
                    continue
                for m in range(len(d)):
                    if (self.nele <= 2 and i == 1 and m == 0):
                        ce1 = 1
                    else:
                        ce1 = 0
                    a = d[m].name
                    na = len(a)
                    if (na == 0):
                        continue
                    if (self.nterms[0] > 0 and na > self.nterms[0]):
                        na = self.nterms[0]
                        a = a[:na]
                    n1 = d[m].nrec
                    if (n1 == 0):
                        n1 = d[m].complex[-1][0]
                    else:
                        a = (a, n1)
                    if (n1 != n2):
                        try:
                            nq = c[j].complex[-1][1]
                        except:
                            nq = 1
                        if (nq > 1):
                            n3 = n2
                        else:
                            try:
                                n3 = c[j].complex[-2][0]
                            except:
                                n3 = 0
                        if (n3 != n1 or n2 > self.n_decay[i]):
                            continue
                    self.run_tr_ce(a, b, a, b, tr=[-1], ce=ce1)

        return


    def run_ci_rr_all(self):
        g = self.grd_complex.name
        for i in range(self.n_shells):
            b = self.ion_complex.cgroup[i].name
            for c in g:
                for d in b:
                    self.run_ci_rr([c], [d], [c], [d], ci=1, rr=1)

        b = self.ion_complex.cgroup[0].name
        for i in range(len(self.exc_complex[0].cgroup)):
            a = self.exc_complex[0].cgroup[i].name
            na = len(a)
            if (na == 0):
                continue
            if (self.nterms[0] > 0 and na > self.nterms[0]):
                na = self.nterms[0]
                a = a[:na]
            n = self.exc_complex[0].cgroup[i].nrec
            for c in a:
                if (n > 0):
                    cp = ([c], n)
                else:
                    cp = [c]
                for d in b[0:na]:
                    self.run_ci_rr(cp, [d], cp, [d], ci=1, rr=1)

        return


    def run_ai_all(self):
        for i in range(len(self.exc_complex)):
            c = self.exc_complex[i].cgroup
            if (i == 0):
                if (self.nele == 2):
                    continue
                if (self.nele == self.nele_max[self.n_shells-1]+1):
                    continue
                p = 1
            elif (i == self.n_shells-1 and
                  self.nele > self.nele_max[self.n_shells-1]+1):
                p = len(self.ion_complex.cgroup)
            else:
                p = 2
            for k in range(p):
                if (self.nele == self.nele_max[self.n_shells-1]+1):
                    t = k
                else:
                    if (p == 2 and k == 1):
                        t = self.n_shells
                    else:
                        t = k
                b = self.ion_complex.cgroup[t].name
                if (len(b) == 0):
                    continue
                if (p > 2 and t > 0):
                    if (i == t):
                        continue
                    nmax = self.ion_complex.cgroup[t].complex[-1][0]
                else:
                    nmax = 0
                for j in range(len(c)):
                    a = c[j].name
                    na = len(a)
                    if (na == 0):
                        continue
                    if (self.nterms[i] > 0 and na > self.nterms[i]):
                        na = self.nterms[i]
                        a = a[:na]
                    n = c[j].nrec
                    if (n > 0):
                        if (nmax > 0):
                            continue
                        a = (a, n)
                    else:
                        if (nmax > 0):
                            nmax1 = c[j].complex[-1][0]
                            if (nmax1 != nmax):
                                continue

                    if (n > self.nrec_max[i] and
                        c[j].nrec_ext == 0):
                        continue
                    m = k*10+i
                    if (c[j].nrec_ext > 0):
                        m = n + (self.nrec_max[i]+1)*1000
                        m = -m
                    if (self.nele > 5 and t > 0):
                        b1 = b[:na]
                    else:
                        b1 = b
                    if (i == 0):
                        if (n > 0):
                            for c1 in a[0]:
                                for c2 in b1:
                                    self.run_ai(([c1], n), [c2], m)
                        else:
                            for c1 in a:
                                for c2 in b1:
                                    self.run_ai([c1], [c2], m)
                    else:
                        self.run_ai(a, b1, m)

        return


    def print_table(self, type = [], v = 0):
        if (v != 0):
            MemENTable(self.bfiles['en'])

        for t in type:
            PrintTable(self.bfiles[t], self.afiles[t], v)

        return


    def set_configs(self):
        n = range(1, self.n_shells+1)
        n.reverse()
        nexc = self.nexc_max
        nrec = self.nrec_max

        for i in n:
            j = self.n_shells - i
            self.add_excited(i, range(nexc[j]+1), -1)
            self.add_ionized(i, -1)
            r = range(nexc[j]+1, nrec[j]+1)
            if (self.nrec_ext > nrec[j] and len(r) > 0):
                r.append(self.nrec_ext)
            self.add_recombined(r, j)
            try:
                self.exc_complex[j].cgroup[-1].nrec_ext = nrec[j]+1
            except:
                pass

        j = self.n_shells
        if (self.nele == self.nele_max[self.n_shells-1]+1):
            i = self.n_shells-1
            mrec = nrec[j+1]
            r = range(nexc[j]+1, mrec+1)
            if (self.nrec_ext > mrec and len(r) > 0):
                r.append(self.nrec_ext)
        else:
            i = self.n_shells
            mrec = nrec[j]
            r = range(nexc[j]+1, mrec+1)
            if (self.nrec_ext > mrec and len(r) > 0):
                r.append(self.nrec_ext)
        if (i > 0):
            self.add_excited(i, range(j+1, nexc[j]+1), 0)
            self.add_ionized(i, 0)
            self.add_recombined(r, j)
            try:
                self.exc_complex[j].cgroup[-1].nrec_ext = mrec+1
            except:
                pass

        if (self.nele > self.nele_max[self.n_shells-1]+1):
            i = self.n_shells
            j = self.n_shells+1
            self.add_excited(i, range(j+1, nexc[j]+1), 1)
            self.add_ionized(i, 1)
            r = range(nexc[j]+1, nrec[j]+1)
            if (self.nrec_ext > nrec[j] and len(r) > 0):
                r.append(self.nrec_ext)
            self.add_recombined(r, j)
            try:
                self.exc_complex[j].cgroup[-1].nrec_ext = nrec[j]+1
            except:
                pass

        # set up configurations
        c = self.grd_complex.name
        t = self.grd_complex.terms
        s = 'adding complex: %s, %s'%(str(c), str(t))
        Print(s)
        for i in range(len(t)):
            Config(t[i], group = c[i])

        for cg in self.exc_complex:
            for i in range(len(cg.cgroup)):
                if (cg.cgroup[i].nrec > 0):
                    continue
                c = cg.cgroup[i].name
                t = cg.cgroup[i].terms
                if (len(c) == 0):
                    continue
                s = 'adding complex: %s, %s'%(str(c), str(t))
                Print(s)
                for j in range(len(t)):
                    Config(t[j], group = c[j])

        cg = self.ion_complex.cgroup
        for i in range(len(cg)):
            c = cg[i].name
            if (len(c) == 0):
                continue
            t = cg[i].terms
            s = 'adding complex: %s, %s'%(str(c), str(t))
            Print(s)
            for j in range(len(t)):
                Config(t[j], group = c[j])

        return


    def run_fac(self):
        start = time.time()
        start0 = start

        # solve the structure
        s = 'EN...'
        Print(s)
        self.run_en()
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print(s)

        # transition rates and collisional excitation
        start = stop
        s = 'TR and/or CE...'
        Print(s)
        self.run_tr_ce_all()
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print(s)

        # CI and/or RR
        start = stop
        s = 'CI and/or RR...'
        Print(s)
        self.run_ci_rr_all()
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print(s)

        start = stop
        s = 'AI...'
        Print(s)
        self.run_ai_all()
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print(s)

        s = 'Total Time: %10.3E s'%(stop-start0)
        print(s)

        return


def atomic_data(nele, asym, iprint=1, dir='', **kw):
    if (type(nele) == IntType):
        n = [nele]
        if (kw.has_key('ecorrections')):
            kw['ecorrections'] = [kw['ecorrections']]
    else:
        n = nele

    asym = asym.capitalize()
    for im in range(len(n)):
        m = n[im]
        atom = ATOM(m)
        if (kw.has_key('no_ce')):
            atom.process['ce'] = not kw['no_ce']
        if (kw.has_key('no_tr')):
            atom.process['tr'] = not kw['no_tr']
        if (kw.has_key('no_rr')):
            atom.process['rr'] = not kw['no_rr']
        if (kw.has_key('no_ci')):
            atom.process['ci'] = not kw['no_ci']
        if (kw.has_key('no_ai')):
            atom.process['ai'] = not kw['no_ai']
        if (kw.has_key('nexc_max')):
            atom.nexc_max = kw['nexc_max']
        if (kw.has_key('nterms')):
            atom.nterms = kw['nterms']
        if (kw.has_key('nfrozen')):
            atom.nfrozen = kw['nfrozen']
        if (kw.has_key('nrec_max')):
            atom.nrec_max = kw['nrec_max']
        if (kw.has_key('nrec_ext')):
            atom.nrec_ext = kw['nrec_ext']
        if (kw.has_key('nexc_rec')):
            atom.nexc_rec = kw['nexc_rec']
        if (kw.has_key('n_decay')):
            atom.n_decay = kw['n_decay']
        if (kw.has_key('rec_pw_max')):
            atom.rec_pw_max = kw['rec_pw_max']
        if (kw.has_key('ecorrections')):
            atom.ecorrections = kw['ecorrections'][im]

        atom.set_configs()
        s = 'NELE = %d'%m
        Print(s)
        s = 'ATOM = %s'%asym
        z = ATOMICSYMBOL.index(asym.capitalize())
        s = '%s, Z = %d'%(s, z)
        Print(s)
        atom.set_atom(asym, dir=dir)
        atom.run_fac()
        if (iprint >= 0):
            t = ['en']
            if (atom.process['tr']):
                t.append('tr')
            if (atom.process['ce']):
                t.append('ce')
            if (atom.process['rr']):
                t.append('rr')
            if (atom.process['ai']):
                t.append('ai')
            if (atom.process['ci']):
                t.append('ci')
            atom.print_table(t, iprint)

        Reinit(0)

    return
