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

from pfac.crm import *
from pfac.table import *
from pfac import const
from math import *
import sys
import time
import copy
import string
import pickle
import pprint

def tabulate_trates(dfile, neles, z=26, pref='Fe'):
    tbl = TABLE(fname=dfile,
                title='Total Ionization and Recombination Rate Coefficients',
                authors=['M. F. Gu'],
                date=time.localtime())
    d = 'Num. of Electrons'
    tbl.add_column(label='NELE', unit='None',
                   description=d, format='I2')
    d = 'Temperature'
    tbl.add_column(label='Temp', unit='[K]',
                   description=d, format='F4.2')
    d = 'Total DR rate coefficients'
    tbl.add_column(label='DR', unit='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total DR Arnaud & Raymond'
    tbl.add_column(label='DR_AR', unit='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total RR rate coefficients'
    tbl.add_column(label='RR', unit='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total RR Arnaud & Raymond'
    tbl.add_column(label='RR_AR', unit='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total DCI rate coefficients'
    tbl.add_column(label='CI', unit='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total DCI Arnaud & Raymond'
    tbl.add_column(label='DCI_AR', unit='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total EA rate coefficients'
    tbl.add_column(label='EA', unit='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    d = 'Total EA Arnaud & Raymond'
    tbl.add_column(label='EA_AR', unit='10^-10^cm^3^/s',
                   description=d, format='E8.2')
    tbl.open('w')
    tbl.write_header()
    for k in neles:
        dir0 = '%s%02d/'%(pref, k)
        with open(dir0+'rates2.sav', 'r') as f:
            rates2 = pickle.load(f)
        with open(dir0+'rates3.sav', 'r') as f:
            rates3 = pickle.load(f)
        logt = rates2['logt']
        nt = len(logt)
        if (k == neles[0]):
            tdr = rates2['tdr'][1][1]
            trr = rates2['trr'][1][1]
            for i in range(nt):
                b = (10**logt[i])*const.kb
                (a0, a1, a2) = Recomb(z, k-1, b)
                (b0, b1, b2) = Ionis(z, k-1, b)
                tbl.write_row(k-1, logt[i],
                              tdr[i], a2,
                              trr[i], a1,
                              0.0, b2,
                              0.0, b1)
        tdr = rates3['tdr'][2][1]
        trr = rates3['trr'][2][1]
        tci = rates2['tci'][0][1]
        tea = rates2['tea'][0][1]     
        for i in range(nt):
            b = (10**logt[i])*const.kb
            (a0, a1, a2) = Recomb(z, k, b)
            (b0, b1, b2) = Ionis(z, k, b)
            tbl.write_row(k, logt[i],
                          tdr[i], a2,
                          trr[i], a1,
                          tci[i], b2,
                          tea[i], b1)
        if (k == neles[-1]):
            tci = rates3['tci'][1][1]
            tea = rates3['tea'][1][1]
            for i in range(nt):
                b = (10**logt[i])*const.kb
                (a0, a1, a2) = Recomb(z, k+1, b)
                (b0, b1, b2) = Ionis(z, k+1, b)
                tbl.write_row(k+1, logt[i],
                              0.0, a2,
                              0.0, a1,
                              tci[i], b2,
                              tea[i], b1)
                
    tbl.close()    
    
def tabulate_rates(dfile, neles, z=26, pref='Fe'):
    tbl = TABLE(fname=dfile,
                title='Line Formation Rate Coefficients',
                authors=['M. F. Gu'],
                date=time.localtime())
    d = 'Num. of electrons'
    tbl.add_column(label='NELE',
                   unit='None',
                   description=d,
                   format='I2')
    d = 'Level Index'
    tbl.add_column(label='Index',
                   unit='None',
                   description=d,
                   format='I3')
    d = 'Temperature'
    tbl.add_column(label='Temp',
                   unit='[K]',
                   description=d,
                   format='F4.2')
    unit = 's^-1^'
    d = 'Total Depletion Rate'
    tbl.add_column(label = 'RT',
                   unit=unit,
                   description = d,
                   format = 'E8.2')
    unit = '10^-10^cm^3^/s'
    d = 'Collisional Excitation' 
    tbl.add_column(label='CE',
                   unit=unit,
                   description=d,
                   format='E8.2')
    d = 'Resonance Excitation'
    tbl.add_column(label='RE',
                   unit=unit,
                   description=d,
                   format='E8.2')
    d = 'Radiative Recombination'
    tbl.add_column(label='RR',
                   unit=unit,
                   description=d,
                   format='E8.2')
    d = 'CE + n=3 Cascades'
    tbl.add_column(label='CE+CS3',
                   unit=unit,
                   description=d,
                   format='E8.2')
    d = 'CE + All Cascades'
    tbl.add_column(label='CE+CS',
                   unit=unit,
                   description=d,
                   format='E8.2')
    d = 'RE + All Cascades'
    tbl.add_column(label='RE+CS',
                   unit=unit,
                   description=d,
                   format='E8.2')
    d = 'DR + RR + n=3 Cascades'
    tbl.add_column(label='DR+RR+CS3',
                   unit=unit,
                   description=d,
                   format='E8.2')
    d = 'DR + RR + All Cascades'
    tbl.add_column(label='DR+RR+CS',
                   unit=unit,
                   description=d,
                   format='E8.2')
    d = 'Collisional Ionization'
    tbl.add_column(label='CI',
                   unit=unit,
                   description=d,
                   format='E8.2')

    tbl.open('w')
    tbl.write_header()
    for k in neles:
        dir0 = '%s%02d/'%(pref, k)
        with open(dir0+'rates3.sav', 'r') as f:
            rates3 = pickle.load(f)
        with open(dir0+'rates2.sav', 'r') as f:
            rates2 = pickle.load(f)
        with open(dir0+'rates1.sav', 'r') as f:
            rates1 = pickle.load(f)
        logt = rates2['logt']
        nt = len(logt)
        rt1 = rates1['rt']
        ce1 = rates1['ce']
        cs1 = rates1['cs']
        cs2 = rates2['cs']
        rr2 = rates2['rr']
        cs3 = rates3['cs']
        ci3 = rates3['ci']
        re3 = rates3['re']
        a1 = rates2['abund']
        b = []
        for i in range(len(rt1)):
            b.append([0.0]*nt)
        c = []
        for i in range(10):
            c.append(copy.deepcopy(b))
        for i in range(nt):
            p1 = a1[i][k-1]
            p2 = a1[i][k]
            p = p1/p2
            for m in range(len(rt1)):
                c[0][m][i] = rt1[m][2][i]
            for m in range(len(ce1)):
                n = ce1[m][1]
                c[1][n][i] = ce1[m][2][i]
            for m in range(len(cs1)):
                n = cs1[m][1]
                c[5][n][i] = cs1[m][2][i]
                c[4][n][i] = cs1[m][3][i]
            for m in range(len(cs2)):
                if (cs2[m][0] != k):
                    continue
                n = cs2[m][1]
                c[8][n][i] = (cs2[m][2][i]-c[5][n][i])/p
                c[7][n][i] = (cs2[m][3][i]-c[4][n][i])/p
                if (c[7][n][i] < 0):
                    c[7][n][i] = 0.0
                if (c[8][n][i] < c[7][n][i]):
                    c[8][n][i] = c[7][n][i]
            for m in range(len(rr2)):
                if (rr2[m][0] != k):
                    continue
                n = rr2[m][1]
                c[3][n][i] = rr2[m][2][i]
            for m in range(len(re3)):
                if (re3[m][0] != k):
                    continue
                n = re3[m][1]
                c[2][n][i] = re3[m][2][i]
            for m in range(len(ci3)):
                if (ci3[m][0] != k):
                    continue
                n = ci3[m][1]
                c[9][n][i] = ci3[m][2][i]
            for m in range(len(cs3)):
                if (cs3[m][0] != k):
                    continue
                n = cs3[m][1]
                c[6][n][i] = (cs3[m][2][i]-c[5][n][i]-p*c[8][n][i])
                if (c[6][n][i] < 0.0):
                    c[6][n][i] = 0.0
        for m in range(len(ce1)):
            for i in range(nt):
                col = [ce1[m][0], ce1[m][1]]
                col.append(logt[i])
                for j in range(10):
                    if (j == 4 or j == 5):
                        c[j][m][i] = c[j][m][i]+c[1][m][i]
                    elif (j == 6):
                        c[j][m][i] = c[j][m][i]+c[2][m][i]
                    elif (j == 7 or j == 8):
                        c[j][m][i] = c[j][m][i]+c[3][m][i]
                    col.append(c[j][m][i])
                col = tuple(col)
                tbl.write_row(*col)
        
    tbl.close()
    
        
def write_trates(f, r, header, nele):
    s = '# %s\n'%(header)
    f.write(s)
    for i in range(len(r)):
        if (r[i][0] == nele):
            continue
        s = '%2d '%(r[i][0])
        for a in r[i][1]:
            s = s + '%9.3E '%(a)
        s = s[:-1] + '\n'
        f.write(s)
    f.write('\n')

def write_rates(f, r, header, nele):
    s = '# %s\n'%(header)
    f.write(s)
    for i in range(len(r)):
        if (r[i][0] == nele):
            continue
        for j in range(2, len(r[i])):
            s = '%2d %4d  '%(r[i][0], r[i][1])
            for a in r[i][j]:
                s = s + '%9.3E '%(a)
            s = s[:-1] + '\n'
            f.write(s)
    f.write('\n')
    

def save_rates(rates, sfile, dfile, **kwd):
    rates.update(kwd)
    rates['tdr'] = copy.deepcopy(rates['tdc'])
    for i in range(1, len(rates['tdc'])):
        rates['tdr'][i][1] = map(lambda x,y: x-y,
                                 rates['tdc'][i][1],
                                 rates['tre'][i-1][1])
    with open(sfile, 'w') as f:
        pickle.dump(rates, f)

    f = open(dfile, 'w')
    if (rates.has_key('temp')):
        temp = rates['temp']
    else:
        temp = []
    if (rates.has_key('logt')):
        logt = rates['logt']
    else:
        logt = []
    if (rates.has_key('abund')):
        abund = rates['abund']
    else:
        abund = []
        
    neles = map(lambda x:x[0], rates['tdc'])
    nt = len(temp)
    s = '#   LogT(K)  Temp(eV)   '
    for k in neles:
        s = s + 'NELE=%-5d '%(k)
    s = s[:-1] + '\n'
    f.write(s)
    for i in range(nt):
        s = '%2d '%i
        if (logt):
            s = s + '%7.3f  '%(logt[i])
        s = s + '%10.4E '%(temp[i])
        if (abund):
            for k in neles:
                s = s + '%10.4E '%(abund[i][k])
        s = s[:-1]+'\n'
        f.write(s)
    f.write('\n')
    
    if (rates.has_key('tdc')):
        write_trates(f, rates['tdc'],
                     'Total Dielectronic Capture', neles[0])
    if (rates.has_key('tdr')):
        write_trates(f, rates['tdr'],
                     'Total Dielectronic Recombination', neles[0])
    if (rates.has_key('trr')):
        write_trates(f, rates['trr'],
                     'Total Radiative Recombination', neles[0])
    if (rates.has_key('tea')):
        write_trates(f, rates['tea'],
                     'Total Excitation Autoionization', neles[-1])
    if (rates.has_key('tci')):
        write_trates(f, rates['tci'],
                     'Total Direct Ionization', neles[-1])
    if (rates.has_key('rt')):
        write_rates(f, rates['rt'],
                    'Total Depletion Rate', -1)
    if (rates.has_key('cs')):
        write_rates(f, rates['cs'],
                    'Radiative Cascades', -1)
    if (rates.has_key('ce')):
        write_rates(f, rates['ce'],
                    'Direct Excitation', -1)
    if (rates.has_key('re')):
        write_rates(f, rates['re'],
                    'Resonance Excitation', neles[-1])
    if (rates.has_key('ci')):
        write_rates(f, rates['ci'],
                    'Direct Ionization', neles[-1])
    if (rates.has_key('rr')):
        write_rates(f, rates['rr'],
                    'Radiative Recombination', neles[0])

    f.close()
    
    
def read_rates(nt, nd, nele, pref='Fe', dir='', nion=2, only_total=0):
    c = get_complexes(nele)
    complexes = [c[1]]
    if (nion > 1):
        c = get_complexes(nele-1)
        complexes.append(c[1])
        if (nion > 2):
            c = get_complexes(nele+1)
            complexes.append(c[1])
    re = []
    ci = []
    rr = []
    cs = []
    ce = []
    rt = []
    tdc = []
    tre = []
    trr = []
    tpi = []
    tci = []
    tea = []
    for t in range(nt):
        for d in range(nd):
            rt_file = dir+'%s%02da_t%02dd%di%d.rt'%(pref, nele, t, d, nion)
            f = open(rt_file, 'r')
            ilev = -1
            den = [-1]*3
            ire = 0
            ici = 0
            irr = 0
            ics = 0
            ice = 0
            irt = 0
            itot = 0
            while (1):
                a = f.readline()
                if (not a):
                    break
                if (len(a) < 4):
                    continue
                if (a[:4] == 'NELE'):
                    a = string.split(a)
                    nel = int(a[2])
                    bsum = 1
                elif (a[:4] == 'NTRA'):
                    a = string.split(a)
                    ntrans = int(a[2])
                elif (a[:4] == 'ILEV'):
                    a = string.split(a)
                    ilev = int(a[2])
                    nilev = 1
                    r3 = 0.0
                elif (a[:4] == ' SUM'):
                    a = string.split(a)
                    b = []
                    for c in a[1:-1]:
                        b.append(float(c))
                    if (den[0] > 0):
                        b[0] = b[0]/den[0]
                        b[1] = b[1]/den[0]
                    if (den[1] > 0):
                        b[2] = b[2]/den[1]
                    if (den[2] > 0):
                        b[3] = b[3]/den[2]
                        b[4] = b[4]/den[2]
                        b[5] = b[5]/den[2]
                    if (t == 0 and d == 0):
                        trr.append([nel, [b[0]]])
                        tdc.append([nel, [b[1]]])
                        tre.append([nel, [b[2]]])
                        tpi.append([nel, [b[3]]])
                        tea.append([nel, [b[4]]])
                        tci.append([nel, [b[5]]])
                    else:
                        trr[itot][1].append(b[0])
                        tdc[itot][1].append(b[1])
                        tre[itot][1].append(b[2])
                        tpi[itot][1].append(b[3])
                        tea[itot][1].append(b[4])
                        tci[itot][1].append(b[5])
                    itot = itot + 1
                elif (ilev >= 0):
                    if (only_total > 0):
                        continue
                    if (a[:4] == 'DENS'):
                        a = string.split(a)
                        d0 = float(a[2])
                        tp = 0.0
                    elif (a[2:4] == '-1'):
                        a = string.split(a)
                        den[0] = float(a[1])
                        if (not den[0]):
                            continue
                        b = float(a[4])
                        if (b > 0):
                            tp = tp + b
                            b = b/den[0]
                            if (t == 0 and d == 0):
                                rr.append([nel, ilev, [b]])
                            else:
                                rr[irr][2].append(b)
                            irr = irr + 1
                    elif (a[2:4] == '-2'):
                        a = string.split(a)
                        den[1] = float(a[1])
                        if (not den[1]):
                            continue
                        b = float(a[2])
                        if (b > 0):
                            tp = tp + b
                            b = b/den[1]
                            r3 = r3/den[1]
                            if (t == 0 and d == 0):
                                cs.append([nel, ilev, [b], [r3]])
                            else:
                                cs[ics][2].append(b)
                                cs[ics][3].append(r3)
                            ics = ics + 1
                        b = float(a[3])
                        if (b > 0):
                            tp = tp + b
                            b = b/den[1]
                            if (t == 0 and d == 0):
                                ce.append([nel, ilev, [b]])
                            else:
                                ce[ice][2].append(b)
                            ice = ice + 1
                        b = float(a[5])
                        if (b > 0):
                            tp = tp + b
                            b = b/den[1]
                            if (t == 0 and d == 0):
                                re.append([nel, ilev, [b]])
                            else:
                                re[ire][2].append(b)
                            ire = ire + 1
                    elif (a[2:4] == '-3'):
                        a = string.split(a)
                        den[2] = float(a[1])
                        if (den[2] > 0):
                            b = float(a[6])
                            if (b > 0):
                                tp = tp + b
                                b = b/den[2]
                                if (t == 0 and d == 0):
                                    ci.append([nel, ilev, [b]])
                                else:
                                    ci[ici][2].append(b)
                                ici = ici + 1
                        if (ntrans > 3):
                            tp = tp/d0
                        if (t == 0 and d == 0):
                            rt.append([nel, ilev, [tp]])
                        else:
                            rt[irt][2].append(tp)
                        irt = irt + 1
                        bsum = 0
                    else:
                        if (not bsum):
                            continue
                        c = a[72:-1]
                        if (len(c) > 1):
                            c = string.strip(c)
                            if c in complexes:
                                a = string.split(a[:72])
                                b = float(a[2])
                                if (b > 0):
                                    r3 = r3 + b
                                
                            
    def compare(x, y):
        if (x[0] < y[0]):
            return -1
        elif (x[0] > y[0]):
            return 1
        else:
            if (x[1] < y[1]):
                return -1
            elif (x[1] > y[1]):
                return 1
            else:
                return 0
    cs.sort(compare)
    ce.sort(compare)
    re.sort(compare)
    ci.sort(compare)
    rr.sort(compare)
    rt.sort(compare)
    return {'tdc': tdc,
            'tre': tre,
            'trr': trr,
            'tpi': tpi,
            'tea': tea,
            'tci': tci,
            'rt':  rt,
            'cs':  cs,
            'ce':  ce,
            'rr':  rr,
            're':  re,
            'ci':  ci}

                    
def get_tgrid(z, nele, dt = 0.15, amin = 1E-2, limits=[]):
    if (len(limits) == 0):
        limits = [5.0, 8.0]
        
    (tmax,a) = MaxAbund(z, nele)
    amax = a[nele-1:nele+1]
    logtm = log10(tmax/const.kb)
    logtm = limits[0]+int((logtm-limits[0])/dt)*dt
    a0 = 1.0
    logt = [logtm]
    t = [const.kb*10**(logtm)]
    a = FracAbund(z, t[0])
    ab = [a]
    while (a0 > amin and logt[-1] < limits[1]):
        logt0 = logt[-1] + dt
        t0 = const.kb*10**(logt0)
        t.append(t0)
        logt.append(logt0)
        a = FracAbund(z, t0)
        ab.append(a)
        a0 = max(map(lambda x,y:x/y, a[nele-1:nele+1], amax))
        
    a0 = 1.0
    while (a0 > amin and logt[0] > limits[0]):
        logt0 = logt[0] - dt
        t0 = const.kb*10**(logt0)
        t.insert(0, t0)
        logt.insert(0, logt0)
        a = FracAbund(z, t0)
        ab.insert(0, a)
        a0 = max(map(lambda x,y:x/y, a[nele-1:nele+1], amax))
        
    return (t, logt, ab)

    
def get_complexes(nelectrons):
    n = 1
    nele = nelectrons
    g = []
    while (nele > 0):
        nqm = 2*n*n
        if (nele < nqm):
            g.append((n,nele))
        else:
            g.append((n,nqm))
        nele = nele - nqm
        n = n + 1
    c0 = ''
    c1 = ''
    for a in g:
        if (a != g[-1]):
            c0 = c0 + '%d*%d '%a
            c1 = c1 + '%d*%d '%a
        else:
            c0 = c0 + '%d*%d'%a
            if (a[1] == 1):
                c1 = c1 + '%d*%d'%(a[0]+1,1)
            else:
                c1 = c1 + '%d*%d %d*%d'%(a[0], a[1]-1, a[0]+1, 1)
    return [c0, c1]

def spectrum(neles, temp, den, population, pref,
             suf='b', osuf='', dir0 = '', dir1= '', nion = 3,
             dist = 0, params=[-1,-1], cascade = 0, rrc = 0, ion0 = 1, 
             abund0 = 1.0, abundm = -1, abundp = -1, iprint=1,
             frr0 = -1, fci0 = -1, frrp = -1, fcip = -1,
             ai = 1, ci = 1, rr = 1, ce = 1, eps = 1E-4, rcomp = [],
             t0=-1, t1=-1, d0=-1, d1=-1,
             mtr='', mce='', mci='', mai='', mrr=''):
    if t0 < 0:
        t0 = 0
    if t1 < 0:
        t1 = len(temp)-1
    if d0 < 0:
        d0 = 0
    if d1 < 0:
        d1 = len(den)-1
    for k in neles:
        rate = get_complexes(k)
        if (nion > 1):
            rate = rate + get_complexes(k-1)
            if (nion > 2):
                rate = rate + get_complexes(k+1)
        if (len(rcomp) > 0):
            rate = rate + rcomp
        print('NELE = %d'%k)
        f1 = '%s%02d%s'%(pref, k-1, suf)
        f2 = '%s%02d%s'%(pref, k, suf)
        f3 = '%s%02d%s'%(pref, k+1, suf)
        AddIon(k, 0.0, dir0+f2)
        if (frr0 > 0):
            SetRateMultiplier(k, 3, frr0)
        if (fci0 > 0):
            SetRateMultiplier(k, 2, fci0)
        if (nion == 3):
            AddIon(k+1, 0.0, dir0+f3)
            if (frrp > 0):
                SetRateMultiplier(k+1, 3, frrp)
            if (fcip > 0):
                SetRateMultiplier(k+1, 2, fcip)
        if (nion > 1 and abundm > -10):
            if (k > 1 and ion0 > 0):
                SetBlocks(0.0, dir0+f1)
            else:
                SetBlocks(0.0)
        else:
            SetBlocks(-1.0)

        print('TR rates...')
        SetTRRates(0)
        if mtr != '':
            ModifyRates(mtr)
        for i in range(t0, t1+1):
            if (abundm < 0 or abundp < 0):
                p1 = population[i][k-1]
                p2 = population[i][k]
                try:
                    p3 = population[i][k+1]
                except:
                    p3 = p2
            if (abundm < 0):
                p1 = abund0*(p1/p2)
            else:
                p1 = abundm
            if (abundp < 0):
                p3 = abund0*(p3/p2)
            else:
                p3 = abundp
            p2 = abund0
            print('Temp = %10.3E'%(temp[i]))
            print('Abund: %10.3E %10.3E %10.3E'%(p1, p2, p3))

            dp = [dist, temp[i]]
            dp[2:] = params
            SetEleDist(*dp)
            
            if (ce > 0):
                print('CE rates...')
                SetCERates(1)
                if (mce != ''):
                    ModifyRates(mce+'.t%02d'%i)
            if (nion > 1):
                if (rr > 0):
                    print('RR rates...')
                    SetRRRates(0)
                    if (mrr != ''):
                        ModifyRates(mrr+'.t%02d'%i)
                if (ci > 0):
                    print('CI rates...')
                    SetCIRates(1)
                    if (mci != ''):
                        ModifyRates(mci+'.t%02d'%i)
                elif (ci < 0):
                    print('CI rates...')
                    SetCIRates(0)        
                    if (mci != ''):
                        ModifyRates(mci+'.t%02d'%i)            
                if (ai > 0):
                    print('AI rates...')
                    SetAIRates(1)
                    if (mai != ''):
                        ModifyRates(mai+'.t%02d'%i)
                elif (ai < 0):
                    print('AI rates...')
                    SetAIRates(0)
                    if (mai != ''):
                        ModifyRates(mai+'.t%02d'%i)
                SetAbund(k-1, p1)
            else:
                SetAbund(k-1, 0.0)
            SetAbund(k, p2)
            if (nion == 3):
                SetAbund(k+1, p3)
                
            for d in range(d0, d1+1):
                print('Density = %10.3E'%den[d])
                SetEleDensity(den[d])
                SetIteration(eps)
                SetCascade(cascade, eps)
                print('Init blocks...')
                InitBlocks()
                s = 't%02dd%di%d%s'%(i, d, nion, osuf)
                if (abundm > 0):
                    s = s + 'm'
                if (abundp > 0):
                    s = s + 'p'
                rt_file = dir1+'%s_%s.rt'%(f2,s)
                sp_file = dir1+'%s_%s.sp'%(f2,s)
                rt_afile = dir1+'%sa_%s.rt'%(f2[:-1],s)
                sp_afile = dir1+'%sa_%s.sp'%(f2[:-1],s)
                dp_afile = dir1+'%sa_%s.d'%(f2[:-1],s)
                LevelPopulation()
                Cascade()
                for dtp in range(6):
                    DumpRates('%s%d'%(dp_afile,dtp), k, dtp, -1, 1)
                RateTable(rt_file, rate)
                SpecTable(sp_file, rrc)
                if iprint:
                    PrintTable(rt_file, rt_afile, 1)
                    PrintTable(sp_file, sp_afile, 1)
                sys.stdout.flush()
                ReinitCRM(3)
            ReinitCRM(2)
        ReinitCRM(0)


def maxwell(e, t):
    x = e/t
    x = 1.12837967*sqrt(x)*exp(-x)/t
    return x

