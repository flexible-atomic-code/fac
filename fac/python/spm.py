from pfac.crm import *
from pfac import const
from math import *
import sys
import string
import biggles

def read_rates(nt, nd, nele, pref='Fe', dir='', nion=2):
    re = []
    ci = []
    rr = []
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
            itot = 0
            while (1):
                a = f.readline()
                if (not a):
                    break
                if (len(a) < 4):
                    continue
                if (a[:4] == 'NELE'):
                    a = string.split(a)
                    nele = int(a[2])
                if (a[:4] == 'ILEV'):
                    a = string.split(a)
                    ilev = int(a[2])
                if (a[:4] == ' Sum'):
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
                        trr.append([nele, [b[0]]])
                        tdc.append([nele, [b[1]]])
                        tre.append([nele, [b[2]]])
                        tpi.append([nele, [b[3]]])
                        tea.append([nele, [b[4]]])
                        tci.append([nele, [b[5]]])
                    else:
                        trr[itot][1].append(b[0])
                        tdc[itot][1].append(b[1])
                        tre[itot][1].append(b[2])
                        tpi[itot][1].append(b[3])
                        tea[itot][1].append(b[4])
                        tci[itot][1].append(b[5])
                    itot = itot + 1
                elif (ilev >= 0):
                    if (a[2:4] == '-1'):
                        a = string.split(a)
                        den[0] = float(a[1])
                        if (not den[0]):
                            continue
                        b = float(a[4])
                        if (b > 0):
                            b = b/den[0]
                            if (t == 0 and d == 0):
                                rr.append([nele, ilev, [b]])
                            else:
                                rr[irr][2].append(b)
                            irr = irr + 1
                    if (a[2:4] == '-2'):
                        a = string.split(a)
                        den[1] = float(a[1])
                        if (not den[1]):
                            continue
                        b = float(a[5])
                        if (b > 0):
                            b = b/den[1]
                            if (t == 0 and d == 0):
                                re.append([nele, ilev, [b]])
                            else:
                                re[ire][2].append(b)
                            ire = ire + 1
                    if (a[2:4] == '-3'):
                        a = string.split(a)
                        den[2] = float(a[1])
                        if (not den[2]):
                            continue
                        b = float(a[6])
                        if (b > 0):
                            b = b/den[2]
                            if (t == 0 and d == 0):
                                ci.append([nele, ilev, [b]])
                            else:
                                ci[ici][2].append(b)
                            ici = ici + 1

    return {'tdc': tdc,
            'tre': tre,
            'trr': trr,
            'tpi': tpi,
            'tea': tea,
            'tci': tci,
            'rr':  rr,
            're':  re,
            'ci':  ci}

                    
def get_tgrid(z, nele, dt = 0.15, amin = 5E-2, limits=[]):
    if (len(limits) == 0):
        limits = [100.0, 7e3]
        
    (tmax,a) = MaxAbund(z, nele)
    amax = a[nele-1:nele+1]
    logtm = log10(tmax/const.kb)
    logtm = int(logtm*10)/10.0
    a0 = 1.0
    logt = [logtm]
    t = [const.kb*10**(logtm)]
    a = FracAbund(z, t[0])
    ab = [a]
    while (a0 > amin and t[-1] < limits[1]):
        logt0 = logt[-1] + dt
        t0 = const.kb*10**(logt0)
        t.append(t0)
        logt.append(logt0)
        a = FracAbund(z, t0)
        ab.append(a)
        a0 = max(map(lambda x,y:x/y, a[nele-1:nele+1], amax))
        
    a0 = 1.0
    while (a0 > amin and t[0] > limits[0]):
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
    return (c0, c1)

def spectrum(neles, temp, den, population, pref,
             suf = 'b', dir0 = '', dir1= '', nion = 3,
             dist = 0, cascade = 0):
    for k in neles:
        rate = get_complexes(k)
        if (nion > 1):
            rate = rate + get_complexes(k-1)
            if (nion > 2):
                rate = rate + get_complexes(k+1)
                
        print 'NELE = %d'%k
        f1 = '%s%02d%s'%(pref, k-1, suf)
        f2 = '%s%02d%s'%(pref, k, suf)
        f3 = '%s%02d%s'%(pref, k+1, suf)
        AddIon(k, 0.0, dir0+f2) 
        if (nion == 3):
            AddIon(k+1, 0.0, dir0+f3)
        if (nion > 1):
            SetBlocks(0.0, dir0+f1)
        else:
            SetBlocks(-1.0)

        for i in range(len(temp)):
            p1 = population[i][k-1]
            p2 = population[i][k]
            p3 = population[i][k+1]
            p1 = p1/p2
            p3 = p3/p2
            p2 = 1.0
            print 'Temp = %10.3E'%(temp[i])
            print 'Abund: %10.3E %10.3E %10.3E'%(p1, p2, p3)

            SetEleDist(dist, temp[i], -1.0, -1.0)
            print 'CE rates...'
            SetCERates(1)
            print 'TR rates...'
            SetTRRates(0)
            if (nion > 1):
                print 'RR rates...'
                SetRRRates(0)
                print 'CI rates...'
                SetCIRates(0)
                print 'AI rates...'
                SetAIRates(1)
                SetAbund(k-1, p1)
            SetAbund(k, p2)
            if (nion == 3):
                SetAbund(k+1, p3)
                
            for d in range(len(den)):
                print 'Density = %10.3E'%den[d]
                SetEleDensity(den[d])
                SetCascade(cascade, 1E-4)
                print 'Init blocks...'
                InitBlocks()
                s = 't%02dd%di%d'%(i, d, nion)
                rt_file = dir1+'%s_%s.rt'%(f2,s)
                sp_file = dir1+'%s_%s.sp'%(f2,s)
                rt_afile = dir1+'%sa_%s.rt'%(f2[:-1],s)
                sp_afile = dir1+'%sa_%s.sp'%(f2[:-1],s)

                LevelPopulation()
                Cascade()

                rt = (rt_file,)+rate
                RateTable(*rt)
                SpecTable(sp_file)
                PrintTable(rt_file, rt_afile, 1)
                PrintTable(sp_file, sp_afile, 1)
                sys.stdout.flush()
                ReinitCRM(2)
            ReinitCRM(1)
        ReinitCRM()

def read_lines(file, nele):
    k = -1
    s = []
    for line in open(file, 'r').readlines():
        a = string.split(line)
        if (len(a) != 5):
            continue
        if (int(a[0]) != nele):
            continue
        e = float(a[3])
        i = int(a[1])
        j = int(a[2])
        s0 = float(a[4])
        s.append([i, j, e, s0])
                
    return s

def id_lines(w, s, name, eps=''):
    n = len(w)
    biggles.configure('screen', 'width', 800)
    biggles.configure('screen', 'height', 500)
    p = biggles.FramedPlot()
    p.aspect_ratio = 0.7
    xmin = min(w)
    xmax = max(w)
    for i in range(n):
        q1 = (w[i], 0)
        q2 = (w[i], s[i])
        l = biggles.Line(q1, q2)
        p.add(l)
        l = biggles.Label(q2[0], q2[1], '%d: %s'%(i,name[i]),
                          color='red', textangle=90)
        p.add(l)
    p.add(biggles.Line((xmin,0), (xmax,0), color='red'))
    if (eps):
        p.write_eps(eps)
    else:
        p.show()
