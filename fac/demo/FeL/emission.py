from pfac.spm import *
from pfac import crm
from pfac import const
import string
import Gnuplot
import copy
from math import *

nt = 9
logt = range(nt)
logt = map(lambda x:x*0.18+6.3, logt)
temp = map(lambda x:(10**x)*const.kb, logt)

trans = (5, 2)
#trans = (4, 2)
#trans = (3, 2)
type = trans[0]*100+trans[0]
tp = '%02d%02d'%trans

z = 26
population = []
for i in range(len(temp)):
    a = FracAbund(z, temp[i])
    population.append(a)

s = [[],[]]
m = 3
for k in range(2,12):
    pref = 'Fe%02d'%k
    if (k == 11):
        nion = 2
    else:
        nion = 3
    for t in range(nt):
        lnf = '%sa_t%dd0i%d.ln%s'%(pref, t, nion, tp)
        a = read_lines(lnf, k)
        if (t == 0):
            b = a[:]
            for b0 in b:
                b0[m] = [b0[m]]
        else:
            for b0 in b:
                for a0 in a:
                    if (b0[0] == a0[0] and b0[1] == a0[1]):
                        b0[m].append(a0[3])
    c = []
    for b0 in b:
        if (len(b0[3]) == nt):
            c.append(b0)               
    s.append(c)

m = 4
for k in range(3, 11):
    pref = 'Fe%02d'%k
    b = s[k]
    for b0 in b:
        b0.append([])
    for t in range(nt):
        lnf = '%sa_t%dd0i%d.ln%s'%(pref, t, 1, tp)
        a = read_lines(lnf, k)
        for b0 in b:
            inew = 0
            for a0 in a:
                if (b0[0] == a0[0] and b0[1] == a0[1]):
                    b0[m].append(a0[3])
                    inew = 1
            if (inew == 0):
                b0[m].append(0.0)

m = 5                   
for k in range(2, 10):
    pref = 'Fe%02d'%k
    if (k == 11):
        nion = 2
    else:
        nion = 3
    b = s[k+1]
    for b0 in b:
        b0.append([])
    for t in range(nt):
        x = population[t][k+1]/population[t][k]
        lnf = '%sa_t%dd0i%d.ln%s'%(pref, t, nion, tp)
        a = read_lines(lnf, k+1)
        for b0 in b:
            inew = 0
            for a0 in a:
                if (b0[0] == a0[0] and b0[1] == a0[1]):
                    b0[m].append(a0[3]/x)
                    inew = 1
            if (inew == 0):
                b0[m].append(0.0)

m = 6
for k in range(4, 12):
    pref = 'Fe%02d'%k
    if (k == 11):
        nion = 2
    else:
        nion = 3
    b = s[k-1]
    for b0 in b:
        b0.append([])
    for t in range(nt):
        x = population[t][k-1]/population[t][k]
        lnf = '%sa_t%dd0i%d.ln%s'%(pref, t, nion, tp)
        a = read_lines(lnf, k-1)
        for b0 in b:
            inew = 0
            for a0 in a:
                if (b0[2] == a0[2]):
                    b0[m].append(a0[3]/x)
                    inew = 1
            if (inew == 0):
                b0[m].append(0.0)

f = open('emiss%s.dat'%tp, 'w')
m = -1
for k in range(3,11):
    enf = '../Fe%02db.en'%k
    s0 = s[k]
    s0.sort(lambda x,y:int((y[2]-x[2])/(0.5*fabs(y[2]-x[2]))))
    smax = max(map(lambda x:x[4][-1], s0))
    f.write('#NELE = %d\n'%(k))
    for i in range(len(s0)):
        if (s0[i][4][-1] < 0.05*smax):
            continue
        if (0.0 in s0[i][3] or 0.0 in s0[i][4] or 0.0 in s0[i][5]):
            continue
        e = s0[i][2]
        w = const.hc/e
        n1 = crm.LevelName(enf, s0[i][0])
        n2 = crm.LevelName(enf, s0[i][1])
        m = m+1
        a = '#%4d\t%10.4E %4d %4d\n'%(m, w, s0[i][0], s0[i][1])
        f.write(a)
        a = '# %-20s %-20s %-s\n'%n1
        f.write(a)
        a = '# %-20s %-20s %-s\n'%n2
        f.write(a)
        for j in range(nt):
            c0 = s0[i][4][j]
            c1 = s0[i][5][j]
            c2 = s0[i][3][j]
            b = (logt[j], temp[j], c0, c1, c1/c0, c2, c2/c0)
            a = '%10.4E %10.4E  %10.4E %10.4E %10.4E %10.4E %10.4E\n'%b
            f.write(a)
        f.write('\n\n')

f.close()
