from pfac.spm import *
from pfac.crm import *

nt = 9
logt = range(nt)
logt = map(lambda x:x*0.18+6.3, logt)
temp = map(lambda x:(10**x)*const.kb, logt)

t = nt-1
b = [[]]*3
emin = [0]*3
emax = [0]*3

z = 26.0
smin = 3E-1
trans = (5, 2)
#smin = 1E-1
#trans = (4, 2)
#smin = 5E-2
#trans = (3, 2)
e1 = 0.5E3
e2 = 2.0E3
de = 0.75
type = trans[0]*100+trans[1]
tp = '%02d%02d'%trans

for k in range(3, 12):
    pref = 'Fe%02d'%k
    if (k == 11):
        nion = 2
    else:
        nion = 3
    spf = '%sb_t%dd0i%d.sp'%(pref, t, 1)
    lnf = '%sa_t%dd0i%d.ln%s'%(pref, t, 1, tp)
    pltf0 = '%sa_t%dd0i%d.plt%st'%(pref, t, 1, tp)
    pltf = '%sa_t%dd0i%d.plt%s'%(pref, t, 1, tp)
    SelectLines(spf, lnf, k, type, e1, e2, smin)
    print 'NELE = %d %s %s'%(k, spf, lnf)
    a = read_lines(lnf, k)
    c = map(lambda x:x[2], a)
    if (len(c) > 0):
        d1 = 0.95*min(c)
        d2 = 1.05*max(c)
    else:
        d1 = emin[-1]
        d2 = emax[-1]
    emin.append(d1)
    emax.append(d2)
    b.append(a)
    PlotSpec(spf, pltf, k, type, emin[k], emax[k], de)
    PlotSpec(spf, pltf0, k, 2000100, emin[k], emax[k], de)

for k in range(2, 12):
    pref = 'Fe%02d'%k
    if (k == 11):
        nion = 2
    else:
        nion = 3
    for t in range(nt):
        spf = '%sb_t%dd0i%d.sp'%(pref, t, nion)
        lnf = '%sa_t%dd0i%d.ln%s'%(pref, t, nion, tp)
        print 'NELE = %d T = %d %s %s'%(k, t, spf, lnf)
        pltf0 = '%sa_t%dd0i%d.plt%st'%(pref, t, nion, tp)
        pltf = '%sa_t%dd0i%d.plt%s'%(pref, t, nion, tp)
        if (k > 3):
            a = b[k-1]
            PlotSpec(spf, pltf+'-', k-1, type, emin[k-1], emax[k-1], de)
            PlotSpec(spf, pltf0+'-', k-1, 2000100, emin[k-1], emax[k-1], de)
            for a0 in a:
                h1 = float(a0[0])
                h2 = float(a0[1])
                SelectLines(spf, lnf, k-1, type, h1, h2, -1.0)
        if (k > 2):
            a = b[k]
            PlotSpec(spf, pltf, k, type, emin[k], emax[k], de)
            PlotSpec(spf, pltf0, k, 2000100, emin[k], emax[k], de)
            if (k > 3):
                for ty in range(3,8):
                    ty1 = ty*10000+type
                    pltf1 = '%s_%d'%(pltf, ty)
                    lnf1 = '%s_%d'%(lnf, ty)
                    if (ty == 7):
                        ty1 = int(1E6) + ty1
                    PlotSpec(spf, pltf1, k, ty1,
                             emin[k-1], emax[k-1], de)
                    SelectLines(spf, lnf1, k, ty1,
                                emin[k-1], emax[k-1], 1E-3)
            for a0 in a:
                h1 = float(a0[0])
                h2 = float(a0[1])
                SelectLines(spf, lnf, k, type, h1, h2, -1.0)  
        if (k < 11):
            a = b[k+1]
            PlotSpec(spf, pltf+'+', k+1, type, emin[k+1], emax[k+1], de)
            PlotSpec(spf, pltf0+'+', k+1, 2000100, emin[k+1], emax[k+1], de)
            if (k > 2):
                for ty in range(3,8):
                    ty1 = ty*10000+type
                    pltf1 = '%s_%d+'%(pltf, ty)
                    lnf1 = '%s_%d'%(lnf, ty)
                    if (ty == 7):
                        ty1 = int(1E6) + ty1
                    PlotSpec(spf, pltf1, k+1, ty1, emin[k], emax[k], de)
                    SelectLines(spf, lnf1, k+1, ty1, emin[k], emax[k], 1E-3)
            for a0 in a:
                h1 = float(a0[0])
                h2 = float(a0[1])
                SelectLines(spf, lnf, k+1, type, h1, h2, -1.0)
                
        if (t < nt-1 and k > 2):
            spf = '%sb_t%dd0i%d.sp'%(pref, t, 1)
            lnf = '%sa_t%dd0i%d.ln%s'%(pref, t, 1, tp)
            pltf0 = '%sa_t%dd0i%d.plt%st'%(pref, t, 1, tp)
            pltf = '%sa_t%dd0i%d.plt%s'%(pref, t, 1, tp)
            a = b[k]
            PlotSpec(spf, pltf, k, type, emin[k], emax[k], de)
            PlotSpec(spf, pltf0, k, 2000100, emin[k], emax[k], de)
            for a0 in a:
                h1 = float(a0[0])
                h2 = float(a0[1])
                SelectLines(spf, lnf, k, type, h1, h2, -1.0)
        

