from pfac.crm import *
from pfac.spm import *
from biggles import *
from Numeric import *
from math import *
from pfac import const
import cPickle

def get_spec(t0, emin, emax, de, ide, nion, z=26.0):
    d = 0.1*de
    d1 = 0.1*d
    n0 = (emax-emin)/d + 1
    nt0 = len(t0)
    logt0 = map(lambda a:log10(a/const.kb), t0)
    abund0 = []
    for t in t0:
        abund0.append(FracAbund(z, t))
        
    x0 = zeros((n0,), Float)
    y0 = zeros((nt0, n0), Float)
    n = 1
    x0[0] = emin
    while (n < n0):
        x0[n] = x0[n-1] + d
        n = n + 1
    for k in range(3, 11):
        pref = 'Fe%02d/Fe%02d'%(k, k)
        rates1 = cPickle.load(open('Fe%02d/rates1.sav'%k))
        temp = rates1['temp']
        logt = rates1['logt']
        abund = rates1['abund']
        nt = len(temp)
        y = []
        for t in range(nt):
            pltf = '%sa_t%02di%dw%d.plt'%(pref, t, nion, ide)
            if (t == 0):
                x = read_column(0, pltf)
            y.append(read_column(1, pltf))
        nx = len(x)
        i0 = 0
        while (1):
            if (fabs(x[0] - x0[i0]) < d1):
                break
            i0 = i0 + 1
        i = 0
        while (i < nx):
            i1 = i0 + i
            yt = map(lambda a,p:a[p], y, [i]*nt)
            y2 = Spline(logt, yt)
            for j in range(nt0):
                a = Splint(logt, yt, y2, logt0[j])
                y0[j][i1] = y0[j][i1] + a*abund0[j][k]
            i = i + 1
        if (k == 10 and nion == 3):
            y = []
            for t in range(nt):
                pltf = '%sa_t%02di%dw%d.plt+1'%(pref, t, nion, ide)
                if (t == 0):
                    x = read_column(0, pltf)
                y.append(read_column(1, pltf))
            nx = len(x)
            i0 = 0
            while (1):
                if (fabs(x[0] - x0[i0]) < d1):
                    break
                i0 = i0 + 1
            i = 0
            while (i < nx):
                i1 = i0 + i
                yt = map(lambda a,p:a[p], y, [i]*nt)
                y2 = Spline(logt, yt)
                for j in range(nt0):
                    a = Splint(logt, yt, y2, logt0[j])
                    y0[j][i1] = y0[j][i1] + a*abund0[j][k]
                i = i + 1
            
    return (x0, y0)

        
emin = 600.0
emax = 1.8e3
t0 = [4e2, 6e2, 8e2, 1e3]
nt0 = len(t0)
de = [1.0, 10.0]
biggles.configure('fontsize_min', 1.0)
for i in [0,1]:
    p = FramedArray(nt0, 3)
    p.xlabel = 'Wavelength ({\AA})'
    p.ylabel = 'Emissivity (s$^{-1}$)'
    p.cellspacing = 0
    p.uniform_limits = 0
    p.aspect_ratio = 0.75
    p.xrange = (10.5, 18)
    for j in [3,2,1]:
        (x0, y0) = get_spec(t0, emin, emax, de[i], i, j)
        x0 = map(lambda a:const.hc/a, x0)
        y0 = map(lambda a:10*a, y0)
        for t in range(nt0):
            c = Curve(x0, y0[t], linewidth=0.6)
            if (i == 0):
                if (t < 2):
                    p[t, j-1].add(LineY(1, type='dot', width=0.8))
                    p[t, j-1].add(LineY(2, type='dot', width=0.8))
                else:
                    p[t, j-1].add(LineY(0.2, type='dot', width=0.8))
                    p[t, j-1].add(LineY(0.4, type='dot', width=0.8))
                    p[t, j-1].add(LineY(0.6, type='dot', width=0.8))
            else:
                if (t == 0):
                    p[t, j-1].add(LineY(0.4, type='dot', width=0.8))
                elif (t == 1):
                    p[t, j-1].add(LineY(0.1, type='dot', width=0.8))
                    p[t, j-1].add(LineY(0.3, type='dot', width=0.8))
                else:
                    p[t, j-1].add(LineY(0.05, type='dot', width=0.8))
                    p[t, j-1].add(LineY(0.1, type='dot', width=0.8))
                
            p[t, j-1].add(c)
            p[t, j-1].add(PlotLabel(0.02, 0.9,
                                    'Te=%3.1f'%(t0[t]/1e3),
                                    size=0.8, halign='left'))
            p[t, j-1].add(PlotLabel(0.02, 0.8,
                                    'Nion=%d'%(j),
                                    size=0.8, halign='left'))
    epsf = 'spec_w%d.eps'%(i)
    p.write_eps(epsf)
    
