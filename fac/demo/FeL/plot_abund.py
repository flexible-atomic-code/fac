import cPickle
import biggles
from pfac import const
from pfac.table import *
from pfac.spm import *
from math import *

biggles.configure('fontsize_min', 1.0)
p = biggles.FramedPlot()
p.ylog = 1
p.ylabel = 'Abundance Ratio'
p.xlabel = 'Log Temperature (K)'
p.aspect_ratio = 0.7
(t, logt, abund) = get_tgrid(26, 6, amin=1E-10, limits=[5,7.5])
nt = len(t)
for k in range(3,12):
    y = []
    for i in range(nt):
        y.append(abund[i][k-1]/abund[i][k])
    p.add(biggles.Curve(logt, y, linetype='solid', linewidth=2))

p.xrange = (6.3, 7.8)
p.yrange = (0.05, 200)
for k in [3,5,7,9,11]:
    (tmax,amax) = MaxAbund(26, k)
    tmax = log10(tmax/const.kb)
    x = [tmax, tmax]
    y = [0.1, 120]
    p.add(biggles.Curve(x, y, linetype = 'dashed', linewidth=2))
    p.add(biggles.Label(tmax, 140, '%d+'%(26-k), size=0.5))
    
    if (k == 11):
        (tmax,amax) = MaxAbund(26, k-1)
        tmax = log10(tmax/const.kb)
        x = [tmax, tmax]
        y = [0.1, 120]
        p.add(biggles.Curve(x, y, linetype = 'dashed', linewidth=2))
        p.add(biggles.Label(tmax, 140, '%d+'%(27-k), size=0.5))
        
    p.add(biggles.Label(7.65, abund[-1][k-1]/abund[-1][k],
                        'Fe$^{%d+}$/Fe$^{%d+}$'%(27-k,26-k), size=1.0))
    
p.write_eps('abund.eps')


p = biggles.FramedPlot()
p.ylog = 1
p.ylabel = 'Abundance Ratio'
p.xlabel = 'Log Temperature (K)'
p.aspect_ratio = 0.7

tbl = TABLE(fname='trates.tbl')
tbl.open('r')
tbl.read_header()
d = tbl.read_columns([0,1,2,4,6,8])
tbl.close()
rr = []
t1 = []
i = 0
while (d[0][i] == 2):
    t1.append(d[1][i])
    rr.append(d[2][i]+d[3][i])
    i = i + 1

ymax = 0
ymin = 1E30
for k in range(3,12):
    t2 = []
    ir = []
    rrp = []
    while (d[0][i] == k):
        t2.append(d[1][i])
        ir.append(d[4][i]+d[5][i])
        rrp.append(d[2][i]+d[3][i])
        i = i + 1
        if (i == len(d[0])):
            break
    t = []
    ir0 = []
    rr0 = []
    nt1 = len(t1)
    nt2 = len(t2)
    q1 = 0
    q2 = 0
    while (q1 < nt1 and q2 < nt2):
        if (fabs(t1[q1] - t2[q2]) < 1E-5):
            t.append(t1[q1])
            ir0.append(ir[q2])
            rr0.append(rr[q1])
            q1 = q1 + 1
            q2 = q2 + 1
        elif (t1[q1] < t2[q2]):
            q1 = q1 + 1
        else:
            q2 = q2 + 1
    y = map(lambda a,b:a/b, ir0, rr0)
    p.add(biggles.Curve(t, y, linetype='solid', linewidth=2))
    if k in [3, 5, 7, 9, 11]:
        p.add(biggles.Label(t[-1], y[-1],
                        'Fe$^{%d+}$/Fe$^{%d+}$'%(27-k,26-k), size=1.0))
    t1 = t2[:]
    rr = rrp[:]
p.yrange = (0.05, 200)
p.xrange = (6.0, 8.2)
p.write_eps('abund0.eps')   
