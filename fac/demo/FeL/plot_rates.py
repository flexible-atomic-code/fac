from pfac.table import *
from pfac.spm import *
from pfac.crm import *
import pprint
import biggles

tbl = TABLE('rates.tbl')
tbl.open('r')
tbl.read_header()
shells = ['3s', '3p', '3d']
pts = ['filled diamond', 'star', 'filled circle']
psize = [2.0, 2.0, 1.5]
labels = ['', 'RE', 'DR+RR', 'CI']
lts = ['', 'solid', 'dashed', 'dotted']

xkey = [0.6, 0.6, 0.6]
ykey = [0.6, 0.6, 0.6]
biggles.configure('fontsize_min', 1.0)

for k in range(3, 11):
    print 'NELE = %d'%k
    efile = 'data/Fe%02db.en'%k
    c = get_complexes(k)
    i = 0
    ish = [[],[],[]]
    while (1):
        s0 = LevelInfor(efile, i)
        if (not s0[3].strip() in c):
            break
        s = s0[4].strip()
        s = s.split()
        if (s[-1][:2] == '3s'):
            ish[0].append(i)
        elif (s[-1][:2] == '3p'):
            ish[1].append(i)
        elif (s[-1][:2] == '3d'):
            ish[2].append(i)
        i = i + 1

    d = tbl.read_columns([1,2,8,9,11,12], filter='c[0]==%d'%k)
    tbl.rewind()
    nt = 0
    while (d[0][nt] == 0):
        nt = nt + 1
    x = d[1][:nt]
    n = 3
    if (k == 10):
        n = 4
    rts = []
    for i in range(3):
        a = []
        for j in range(n):
            a.append([0]*nt)
        rts.append(a)
        
    i = 0
    nd = len(d[0])
    while (i < nd):
        for q in range(3):
            if (d[0][i] in ish[q]):
                for t in range(nt):
                    m = i + t
                    for j in range(n):
                        rts[q][j][t] = rts[q][j][t] + d[j+2][m]
        i = i + nt
        
    p = biggles.FramedPlot()
    p.aspect_ratio = 0.7
    p.title = 'Fe$^{%d+}$'%(26-k)
    p.ylog = 1
    p.xlabel = 'Log Temperature (K)'
    p.ylabel = 'Ratio to CE'
    bp = [0]*(n+2)
    ymax = 0
    ymin = 1E10
    for i in range(3):
        for j in range(1, n):
            if (k == 10 and i > 0 and j == n-1):
                continue
            y = map(lambda a,b:a/b, rts[i][j], rts[i][0])
            b0 = biggles.Curve(x, y, linetype=lts[j],
                               linewidth=2)
            b1 = biggles.Points(x, y, symboltype=pts[i],
                                   symbolsize=psize[i])
            if (i == 0):
                b0.label = labels[j]
                bp[2+j] = b0
            if (j == 1):
                b1.label = shells[i]
                bp[i] = b1
            p.add(b0)
            p.add(b1)
            if (k == 10 and j == n-1):
                continue
            ym = max(y)
            if (ymax < ym):
                ymax = ym
            ym = min(y)
            if (ymin > ym):
                ymin = ym
    p.yrange = (ymin*0.85, ymax*1.25)
    x0 = 0.15
    if (k < 10):
        y0 = 0.3
    else:
        y0 = 0.25
    p.add(biggles.PlotKey(x0, y0, bp[:3]))
    x0 = 0.7
    y0 = 0.85
    p.add(biggles.PlotKey(x0, y0, bp[3:]))
    p.write_eps('rates%02d.eps'%k)
        
            
