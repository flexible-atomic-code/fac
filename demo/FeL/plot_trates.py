from pfac.table import *
from pfac.spm import *
import pprint
import biggles
from math import *
from pfac import const

tbl = TABLE(fname='trates.tbl')
tbl.open('r')
tbl.read_header()
biggles.configure('fontsize_min', 1.0)
for k in range(3,11):
    print 'NELE %d'%k
    d = tbl.read_columns(range(1,10), filter='c[0] == %d'%k)
    tbl.rewind()
    p = biggles.FramedPlot()
    p.title = 'Fe$^{%d+}$'%(26-k)
    p.xlabel = 'Log Temperature (K)'
    p.ylabel = 'Rate Coeff. (10$^{-10}$ cm$^{3}$ s$^{-1}$)'
    p.ylog = 1
    p.aspect_ratio = 0.7
    ymax = 0
    ymin = 1E10
    n = 7
    if (k == 11):
        n = 9
    b = [0,0]
    for i in range(1,n,2):
        if (k == 11 and i < 5):
            continue
        c = biggles.Curve(d[0], d[i], linetype='solid', linewidth=2)
        if (i == 5):
            c.label = 'Present'
            b[0] = c
        p.add(c)
        ym = max(d[i])
        if (ym > ymax):
            ymax = ym
        ym = min(d[i])
        if (ym < ymin):
            ymin = ym
    for i in range(2, n, 2):
        c = biggles.Curve(d[0], d[i], linetype='dashed', linewidth=2)
        if (i == 2):
            c.label = 'AR92'
            b[1] = c
        p.add(c)
        ym = max(d[i])
        if (ym > ymax):
            ymax = ym
        ym = min(d[i])
        if (ym < ymin):
            ymin = ym
    x = d[0][-1] + 0.05
    p.xrange = (min(d[0])-0.1, x+0.1)
    p.yrange = (ymin*0.75, ymax*1.25)
    p.add(biggles.Label(x, d[2][-1], 'DR'))
    p.add(biggles.Label(x, d[4][-1], 'RR'))
    p.add(biggles.Label(x, d[6][-1]*0.9, 'CI'))
    if (k == 11):
        p.add(biggles.Label(x, d[7][-1], 'EA'))

    if (k == 5):
        x = [2e2,3e2,4e2,5e2,6e2,7e2,8e2,9e2,1e3,1.4e3,
             2e3,2.4e3,3e3,3.4e3,4e3,4.4e3,5e3]
        y = [5.208e-1,4.228e-1,3.719e-1,3.33e-1,2.998e-1,2.708e-1,
             2.456e-1,2.237e-1,2.043e-1,1.483e-1,1.009e-1,8.168e-2,
             6.242e-2,5.345e-2,4.35e-2,3.848e-2,3.256e-2]
        x = map(lambda a:log10(a/const.kb),x)
        c = biggles.Points(x, y, symboltype='filled circle',symbolsize=1.5)
        p.add(c)
        c.label = 'Chen98'
        b.append(c)
    if (k == 9):
        x = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0,4.0,5.0,6.0]
        y = [2e-1,2.71e-1,2.87e-1,2.78e-1,2.6e-1,2.4e-1,2.2e-1,
             2.01e-1,1.84e-1,9.02e-2,5.48e-2,3.76e-2,2.78e-2,2.16e-2]
        y1 = [9.86e-2,6.41e-2,4.56e-2,3.44e-2,2.71e-2,2.21e-2,1.84e-2,
              1.57e-2,1.36e-2,5.06e-3,2.8e-3,1.84e-3,1.32e-3,1.01e-3]
        y = map(lambda a,b:a+b, y, y1)
        x = map(lambda a:log10(a*1e3/const.kb), x)
        c = biggles.Points(x[:-5], y[:-5],
                           symboltype='filled circle',symbolsize=1.5)
        p.add(c)
        c.label = 'Roszman87a'
        b.append(c)
    if (k == 8):
        x = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0,4.0,5.0]
        y = [1.1e-1,1.84e-1,2.14e-1,2.19e-1,2.11e-1,1.99e-1,1.86e-1,
             1.72e-1,1.60e-1,8.15e-2,5.02e-2,3.47e-2,2.58e-2]
        y1 = [1.52e-1,9.57e-2,6.7e-2,5.01e-2,3.93e-2,3.19e-2,2.65e-2,
             2.25e-2,1.94e-2,7.18e-3,3.97e-3,2.60e-3,1.87e-3]
        y = map(lambda a,b:a+b, y, y1)
        x = map(lambda a:log10(a*1e3/const.kb), x)
        c = biggles.Points(x[:-3], y[:-3],
                           symboltype='filled circle',symbolsize=1.5)
        p.add(c)
        c.label = 'Roszman87b'
        b.append(c)
        
    p.add(biggles.PlotKey(0.3, 0.2, b))
    
    p.write_eps('trates%02d.eps'%(k))
        
tbl.close()
