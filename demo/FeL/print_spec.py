from pfac.crm import *
from pfac.spm import *

try:
    import cPickle
except ImportError:  # Python3 does not have cPickle
    import _pickle as cPickle


de = [1.0, 10.0]
elimits = [[],[],[],
           [1.e3, 1.8e3], #03
           [950., 1.8e3], #04
           [850., 1.7e3], #05
           [800., 1.6e3], #06
           [750., 1.5e3], #07
           [700., 1.4e3], #08
           [600., 1.3e3], #09
           [600., 1.3e3], #10
           [600., 1.3e3]] #11

type = 2000100

for k in range(3, 11):
    pref = 'Fe%02d/Fe%02d'%(k, k)
    rates1 = cPickle.load(open('Fe%02d/rates1.sav'%k))
    temp = rates1['temp']
    nt = len(temp)
    emin = elimits[k][0]
    emax = elimits[k][1]
    for t in range(nt):
        for j in range(1, 4):
            spf = '%sb_t%02dd0i%d.sp'%(pref, t, j)
            for i in range(len(de)):
                pltf = '%sa_t%02di%dw%d.plt'%(pref, t, j, i)
                print('NELE=%02d T=%d NION=%d DE=%d %s'%(k,t,j,i,pltf))
                PlotSpec(spf, pltf, k, type, emin, emax, de[i])
                if (k == 10 and j == 3):
                    pltf = pltf+'+1'
                    PlotSpec(spf, pltf, k+1, type, emin, emax, de[i])
