from pfac.spm import *
from pfac import const
import os

z = 26
neles = range(3,4)

dir0 = 'data/'
den = [1.0]

for nele in neles:
    dir1 = 'Fe%02d/'%nele
    try:
        os.system('mkdir %s'%dir1)
    except:
        pass
    (temp, logt, population) = get_tgrid(z, nele)
    nt = len(temp)
    nd = len(den)
    print 'NION = 3'
    spectrum([nele], temp, den, population,
             'Fe', dir0=dir0, dir1=dir1, nion=3)
    print 'NION = 2'
    spectrum([nele], temp, den, population,
             'Fe', dir0=dir0, dir1=dir1, nion=2)
    print 'NION = 1'
    spectrum([nele], temp, den, population,
             'Fe', dir0=dir0, dir1=dir1, nion=1)

    
