from pfac.spm import *
from pfac import const
import os

z = 26
neles = range(1,3)

dir0 = 'data/'
den = [1.0]

for nele in neles:
    dir1 = 'Li%02d/'%nele
    try:
        os.system('mkdir %s'%dir1)
    except:
        pass
    (temp, logt, population) = get_tgrid(z, nele)
    nt = len(temp)
    nd = len(den)
    print 'NION = 2'
    spectrum([nele], temp[:2], den, population,
             'Li', dir0=dir0, dir1=dir1, nion=2)
    print 'NION = 1'
    spectrum([nele], temp[:2], den, population,
             'Li', dir0=dir0, dir1=dir1, nion=1)
