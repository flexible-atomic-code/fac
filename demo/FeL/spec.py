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
    spectrum([nele], den, population,
             'Fe', dir0=dir0, dir1=dir1, nion=3, dist=[0,temp,-1,-1])
    rates3 = read_rates(nt, nd, nele, dir=dir1, nion=3)
    save_rates(rates3, dir1+'rates3.sav', dir1+'rates3.dat',
               temp=temp, logt=logt, abund=population)
    
    print 'NION = 2'
    spectrum([nele], den, population,
             'Fe', dir0=dir0, dir1=dir1, nion=2, dist=[0,temp,-1,-1])
    rates2 = read_rates(nt, nd, nele, dir=dir1, nion=2)
    save_rates(rates2, dir1+'rates2.sav', dir1+'rates2.dat',
               temp=temp, logt=logt, abund=population)
    
    print 'NION = 1'
    spectrum([nele], den, population,
             'Fe', dir0=dir0, dir1=dir1, nion=1, dist=[0,temp,-1,-1])
    rates1 = read_rates(nt, nd, nele, dir=dir1, nion=1)
    save_rates(rates1, dir1+'rates1.sav', dir1+'rates1.dat',
               temp=temp, logt=logt, abund=population)
    
