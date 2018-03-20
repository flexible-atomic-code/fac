from pfac.spm import *

z = 26
neles = range(3,11)

for nele in neles:
    dir1 = 'Fe%02d/'%nele
    print('%d %s'%(nele, dir1))
    (temp, logt, population) = get_tgrid(z, nele)
    nt = len(temp)
    nd = 1

    rates3 = read_rates(nt, nd, nele, dir=dir1, nion=3)
    save_rates(rates3, dir1+'rates3.sav', dir1+'rates3.dat',
               temp=temp, logt=logt, abund=population)

    rates2 = read_rates(nt, nd, nele, dir=dir1, nion=2)
    save_rates(rates2, dir1+'rates2.sav', dir1+'rates2.dat',
               temp=temp, logt=logt, abund=population)

    rates1 = read_rates(nt, nd, nele, dir=dir1, nion=1)
    save_rates(rates1, dir1+'rates1.sav', dir1+'rates1.dat',
               temp=temp, logt=logt, abund=population)
