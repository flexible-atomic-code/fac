import sys
from pfac.fac import *

# atomic number, number of electrons, number of excitation levels
z = 10
nele = 3
nmax = 5

a = ATOMICSYMBOL[z]  # atomic symbol (Ne)
p = '%s%02d'%(a, nele)

SetAtom(a)

Config('g2', '1s2 2*%d'%(nele-2))
for n in range(3, nmax+1):
    Config('g%d'%n, '1s2 2*%d %d*1'%(nele-3, n))

# radial potential
ConfigEnergy(0)
OptimizeRadial('g2')
ConfigEnergy(1)

