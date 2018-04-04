"""
This script generates atomic data for Z=10 neon, number of electrons =3,
maximum n=5 for excitation.
https://www-amdis.iaea.org/FAC/
"""
from pfac.fac import *

# atomic number
z = 10
# number of electrons
nele = 3
# number of excitation levels to be considered
nmax = 5

# a is the atomic symbol corresponding to z, i.e. 'Ne'
a = ATOMICSYMBOL[z]
p = '%s%02d'%(a, nele)

SetAtom(a)

Config('g2', '1s2 2*%d'%(nele-2))
for n in range(3, nmax+1):
    Config('g%d'%n, '1s2 2*%d %d*1'%(nele-3, n))

# radial potential
ConfigEnergy(0)
OptimizeRadial('g2')
ConfigEnergy(1)

# modify the energy levels slightly
##CorrectEnergy([0,1,2,3,4],[0.0,722.8,724.8,735.4,753.1],1)

# atomic structure
Structure(p+'b.en', ['g2', 'g3'])

for n in range(4, nmax+1):
    Structure(p+'b.en', ['g%d'%n])

MemENTable(p+'b.en')
PrintTable(p+'b.en', p+'a.en')

# transition rates
for n in range(2, nmax+1):
    for m in range(n, nmax+1):
        Print('TR: g%d -> g%d'%(n, m))
        TRTable(p+'b.tr', ['g%d'%n], ['g%d'%m])
PrintTable(p+'b.tr', p+'a.tr')

# excitation
for n in range(2, nmax+1):
    if n > 2:
        continue
    for m in range(n, nmax+1):
        Print('CE: g%d -> g%d'%(n, m))
        CETable(p+'b.ce', ['g%d'%n], ['g%d'%m])
PrintTable(p+'b.ce', p+'a.ce')
