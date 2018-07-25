"""
This script generates atomic data for Z=10 neon, number of electrons =3,
maximum n=5 for excitation.
https://www-amdis.iaea.org/FAC/
"""
import sys
from pfac.fac import *

use_openmp = False
if len(sys.argv) == 2 and sys.argv[1] == 'openmp':
    use_openmp = True


if use_openmp:
    # enable openmp with 2 cores
    InitializeMPI(2)

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

if use_openmp:
    FinalizeMPI()
