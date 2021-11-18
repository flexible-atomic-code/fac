"""
s.py runs the collisional radiative module. sel.py prints the line intensities
in ascii format.
"""
import sys
from pfac.crm import *
from pfac import fac
from pfac import spm

use_openmp = False
if len(sys.argv) == 2 and sys.argv[1] == 'openmp':
    use_openmp = True


if use_openmp:
    # enable openmp with 2 cores
    InitializeMPI(2)

# atomic number
z = 10
# number of electrons
nele = 3

temp = 500.0*((1+z-nele)/16.0)**2
dens = 1.0

# a is the atomic symbol corresponding to z
a = fac.ATOMICSYMBOL[z]

p = '%s%02d'%(a, nele)

# add the ion with nele number of electrons
AddIon(nele, 0.0, '%sb'%p)
SetBlocks(-1)

SetEleDist(0, temp, -1, -1)
# radiative transition rates
SetTRRates(0)
# collisional excitation rate coefficients.
SetCERates(1)

# abundance of the ion with nele number of electrons
SetAbund(nele, 1.0)
SetEleDensity(dens)

InitBlocks()

SetIteration(1e-6, 0.5)

LevelPopulation()

# spectral data table
SpecTable(p+'b.sp', 0)

PrintTable(p+'b.sp', p+'a.sp')

# save rate matrix
rate = spm.get_complexes(nele)
RateTable(p+'b.rt', rate)
PrintTable(p+'b.rt', p+'a.rt', 1)

if use_openmp:
    FinalizeMPI()
