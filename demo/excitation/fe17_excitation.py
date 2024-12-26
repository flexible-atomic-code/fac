"""calculate the electron impact excitation cross sections
"""
import sys
from pfac import fac


use_openmp = False
if len(sys.argv) == 2 and sys.argv[1] == 'openmp':
    use_openmp = True

if use_openmp:
    # enable openmp with 2 cores
    fac.InitializeMPI(2)

fac.SetAtom('Fe')
# 1s shell is closed
fac.Closed('1s')
fac.Config('2*8', group = 'n2')
fac.Config('2*7 3*1', group = 'n3')

# Self-consistent iteration for optimized central potential
fac.ConfigEnergy(0)
fac.OptimizeRadial('n2')
fac.ConfigEnergy(1)
fac.Structure('ne.lev.b')
fac.MemENTable('ne.lev.b')
fac.PrintTable('ne.lev.b', 'ne.lev', 1)

fac.CETable('ne.ce.b', ['n2'], ['n3'])
fac.PrintTable('ne.ce.b', 'ne.ce', 1)

if use_openmp:
    fac.FinalizeMPI()
