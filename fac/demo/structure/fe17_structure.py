
""" calculation of the energy levels. and radiative transition rates
"""

from pfac import fac

fac.SetAtom('Fe')
# 1s shell is closed
fac.Closed('1s')
fac.Config('2*8', group = 'n2')
fac.Config('2*7 3*1', group = 'n3')

# Self-consistent iteration for optimized central potential
fac.OptimizeRadial('n2', 'n3')
fac.Structure('ne.lev.b')
fac.PrintTable('ne.lev.b', 'ne.lev')

fac.TransitionTable('ne.tr.b', ['n2'], ['n3'])
fac.PrintTable('ne.tr.b', 'ne.tr')

