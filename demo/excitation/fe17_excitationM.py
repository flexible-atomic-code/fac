"""calculate the electron impact excitation cross sections
"""

# import the modules
from pfac import fac

fac.SetAtom('Fe')
# 1s shell is closed
fac.Closed('1s')
fac.Config('2*8', group = 'n2')
fac.Config('2*7 3*1', group = 'n3')

# Self-consistent iteration for optimized central potential
fac.SetScreening([10])
fac.ConfigEnergy(0)
fac.OptimizeRadial('n2', 'n3')
fac.ConfigEnergy(1)
fac.Structure('ne.lev.b')
fac.MemENTable('ne.lev.b')
fac.PrintTable('ne.lev.b', 'ne.lev', 1)

fac.CETableMSub('ne.ceM.b', [0], [2,4,16,22,26,30,32])
fac.PrintTable('ne.ceM.b', 'ne.ceM', 1)

