""" calculate the electron impact ionization cross sections
"""

# import the modules
from pfac import fac

fac.SetAtom('Fe')
# 1s shell is closed
fac.Closed('1s')
# Ne-like ground state
fac.Config('2*8', group='fe17')
# F-like configuations
fac.Config('2*7', group='fe18')

# solve the structure problem
fac.ConfigEnergy(0)
fac.OptimizeRadial(['fe17'])
fac.ConfigEnergy(1)
fac.Structure('ne_f.lev.b', ['fe17'])
fac.Structure('ne_f.lev.b', ['fe18'])
fac.MemENTable('ne_f.lev.b')
fac.PrintTable('ne_f.lev.b', 'ne_f.lev', 1)

# set the output collision energies
e = [500.0, 900.0, 1.3e3, 1.7e3, 2.1e3, 4.2e3, 6.0e3, 8.0e3]
fac.SetUsrCIEGrid(e)
fac.CITable('ne.ci.b', ['fe17'], ['fe18'])
fac.PrintTable('ne.ci.b', 'ne.ci', 1)
