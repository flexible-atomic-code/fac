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
fac.OptimizeRadial(['fe17', 'fe18'])
fac.Structure('ne_f.lev.b', ['fe17'])
fac.Structure('ne_f.lev.b', ['fe18'])
fac.PrintTable('ne_f.lev.b', 'ne_f.lev')

# set the collision energies
# 21 points with 400 eV step starting from 500 eV.
e = [500.0]
for i in range(20):
    e.append(e[i]+400.0)
fac.SetUsrCIEGrid(e)

fac.CITable('ne.ci.b', ['fe17'], ['fe18'])
fac.PrintTable('ne.ci.b', 'ne.ci')
