""" calculate the electron impact ionization cross sections
"""

# import the modules
import fac
import config

fac.SetAtom('Fe')
# 1s shell is closed
config.closed('1s')
# Ne-like ground state
config.config('2*8', group='fe17')
# F-like configuations
config.config('2*7', group='fe18')

# solve the structure problem
fac.OptimizeRadial(['fe17', 'fe18'])
fac.Structure('fe17')
fac.Structure('fe18')
fac.LevelTable('ne_f.lev')

# set the collision energies
# 21 points with 400 eV step starting from 500 eV.
e = [500.0]
for i in range(20):
    e.append(e[i]+400.0)
fac.SetUsrCIEGrid(e)

fac.CITable(['fe17'], ['fe18'], 'ne.ci')

