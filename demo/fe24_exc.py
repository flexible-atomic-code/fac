"""
calculation of collisional strength. 
"""

import config
import fac

fac.SetAtom('Fe',26)
config.closed('1s')
config.config('2*1', group = 'n2')
config.config('3*1', group = 'n3')
fac.OptimizeRadial('n2', 'n3')
fac.Structure()

fac.LevelTable('li.lev')

fac.SetUsrCEGrid([1e2, 1e3, 2e3, 5e3])
fac.CETable(['n2'], ['n3'], 'li.ce')
