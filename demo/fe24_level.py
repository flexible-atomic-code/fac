"""
calculation of the energy levels. and radiative transition rates
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
fac.TransitionTable(['n2'], ['n3'], 'li.tr', -1)

