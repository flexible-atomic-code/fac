"""
calculation of photoionization and radiative recombination. 
"""

import config
import fac

fac.SetAtom('Fe',26)
config.config('1s2 2*1', group = 'n2')
#this is the Be-like configuration after recombination
config.config('1s2 2*2', group = 'be_n2')

fac.OptimizeRadial('n2')
fac.Structure('n2')
fac.Structrue('be_n2')

fac.LevelTable('li.lev')

fac.SetPEGrid([1e2, 5e2, 1e3, 2e3])
fac.RRTable(['be_n2'], ['n2'], 'li.rr')
