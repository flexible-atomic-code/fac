"""
calculation of ionization cross sections. 
"""

import config
import fac

fac.SetAtom('Fe',26)
#this is He-like configulation after ionization 
config.config('1s2', group = 'n1')
#this is the Li-like configuration
config.config('1s2 2*1', group = 'n2')

#the radial optimzation is carried out in terms of ionized configurations
fac.OptimizeRadial('n1')
fac.Structure('n1')
fac.Structrue('n2')

fac.LevelTable('li.lev')

fac.SetUsrCIEGrid([1e2, 1e3, 2e3])
fac.CITable(['n2'], ['n1'], 'li.ci')
