""" calculate the autoionization rates for Ne-like Se.
"""

# import the modules
from pfac import config, fac

fac.SetAtom('Se')

# configurations for the F-like ion
config.closed('1s')
config.closed('2s')
config.config('2p5', group = 'n2')

# configurations of doubly excited Ne-like ion
config.config('2p4 3s2', '2p4 3s1 3p1', group = 'n33')

fac.OptimizeRadial('n33')
fac.Structure('n2')
fac.Structure('n33')

fac.LevelTable('se.lev')
fac.AITable(['n33'], ['n2'], 'se.ai')

