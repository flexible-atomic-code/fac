""" calculate the autoionization rates for Ne-like Se.
"""

# import the modules
from pfac import fac

fac.SetAtom('Se')

# configurations for the F-like ion
fac.Closed('1s')
fac.Closed('2s')
fac.Config('2p5', group = 'n2')

# configurations of doubly excited Ne-like ion
fac.Config('2p4 3s2', '2p4 3s1 3p1', group = 'n33')

fac.OptimizeRadial('n33')
fac.Structure('se.lev.b', ['n2'])
fac.Structure('se.lev.b', ['n33'])
fac.PrintTable('se.lev.b', 'se.lev')

fac.AITable('se.ai.b', ['n33'], ['n2'])
fac.PrintTable('se.ai.b', 'se.ai')


