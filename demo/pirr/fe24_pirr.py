""" calculate the photoionization and 
    radiative recombination cross sections
"""

# import the modules
from pfac import fac 

fac.SetAtom('Fe')

# specify the configurations for both recombining
# and recombined ions.
fac.Config('1s2', group = 'n1')
fac.Config('1s1 2*1', group = 'n2')
fac.Config('1s2 2*1', group = 'rn2')

# since the recombined electron is in n=2 shell,
# set the appropriate screening
fac.SetScreening([2])

fac.OptimizeRadial('n1')

# configuration interaction between n=1 and n=2
# complexes are included for the recombining ion.
fac.Structure('li.lev.b', ['n1', 'n2'])
fac.Structure('li.lev.b', ['rn2'])
fac.MemENTable('li.lev.b')
fac.PrintTable('li.lev.b', 'li.lev', 1)

fac.RRTable('li.rr.b', ['rn2'], ['n1'])
fac.PrintTable('li.rr.b', 'li.rr', 1)

