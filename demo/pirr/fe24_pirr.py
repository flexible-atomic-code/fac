""" calculate the photoionization and 
    radiative recombination cross sections
"""

# import the modules
from pfac import config, fac 

fac.SetAtom('Fe')

# specify the configurations for both recombining
# and recombined ions.
config.config('1s2', group = 'n1')
config.config('1s1 2*1', group = 'n2')
config.config('1s2 2*1', group = 'rn2')

# since the recombined electron is in n=2 shell,
# set the appropriate screening
fac.SetScreening([2])

fac.OptimizeRadial('n1')

# configuration interaction between n=1 and n=2
# complexes are included for the recombining ion.
fac.Structure(['n1', 'n2'])
fac.Structure('rn2')

fac.LevelTable('li.lev')

# photo-electron energy grid
fac.SetUsrPEGrid([20.0, 400.0, 2E3, 4E3, 6E3, 10E3])

fac.RRTable(['rn2'], ['n1'], 'li.rr')

