"""calculate the electron impact excitation cross sections
"""

# import the modules
from pfac import fac


def test_li_excitation():
    fac.SetAtom('Li')
    # 1s shell is closed
    fac.Closed('1s')
    fac.Config('2*1', group='n2')

    # Self-consistent iteration for optimized central potential
    fac.ConfigEnergy(0)
    fac.OptimizeRadial('n2')
    fac.ConfigEnergy(1)
    fac.Structure('ne.lev.b')
    fac.MemENTable('ne.lev.b')
    fac.PrintTable('ne.lev.b', 'ne.lev', 1)

    fac.CETable('ne.ce.b', ['n2'])
    fac.PrintTable('ne.ce.b', 'ne.ce', 1)
