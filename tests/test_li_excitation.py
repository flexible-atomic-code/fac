"""calculate the electron impact excitation cross sections
"""
import os
from pfac import fac


THIS_DIR = os.path.abspath(os.path.join(__file__, os.pardir))
output_dir = THIS_DIR + '/data/'


def test_li_excitation():
    fac.SetAtom('Li')
    # 1s shell is closed
    fac.Closed('1s')
    fac.Config('2*1', group='n2')

    # Self-consistent iteration for optimized central potential
    fac.ConfigEnergy(0)
    fac.OptimizeRadial('n2')
    fac.ConfigEnergy(1)
    fac.Structure(output_dir + 'ne.lev.b')
    fac.MemENTable(output_dir + 'ne.lev.b')
    fac.PrintTable(output_dir + 'ne.lev.b', 'ne.lev', 1)

    fac.CETable(output_dir + 'ne.ce.b', ['n2'])
    fac.PrintTable(output_dir + 'ne.ce.b', 'ne.ce', 1)
