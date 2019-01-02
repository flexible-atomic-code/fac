"""
This script generates atomic data for neutral oxygen without any CI.
"""
from pfac.fac import *

SetAtom('O')
Closed('1s')
Config('g1', '2*6')
Config('g3', '2*5 3*1')
Config('g4', '2*5 4*1')
Config('ion', '2*5')
ConfigEnergy(0)
OptimizeRadial('g1')
ConfigEnergy(1)
SetCILevel(-1)
Structure('diagh_O_b.en', ['g1', 'g3', 'g4', 'ion'])
PrintTable('diagh_O.en', 'diagh_O_a.en')
