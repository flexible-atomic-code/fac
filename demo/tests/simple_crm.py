import os
import sys
from pfac.fac import *
from pfac.crm import *

print('running crm_OI.py')
if not os.path.exists('O_crm'):
    os.mkdir('O_crm')

use_openmp = False
if len(sys.argv) == 2 and sys.argv[1] == 'openmp':
    use_openmp = True

if use_openmp:
    # enable openmp with 2 cores
    InitializeMPI(2)

# For Structure and rates
SetAtom('O')
Closed('1s')
Config('g1', '2s2 2p4')
Config('g2', '2s2 2p3 3s1')
# Config('ion', '2s2 2p3')

ConfigEnergy(0)
OptimizeRadial('g1')
ConfigEnergy(1)

Structure('O_crm/O.en', ['g1', 'g2'])
# Structure('O_crm/O.en', ['ion'])
MemENTable('O_crm/O.en')
TransitionTable('O_crm/O.tr', ['g1', 'g2'], ['g1', 'g2'])
CETable('O_crm/O.ce', ['g1', 'g2'], ['g1', 'g2'])

# For CRmodel
AddIon(8, 1, 'O_crm/O')
SetBlocks(-1)
SetEleDist(0, 1.0, -1, -1)
SetTRRates(0)
SetCERates(1)
SetEleDensity(1.0)
InitBlocks()
SetIteration(1e-6, 0.5)
LevelPopulation()
SpecTable('O.b.sp', 0)
PrintTable('O.b.sp', 'O.sp')

if use_openmp:
    FinalizeMPI()
