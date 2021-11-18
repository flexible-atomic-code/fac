from pfac.fac import *

SetAtom('Fe')

Config('g0', '1s2')
Config('g1', '1s1 2*1')

ConfigEnergy(0)
OptimizeRadial('g0')
ConfigEnergy(1)

Structure('tb.en', ['g0','g1'])
MemENTable('tb.en')
PrintTable('tb.en', 'ta.en', 1)

for m in [-1,1,-2,2]:
    TransitionTable('tb.tr', ['g0'], ['g1'], m)
    TransitionTable('tb.tr', ['g1'], ['g1'], m)
PrintTable('tb.tr', 'ta.tr', 1)

CETableMSub('tb.ce', ['g0'], ['g1'])
CETableMSub('tb.ce', ['g1'], ['g1'])
PrintTable('tb.ce', 'ta.ce', 1)
