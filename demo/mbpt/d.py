#an exmaple of calculating He-like energy levels in mbpt
from pfac.fac import *
import sys

# job index for split n1 array
i = int(sys.argv[1])
# number of n1 splits
n = int(sys.argv[2])
# number of processors in parallel mode
np = int(sys.argv[3])
if i < 0:
    np = 0
    
if (np > 1):
    InitializeMPI(np)

if (n == 1):
    n1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 24, 28, 32, 38, 43, 50, 65, 80, 100, 125]
else:
    if ( i == 0 ):
        n1 = [1, 2, 3, 4, 5, 6, 7]
    elif ( i == 1 ):
        n1 = [8, 9, 10, 11]
    elif ( i == 2 ):
        n1 = [12, 13, 15, 17, 20]
    else:
        n1 = [24, 28, 32, 38, 43, 50, 65, 80, 100, 125]

n2 = range(8)+[8, 13, 20, 28, 38, 50, 65]

SetAtom('Fe')

SetVP(103)
SetMS(1,1)
SetSE(-1,161)
SetBreit(-1,1,5,1e-5)

Config('g0', '1s2')
Config('g1', '1*1 2*1')
Config('g2', '1*1 3*1')
g = ['g0','g1','g2']
ng = len(g)
if i>=0:
    OptimizeRadial('g0')
    SetBoundary(3, 1e-5, 1e30)
    ReinitRadial(0)
    SetRadialGrid(3000, 1.1, -1e30, 0)
OptimizeRadial('g0')
if i>=0:
    SetBoundary(3, 1e-5, 1e30)

GetPotential('pot.txt')

TransitionMBPT(3, 3)
if i>=0:
    StructureMBPT('tb%d.en'%i, 'tb%d.ham'%i, g,
                  [-1,-1,-1,-1,-1,-1,-1,-1,12,12,12,12],
                  n1, n2, ng)
    MemENTable('tb%d.en'%i)
    PrintTable('tb%d.en'%i, 'ta%d.en'%i)
else:
    h = ['tb%d.ham'%i for i in range(n)]
    TransitionMBPT('tb.tr', ['g0'], ['g1','g2'])
    StructureMBPT('tb.en', 'ta.ham', h, g, ng)
    MemENTable('tb.en')
    PrintTable('tb.en', 'ta.en')
    PrintTable('tb.tr', 'ta.tr')

if np > 1:
    FinalizeMPI()
    
