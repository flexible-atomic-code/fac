from pfac.fac import *

"""
Calculate R-matrix basis and surface amplitudes for He-like O.
Targets include 1s2, 1s1 2*1, and 1s1 3*1 configurations.
"""

#ConvertToSFAC('t.sf')
#InitializeMPI(8)

SetAtom('O')

# Target configurations
Config('g0', '1s2')
Config('g1', '1s1 2*1')
Config('g2', '1s1 3*1')

# Correlation configurations
Config('g3', '1s2 2*1')
Config('g3', '1s2 3*1')
Config('g3', '1s1 2*2')
Config('g3', '1s1 2*1 3*1')
Config('g3', '1s1 3*2')

# disable Breit interaction, and nuclear recoil effects.
SetBreit(0)
SetMS(0,0)
SetSE(0)
SetVP(0)

# use maximum number of grid points, important.
SetRadialGrid(3000,1.4,300,0)

OptimizeRadial('g0')

Structure('tb.en', ['g0','g1','g2'])
Structure('tb.en', ['g3'])

MemENTable('tb.en')
PrintTable('tb.en', 'ta.en')
BasisTable('bas.d')

# specify targets correlation states.
RMatrixTargets(['g0','g1','g2'], ['g3'])

nb = 25
kmax = 30
# slater integrals with l >= 12 are treated with no exchange.
SetSlaterCut(12,12)

# boundary radius of the first R-matrix zone.
# it enclose the 1s,2s,2p,3s,3p,3d wavefunctions
# with amplitude > 1E-6
RMatrixBoundary(0, 1e-6, 0.0)
GetPotential('pot.d')

# calculate basis functions, with maximum angular momemtum kmax,
# and nb functions per-kappa.
RMatrixBasis('rmx.b0', kmax, nb)

# R-matrix surface amplitudes.
RMatrixSurface('rmx.d0')

# boundary radius of the second zone.
# it starts from the end of first zone, to the radius 13.0 a.u.
RMatrixBoundary(-1, 13.0, 0.0)

# calculate basis functions in the second zone.
RMatrixBasis('rmx.b1', kmax, nb)
# R-matrix surface amplitudes.
RMatrixSurface('rmx.d1')

#FinalizeMPI()

#CloseSFAC()

