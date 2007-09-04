from pfac.pol import *

# setup the magnetic sublevel table and transition rates.
SetMLevels('tb.en', 'tb.tr')
# setup the magnetic sublevel excitation rates at E=12.24 keV (900 Ryd)
SetEnergy(8E3, 12.0)
SetMCERates('tb.ce')
# calculates the magnetic sublevel populations.
# results are output in pop.txt.
SetDensity(1.0)
PopulationTable('pop.txt')
# calculate the orientation parameters.
Orientation()
# calculate polarizations of all polarizations.
# results are output in pol.txt.
PolarizationTable('pol.txt')
