from pfac.pol import *

# setup the magnetic sublevel table and transition rates.
SetMLevels('tb.en', 'tb.tr')
# setup the magnetic sublevel excitation rates at E=12.24 keV (900 Ryd)
SetMCERates('tb.ce', 1.224E4)
# calculates the magnetic sublevel populations.
# results are output in pop.txt.
PopulationTable('pop.txt', 1.0)
# calculate polarizations of all polarizations.
# results are output in pol.txt.
PolarizationTable('pol.txt')
