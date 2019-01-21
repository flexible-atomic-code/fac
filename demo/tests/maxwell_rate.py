from pfac.fac import *

# atomic number, number of electrons, number of excitation levels
z = 10
nele = 3
nmax = 5

a = 'Ne'  # atomic symbol (Ne)

SetAtom(a)

Config('g2', '1s2 2*1')
Config('g3', '1s2 3*1')
Config('ion', '1s2')

# radial potential
ConfigEnergy(0)
OptimizeRadial('g2')
ConfigEnergy(1)

# atomic structure
# also write Hamiltonian to file
Structure(a+'b.en', a+'b.ham', ['g2', 'g3'])
Structure(a+'b.en', a+'b.ham', ['ion'])

MemENTable(a+'b.en')
PrintTable(a+'b.en', a+'a.en')

# transition rates
TRTable(a+'b.tr', ['g2'], ['g3'])

# excitation
CETable(a+'b.ce', ['g2'], ['g3'])

# ionization
CITable(a+'b.ci', ['g2'], ['ion'])
CITable(a+'b.ci', ['g3'], ['ion'])

# write maxwell rate
MaxwellRate('Neb.ci', 'Neb.cir', -1, -1, 1000)
MaxwellRate('Neb.ce', 'Neb.cer', -1, -1, 1000)
