from pfac import crm
from pfac.crm import *

# by uncommenting the following line and the one at the end of this file,
# this script will be converted to the input file for the SCRM interface.
# ConvertToSCRM('crm.sf')

# This script must be run after fe.py, which generates the
# necessary atomic data.

# include Fe XVI -- XXII, with the relative abundance
# specified by population.
neles = [5, 6, 7, 8, 9, 10, 11]
population = [0.05, 0.06, 0.11, 0.13, 0.25, 0.34, 0.06]
# the population of the ion NELE=4, which is needed 
# to construct the spectrum for NELE=5. leave it 0.0 to 
# let CRM determine its abundance.
p2 = 0.0

# temperature = 600 eV, electron density = 10^12 cm^{-3}.
# the code uses the unit of 10^10 cm^{-3}
temp = 600.0
den = [1e2]

for k in range(len(neles)):
    nele = neles[k]
    p1 = population[k]
    s = 'NELE = %d, Population = %10.3E'%(nele, p1)
    Print(s)
    if (k > 0):
        p2 = population[k-1]
    
    f1 = 'Fe%02db'%nele
    f2 = 'Fe%02db'%(nele-1)
    rt_file = '%s.rt'%f1
    sp_file = '%s.sp'%f1
    rt_afile = '%sa.rt'%f1[:-1]
    sp_afile = '%sa.sp'%f1[:-1]
    plot_file = '%ssp.d'%f1[:-1]

    AddIon(nele, p1, f1)
    if (nele == 1):
        SetBlocks(p2)
    else:
        SetBlocks(p2, f2)

    SetEleDist(0, temp)
    SetCERates(1)
    SetTRRates(0)
    SetRRRates(0)
    SetCIRates(0)
    SetAIRates(1)

    for i in range(len(den)):
        s = 'Electron Density = %10.3E\n'%(den[i])
        Print(s)
        SetEleDensity(den[i])
        InitBlocks()
        LevelPopulation()
        RateTable(rt_file)
        SpecTable(sp_file)

    PrintTable(rt_file, rt_afile, 1)
    PrintTable(sp_file, sp_afile, 1)

    PlotSpec(sp_file, plot_file, 1, 500.0, 2.0e3, 1.0)

    ReinitCRM()
    
# CloseSCRM()
