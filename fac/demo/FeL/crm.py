from pfac.crm import *

# by uncommenting the following line and the one at the end of this file,
# this script will be converted to the input file for the SCRM interface.
# ConvertToSCRM('crm.sf')

# This script construct a three-ion model for Fe L-shell ions.
# it must be run after fe.py, which generates the necessary atomic data.

r = [1.0]#, 1.5, 2.0, 2.5, 3.0]
den = 1e1
z = 26
k = 10
(t1, a1) = MaxAbund(z, k-1)
(t2, a2) = MaxAbund(z, k)
(t3, a3) = MaxAbund(z, k+1)
temp = [t1, t2, t3]
population = [a1, a2, a3]

f1 = 'Fe%02db'%(k-1)
f2 = 'Fe%02db'%k
f3 = 'Fe%02db'%(k+1)
AddIon(k, 0.0, f2) 
AddIon(k+1, 0.0, f3) 
SetBlocks(0.0, f1) 

for i in range(len(temp)):
    print 'Temp = %10.3E'%(temp[i])
    p1 = population[i][k-1]
    p2 = population[i][k]
    p3 = population[i][k+1]
    
    SetEleDist(0, temp[i])
    print 'CE rates...'
    SetCERates(1)
    print 'TR rates...'
    SetTRRates(0)
    print 'RR rates...'
    SetRRRates(0)
    print 'CI rates...'
    SetCIRates(0)
    print 'AI rates...'
    SetAIRates(1)
    
    SetEleDensity(den)
    print 'Init blocks...'
    InitBlocks()
    for j in range(len(r)):
        print 'Abund: %10.3E %10.3E %10.3E'%(p1*r[j], p2, p3)
        s = 't%da%d'%(i,j)
        rt_file = '%s_%s.rt'%(f2,s)
        sp_file = '%s_%s.sp'%(f2,s)
        rt_afile = '%sa_%s.rt'%(f2[:-1],s)
        sp_afile = '%sa_%s.sp'%(f2[:-1],s)
        plot_file = '%ssp_%s.d'%(f2[:-1],s)
    
        SetAbund(k-1, p1*r[j])
        SetAbund(k, p2)
        SetAbund(k+1, p3)

        LevelPopulation()
        Cascade()
        RateTable(rt_file)
        SpecTable(sp_file)
    
        PrintTable(rt_file, rt_afile, 1)
        PrintTable(sp_file, sp_afile, 1)

        PlotSpec(sp_file, plot_file, 1, 500.0, 2.0e3, 1.0)
        ReinitCRM(2)

    ReinitCRM(1)
    
# CloseSCRM()

