from pfac.atom import *
from pfac import fac
import os

# by uncommenting the following line and the one at the end of this file,
# this script will be converted to the input file for the SFAC interface.
#fac.ConvertToSFAC('fe.sf')

asym = 'Fe'

# generate atomic data for H-like to Na-like ions.
os.system('mkdir data')
nele = range(1, 2)
atomic_data(nele, asym, iprint=1, dir='data/')

#fac.CloseSFAC()
