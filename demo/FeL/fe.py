from pfac.atom import *
from pfac import fac

# by uncommenting the following line and the one at the end of this file,
# this script will be converted to the input file for the SFAC interface.
#fac.ConvertToSFAC('fe.sf')

asym = 'Fe'

# generate atomic data for H-like to Na-like ions.
nele = range(1,12)
atomic_data(nele, asym, iprint=1)

#fac.CloseSFAC()
