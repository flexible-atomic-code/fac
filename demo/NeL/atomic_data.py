from pfac.atom import *
from pfac import fac
import os

asym = 'Ne'

# generate atomic data for He-like ions.
nele = list(range(2, 3))

# nterms: number 
atomic_data(nele, asym, iprint=1, dir='./', nterms=[0, 4],
            nexc_max=[3, 3, 2], nrec_max=[3, 3, 2])
