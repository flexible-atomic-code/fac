from pfac.fac import *
import time

"""
Calculate the collision strengths using the R-matrix method.
"""

#ConvertToSFAC('trmx.sf')

t0 = time.time()

# list of R-matrix basis files, calculated with rmx.py
b = ['rmx.b0', 'rmx.b1']

# list of R-matrix surface files.
r = ['rmx.d0', 'rmx.d1']

# call RMatrixCE, calculate collsion strength at E=750--1.5E3 eV,
# in 100 eV step. Results are in r1.d
RMatrixCE('r1.d', b, r, 563, 583, 2.0)

t1 = time.time()
print 'Time = %12.5E'%(t1-t0)

#CloseSFAC()    
