#
#   FAC - Flexible Atomic Code
#   Copyright (C) 2001-2015 Ming Feng Gu
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <http://www.gnu.org/licenses/>.
#

# this provides an interface to the average atom model

from pfac.fac import *
from pfac import util
import numpy as np
from multiprocessing import Pool
import time, os

class AA:
    """
    for single component plasma, z is scaler, no wm
    for mixtures, z is array of atomic numbers, wm is array of weights by number
    d -- density in g/cc
    t -- temperature in eV
    dd -- directory for the output files
    example:
    Cu at 18 g/cc, 100 eV
    AA(29, 18.0, 100.0).run() 
    C(0.42)H(0.56)Cu(0.02) at 10 g/cc, 200 eV
    AA([1,6,29],10.0,200.0,wm=[0.56,0.42,0.02]).run()
    """
    def __init__(self, z=1, d=1.0, t=10.0, wm=None, dd=None, pref='', nr=8, nc=None):
        if nc is None:
            nc = nr
        if nc > nr:
            nc = nr
        if wm is None:
            self.z = z
            self.d = d
            self.asym = ATOMICSYMBOL[z]
            self.wm = None
        else:
            self.zm = z
            self.dm = d
            wm = np.array(wm)
            self.wm = wm/wm.sum()
            self.ms = np.array([ATOMICMASS[x] for x in z])
            self.mm = np.sum(self.ms*self.wm)
            self.ds = self.dm*(self.ms/self.mm)**(1./3.)
        self.t = t
        self.pref = pref
        self.dd = dd
        self.nr = nr
        self.nc = nc
        if not dd is None:
            if not os.path.exists(dd):
                os.system('mkdir %s'%dd)
        else:
            self.dd = '.'
            
    def aa1p(self, rbf, pref):
        ReinitRadial(0)
        SetAtom(self.asym)
        SetOption('orbital:sc_rbf', rbf)
        AverageAtom(pref, 4, self.d, self.t)

    def ploop(self, i0):
        for i in range(i0, self.nr, self.nc):
            self.aa1p(self.ra[i], self.rs[i])
            
    def aaip(self):
        na = len(self.ra)
        a = self.asym
        b = np.zeros(na)
        c = np.zeros(na)
        for i in range(na):
            r = np.loadtxt('%s/%s%s_aanp_%d.den'%(self.dd,self.pref,a,i), unpack=1)
            d = (r[7]+r[6])/(4*np.pi*r[1]**2)
            d0 = r[8]/(4*np.pi*r[1]**2)
            w = np.where(d > 0)[0]
            i0 = w[-1]
            i1 = w[-2]
            i2 = w[-3]
            dr = r[1][i0]-r[1][i1]
            b[i] = (d[i0]-d[i1])/(0.5*(d[i0]+d[i1])*dr)
            c[i] = (d0[i0]-d0[i1])/(0.5*(d0[i0]+d0[i1])*dr)
        self.rb = b
        self.rc = c
        
    def aanp(self, rmin, rmax):        
        a = self.asym
        nr = self.nr
        nc = self.nc
        self.rmin = rmin
        self.rmax = rmax
        self.ra = np.linspace(rmin, rmax, nr)
        self.rs = ['%s/%s%s_aanp_%d'%(self.dd,self.pref,a,i) for i in range(nc)]
        p = Pool(processes=nc)
        p.map(self.ploop, range(nc))
        p.close()
        p.join()
        self.aaip()

    def run1z(self, rmin=1.0, rmax=2.0, rtol=0.01, rf=1.5):
        niter = 0
        while (True):
            t0 = time.time()
            print('iter beg: %3d %3d %12.5E %12.5E %12.5E %12.5E'%(niter, self.z, self.d, self.t, rmin, rmax))
            self.aanp(rmin, rmax)
            y = self.rb-self.rc
            t1 = time.time()
            if rmax-rmin < rtol*self.nr:
                break
            if y[-1] < 0:
                rmin = rmax
                rmax *= rf
            else:
                w = np.where(y > 0)[0]
                rmin = self.ra[w[0]-1]
                rmax = self.ra[w[0]]
                print('iter end: %3d %3d %12.5E %12.5E %12.5E %12.5E %12.5E'%(niter, w[0], rmin, rmax, y[w[0]-1],y[w[0]], t1-t0))
            niter += 1
        x0 = util.UVIP3P(list(y), list(self.ra), 0.0)
        Print('best rbf: %3d %3d %12.5E %12.5E %12.5E'%(niter, self.z, self.d, self.t, x0))
        self.aa1p(float(x0), '%s/%s%s'%(self.dd,self.pref,self.asym))

    def rden(self, pref, header=None):
        fn = '%s.den'%pref
        if header is None:
            return np.loadtxt(fn, unpack=1)
        rs = np.loadtxt(fn, unpack=1, max_rows=11, comments='@', usecols=1, dtype=str)
        rd = np.loadtxt(fn, unpack=1, max_rows=11, comments='@', usecols=2)
        r = {}
        for i in range(len(rs)):
            r[rs[i][:-1]] = rd[i]
        if len(header) == 0:
            return r
        return r[header]

    def rpot(self, pref, cfg=None):
        fn = '%s.pot'%pref
        if cfg is None:
            return np.loadtxt(fn, unpack=1)
        nc = np.int32(np.loadtxt(fn, unpack=1, skiprows=51, max_rows=1, comments='@', usecols=3))
        d = np.loadtxt(fn, unpack=1, comments='@', skiprows=52, max_rows=nc, usecols=range(1,10))
        if len(cfg) == 0:
            return d
        k = SPECSYMBOL.index(cfg[-1:])
        n = int(cfg[:-1])
        dn = np.int32(d[1])
        dk = np.int32(d[2])

        w1 = np.where((dn == n) & (dk == k))[0]
        w2 = np.where((dn == n) & (dk == -(k+1)))[0]
        
        fb = np.sum(d[3][w1])+np.sum(d[3][w2])
        eb = -(np.sum(d[3][w1]*d[8][w1])+np.sum(d[3][w2]*d[8][w2]))/fb * 27.2
        return fb,eb
            
    def run(self, dtol=0.01, init=True):
        if self.wm is None:
            self.run1z()
            return
        nm = len(self.wm)
        SetOption('orbital:sc_rsf', 1.0)
        SetOption('orbital:sc_rbf', 1.0)
        if init:
            for i in range(nm):
                self.z = self.zm[i]
                self.asym = ATOMICSYMBOL[self.z]
                self.d = self.ds[i]
                self.run1z()
        self.eden = 0.0
        niter = 0
        while (True):
            niter += 1
            eden = 0.0
            for i in range(nm):
                self.z = self.zm[i]
                self.asym = ATOMICSYMBOL[self.z]
                r = self.rden('%s/%s%s'%(self.dd,self.pref,self.asym), header='')
                eden += (self.dm/(1.67*self.mm))*self.wm[i]*r['zb']
            if (abs(eden-self.eden)/eden < dtol):
                break
            eden0 = self.eden
            self.eden = eden
            print('eden beg: %3d %12.5E %12.5E %12.5E %12.5E'%(niter, self.dm, self.t, eden0, eden))
            if niter > 3:
                eden = 0.5*(eden0+eden)
            for i in range(nm):
                self.z = self.zm[i]
                self.asym = ATOMICSYMBOL[self.z]
                self.d = self.ds[i]
                while (True):
                    r = self.rden('%s/%s%s'%(self.dd,self.pref, self.asym), header='')
                    db = r['zb']*r['dn']
                    print('eden itr: %3d %3d %12.5E %12.5E %12.5E %12.5E'%(niter, self.z, eden, db, self.d, self.t))
                    if abs(eden-db)/eden < dtol:
                        break
                    self.d *= eden/db
                    self.run1z()
                    
            
            
