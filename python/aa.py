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

# decode zs and ws arrays for compound.
# SiO2 -> zs=[14,8],ws=[1,2]
# H0.25C0.7Cu0.05 -> zs=[1,6,29], ws=[0.25,0.7,0.05]
def zw4c(s):
    a = ATOMICSYMBOL
    s = s+'X'
    n = len(s)
    zs = []
    ws = []
    i = 0
    j = 0
    for k in range(i+1,n):
        if j == i and s[k].isdigit():
            j = k
            continue
        if s[k].isupper():
            if j == i:
                zs.append(a.index(s[i:k]))
                ws.append(1.0)
            else:
                zs.append(a.index(s[i:j]))
                ws.append(float(s[j:k]))
            i = k
            j = k
    return zs,ws

class AA:
    """
    for single component plasma, z is integer for atomic number, no wm
    for mixtures, z is array of atomic numbers, 
                  wm is array of weights by number
                  or with no wm, but z is a chemical formula, 
                  e.g., CO2, H2O, H0.56C0.42Cu0.02
    d -- density in g/cc
    t -- temperature in eV
    dd -- directory for the output files
    example:
    Cu at 18 g/cc, 100 eV
    AA(29, 18.0, 100.0).run() 
    C(0.42)H(0.56)Cu(0.02) at 10 g/cc, 200 eV
    AA([1,6,29],10.0,200.0,wm=[0.56,0.42,0.02]).run()
    or
    AA('H0.56C0.42Cu0.02', 10.0, 200).run()
    """
    def __init__(self, z=1, d=1.0, t=1.0, wm=None, dd=None, pref='',
                 nr=8, nc=None, sc=0, bqp=-1E12, vxf=2, hxs=0.67):
        if nc is None:
            nc = nr
        if nc > nr:
            nc = nr
        if type(z) == type(''):
            z,wm = zw4c(z)
            if len(z) == 1:
                z = z[0]
                wm = None
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
            self.ds = np.repeat(self.dm, len(self.wm))
        self.t = t
        self.pref = pref
        self.dd = dd
        self.nr = nr
        self.nc = nc
        self.sc = sc
        self.bqp = bqp
        self.vxf = vxf
        self.hxs = hxs
        if not dd is None:
            if not os.path.exists(dd):
                os.system('mkdir %s'%dd)
        else:
            self.dd = '.'
            
    def aa1p(self, rbf, pref):
        ReinitRadial(0)
        SetAtom(self.asym)
        if self.sc > 0:
            SetOption('radial:sc_print', 1)
        SetOption('radial:sc_vxf', self.vxf)
        SetOption('orbital:sc_bqp', self.bqp)
        SetOption('orbital:sc_rbf', rbf)
        if (not self.hxs is None):
            SetPotentialMode(0, 1e11, 1e11, -1, self.hxs, 0.0)
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
            fp = '%s/%s%s_aanp_%d'%(self.dd,self.pref,a,i)
            r = self.rden(fp)
            d = (r[4]+r[6]+r[7])/(4*np.pi*r[1]**2)
            d0 = r[8]/(4*np.pi*r[1]**2)
            rs = self.rpot(fp, header='rps')
            w = np.where(r[1] <= rs)[0]
            i0 = w[-1]
            i1 = w[-2]
            dr = (r[1][i0]-r[1][i1])/rs
            if (d[i0] > 0 and d[i1] > 0):
                b[i] = (d[i0]-d[i1])/(0.5*(d[i0]+d[i1])*dr)
            else:
                b[i] = 0.0
            if (d0[i0] > 0 and d0[i1] > 0):
                c[i] = (d0[i0]-d0[i1])/(0.5*(d0[i0]+d0[i1])*dr)
            else:
                c[i] = 0.0
        self.rb = b
        self.rc = c
        
    def aanp(self, rmin, rmax):        
        a = self.asym
        nr = self.nr
        nc = self.nc
        self.rmin = rmin
        self.rmax = rmax
        self.ra = np.linspace(rmin, rmax, nr)
        self.rs = ['%s/%s%s_aanp_%d'%(self.dd,self.pref,a,i) for i in range(nr)]
        p = Pool(processes=nc)
        p.map(self.ploop, range(nc))
        p.close()
        p.join()
        self.aaip()

    def run1z(self, rmin=1.0, rmax=2.0, rtol=0.01, rf=2.0, nmr=0):
        pf = '%s/%s%s'%(self.dd,self.pref,self.asym)
        if nmr == 0:
            self.aa1p(1.0, pf)
            return
        niter = 0
        nm = 1
        while (nm <= nmr):
            t0 = time.time()
            print('iter beg: %3d %3d %12.5E %12.5E %12.5E %12.5E'%(niter, self.z, self.d, self.t, rmin, rmax))
            self.aanp(rmin, rmax)
            y = self.rb-self.rc
            t1 = time.time()
            if rmax-rmin < rtol*self.nr:
                break
            if y[-1] < 0:
                nm += 1
                w = len(y)-1
                if (nm < nmr):
                    rmin = rmax
                    rmax *= rf
            else:
                w = np.where(y >= 0)[0][0]
                rmin = self.ra[w-1]
                rmax = self.ra[w]
            print('iter end: %3d %3d %12.5E %12.5E %12.5E %12.5E %12.5E'%(niter, w, rmin, rmax, y[w-1],y[w], t1-t0))
            niter += 1
        w = np.where(y >= 0)[0]
        if len(w) > 0:
            w = w[0]
            k = w-1
            xf = y[w]/(y[w]-y[k])
            x0 = self.ra[k]*xf + self.ra[w]*(1-xf)
            nx = -1
        else:
            nx = self.nr-1
            x0 = self.ra[nx]
        Print('best rbf: %3d %3d %12.5E %12.5E %12.5E'%(niter, self.z, self.d, self.t, x0))
        if nx < 0:
            self.aa1p(float(x0), pf)
        else:
            os.system('cp %s_aanp_%d.pot %s.pot'%(pf,nx,pf))
            os.system('cp %s_aanp_%d.den %s.den'%(pf,nx,pf))

    def rden(self, pref, header=None):
        fn = '%s.den'%pref
        if header is None:
            return np.loadtxt(fn, unpack=1)
        rs = np.loadtxt(fn, unpack=1, max_rows=100, comments='@', usecols=0, dtype=str)
        w = np.where(rs == '#')[0]
        nw = len(w)
        rs = np.loadtxt(fn, unpack=1, max_rows=nw, comments='@', usecols=1, dtype=str)
        rd = np.loadtxt(fn, unpack=1, max_rows=nw, comments='@', usecols=2)
        r = {}
        for i in range(len(rs)):
            r[rs[i][:-1]] = rd[i]
        if len(header) == 0:
            return r
        return r[header]

    def rpot(self, pref, cfg=None, header=None):
        fn = '%s.pot'%pref
        if cfg is None and header is None:
            return np.loadtxt(fn, unpack=1)
        rs = np.loadtxt(fn, unpack=1, max_rows=100, comments='@', usecols=1, dtype=str)
        w = np.where(rs == 'Mean')[0]
        nw = w[0]
        if not header is None:
            rs = np.loadtxt(fn, unpack=1, max_rows=nw, comments='@', usecols=1, dtype=str)
            rd = np.loadtxt(fn, unpack=1, max_rows=nw, comments='@', usecols=3)
            r = {}
            for i in range(len(rs)):
                r[rs[i]] = rd[i]
            if len(header) == 0:
                return r
            return r[header]
        nc = np.int32(np.loadtxt(fn, unpack=1, skiprows=nw, max_rows=1, comments='@', usecols=3))
        d = np.loadtxt(fn, unpack=1, comments='@', skiprows=nw+1, max_rows=nc, usecols=range(1,10))
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

    def wden(self, pref, nr, ofn, rmin=None):
        hd = self.rden(pref, header='')
        d = self.rden(pref)
        z = int(hd['zps'])
        te = hd['T']
        de = hd['d0']
        r = d[1]
        rs = hd['rps']
        if rmin is None:
            rmin = r[0]
        h = np.log(rs/rmin)/(nr-1)
        xr = np.arange(nr)
        rr = list(rmin*np.exp(xr*h))
        r = list(r)        
        da = util.UVIP3P(r, list(d[4]+d[9]), rr);
        df = util.UVIP3P(r, list(d[4]+d[7]), rr);
        db = util.UVIP3P(r, list(d[10]), rr);
        with open(ofn, 'w') as f:
            f.write('%5d %6d %6d\n'%(z, 1, 1))
            f.write('%12.5E %12.5E\n'%(te, de))
            f.write('%4d %12.5E %12.5E\n'%(nr, h, rs))
            for i in range(nr):
                f.write('%13.7E '%rr[i])
                if (i+1)%5 == 0:
                    f.write('\n')
            for i in range(nr):
                f.write('%13.7E '%1.0)
                if (i+1)%5 == 0:
                    f.write('\n')
            for i in range(nr):
                f.write('%13.7E '%da[i])
                if (i+1)%5 == 0:
                    f.write('\n')
            for i in range(nr):
                f.write('%13.7E '%df[i])
                if (i+1)%5 == 0:
                    f.write('\n')
            for i in range(nr):
                f.write('%13.7E '%db[i]);
                if (i+1)%5 == 0:
                    f.write('\n')
                    
    def run(self, dtol=0.05, init=True):
        if self.wm is None:
            self.run1z()
            return
        t0 = time.time()
        nm = len(self.wm)
        SetOption('orbital:sc_rsf', 1.0)
        SetOption('orbital:sc_rbf', 1.0)
        if init:
            for i in range(nm):
                self.z = self.zm[i]
                self.asym = ATOMICSYMBOL[self.z]
                self.d = self.ds[i]
                print('init run: %3d %12.5E %12.5E %12.5E %10.3E'%(self.z, self.dm, self.d, self.t, time.time()-t0))
                self.run1z()
        self.eden = 0.0
        niter = 0
        wst0 = 0.9
        wst1 = 0.1
        while (True):
            x = np.exp(-niter/5.)
            wst = wst0*x + wst1*(1-x)
            niter += 1
            eden = 0.0
            for i in range(nm):
                self.z = self.zm[i]
                self.asym = ATOMICSYMBOL[self.z]
                r = self.rden('%s/%s%s'%(self.dd,self.pref,self.asym), header='')
                eden += (self.dm/(1.67*self.mm))*self.wm[i]*max(0.0,r['zf'])
            done = 0
            if (abs(eden-self.eden) < dtol*max(1E-24,eden)):
                done = 1
            eden0 = self.eden
            self.eden = eden
            print('eden beg: %3d %12.5E %12.5E %12.5E %12.5E %10.3E %10.3E'%(niter, self.dm, self.t, eden0, eden, wst, time.time()-t0))
            if done:
                break
            if niter > 1:
                eden = eden0*(1-wst) + eden*wst
            for i in range(nm):
                self.z = self.zm[i]
                self.asym = ATOMICSYMBOL[self.z]
                self.d = self.ds[i]
                ni = 0
                while (True):
                    x = np.exp(-ni/5.)
                    wst = wst0*x + wst1*(1-x)
                    ni += 1                    
                    r = self.rden('%s/%s%s'%(self.dd,self.pref, self.asym),
                                  header='')
                    db = abs(r['zf'])*r['dn']
                    z0 = eden/r['dn']
                    zb = db/r['dn']
                    print('eden itr: %3d %3d %3d %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E'%(niter, ni, self.z, eden, db, z0, zb, self.d, self.t, wst, time.time()-t0))
                    if abs(z0-zb)/max(1e-3,z0) < dtol/2:
                        break
                    db = max(1e-3*eden,db)
                    if db > 0:
                        dn = self.d*(eden/db)
                        self.d = self.d*(1-wst) + dn*wst
                        self.run1z()                    
                self.ds[i] = self.d
            
