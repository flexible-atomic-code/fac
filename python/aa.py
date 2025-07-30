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
from pfac import util, const, rfac
import numpy as np
from multiprocessing import Pool, cpu_count
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
    def __init__(self, z=1, d=1.0, t=1.0,
                 wm=None, dd=None, pref='',
                 cc=None, znb=0, nr=6, nc=0, sc=0, pmi=0, bqp=-1E12,
                 vxf=2, vxm=2, hxs=-10.0, ngrid=0, maxiter=512,
                 ewm=0, ewf=1.0, vmin=0.2, ztol=-1.0, ids='',
                 ed=0, mmiter=10, mniter=5):
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
            self.mm = ATOMICMASS[z]
            self.ms = [self.mm]
            self.nc = 0
            self.nm = 0
            self.nmr = 0
            self.cc = cc
            if znb == 1 and ed == 1:
                nbe = 0.0
                if not cc is None:
                    nqs = rfac.nlqs(cc)
                    for a in nqs:
                        nbe = nbe + a[2]
                self.d = (d/(z-nbe))*self.mm*1.673
        else:
            self.zm = z
            self.dm = d
            wm = np.array(wm)
            self.wm = wm/wm.sum()
            self.ms = np.array([ATOMICMASS[x] for x in z])
            self.asym = [ATOMICSYMBOL[x] for x in z]
            self.mm = np.sum(self.ms*self.wm)
            self.vt = (1.67e-24*self.mm/self.dm)
            self.ds = np.repeat(self.dm, len(self.wm))
            self.nm = len(self.wm)
            self.nmr = nr*self.nm
            self.ncpu = cpu_count()
            if nc > 0:
                self.nc = min(min(self.ncpu,self.nmr), nc)
            else:
                self.nc = min(self.nmr,self.ncpu)
            if cc is None:
                self.cc = [None for i in range(self.nm)]
            elif type(cc) == type(''):
                self.cc = [cc for i in range(self.nm)]
            else:
                self.cc = cc
            if znb == 1 and ed == 1:
                nfe = 0.0
                for i in range(self.nm):
                    c = self.cc[i]
                    if c is None:
                        nfe = nfe + self.wm[i]*self.zm[i]
                    else:
                        nqs = rfac.nlqs(c)
                        nbe = 0.0
                        for a in nqs:
                            nbe = nbe + a[2]
                        nfe = nfe + self.wm[i]*(self.zm[i]-nbe)
                self.dm = (d/nfe)*self.mm*1.673
        self.t = t
        self.pref = pref
        self.dd = dd
        self.nr = nr
        self.sc = sc
        self.pmi = pmi
        self.bqp = bqp
        self.znb = znb
        self.vxf = vxf
        self.vxm = vxm
        self.hxs = hxs
        self.ewm = ewm
        self.ewf = ewf
        self.maxiter = maxiter
        self.vmin = vmin
        self.ztol = ztol
        self.ngrid = ngrid
        self.ids = ids
        self.miter = 0
        self.niter = 0
        self.mmiter = mmiter
        self.mniter = mniter
        if not dd is None:
            if not os.path.exists(dd):
                os.system('mkdir -p %s'%dd)
        else:
            self.dd = '.'
            
    def aa1p(self, asym, d, t, pref, cc):
        ReinitRadial(0)
        SetAtom(asym)
        if not cc is None:
            if len(cc) > 1:
                AvgConfig(cc)
        if (self.ngrid > 0):
            SetRadialGrid(self.ngrid, -1, -1, -1, -1)
        SetOption('radial:sc_print', self.sc)
        SetOption('radial:print_maxiter', self.pmi)
        SetOptimizeMaxIter(self.maxiter)
        SetOption('radial:sc_vxf', self.vxf)
        SetOption('radial:znbaa', self.znb)
        SetOption('radial:vxm', self.vxm)
        SetOption('orbital:sc_bqp', self.bqp)
        SetOption('orbital:sc_ewm', self.ewm)
        SetOption('orbital:sc_ewf', self.ewf)
        if (not self.hxs is None):
            SetPotentialMode(0, 1e11, 1e11, -1, self.hxs, 0.0)
        AverageAtom(pref, 4, d, t, self.ztol)

    def ploop(self, i0):
        nc = min(self.nmr,max(1,self.nc))
        for i in range(i0, self.nmr, nc):
            print('aa1p: %s %2d %2d %2s %10.3E %10.3E %s %s'%(self.ids,
                                                              self.niter,
                                                              self.miter,
                                                              self.xs[i][0],
                                                              self.xs[i][1],
                                                              self.t,
                                                              self.xs[i][2],
                                                              self.xs[i][3]))
            self.aa1p(self.xs[i][0], self.xs[i][1], self.t,
                      self.xs[i][2], self.xs[i][3])

    def run1z(self, asym):
        pf = '%s/%s%s'%(self.dd,self.pref,asym)
        self.aa1p(asym, self.d, self.t, pf, self.cc)

    def rdos(self, pref, nm=0, emin=-1E31, emax=1E31, emde=0.0, bs=1.0):
        d = np.loadtxt('%s.dos'%pref, unpack=1)
        if nm <= 0 and emin < -1E30 and emax > 1E30:
            return np.array((d[1],d[3]/d[2]))
        
        de = self.rpot(pref, header='ewd')
        if emde > 0:
            emax = emde*de
        c = self.rpot(pref, cfg='')
        if nm > 0:
            w = np.where(c[1] >= nm)[0]
            c = c[:,w]
        if emin > -1E30:
            w = np.where(c[-1] >= emin)[0]
            c = c[:,w]
        if emax < 1E30:
            w = np.where(c[-1] <= emax)[0]
            c = c[:,w]
        if len(c) > 0:    
            xb,yb = rfac.convd(c[-1], c[3], de*bs)
            x = np.arange(xb[0], min(emax,d[1][-1]), min(d[2][0],xb[1]-xb[0]))
            y0 = np.interp(x, xb, yb)
            y1 = np.interp(x, d[1], d[3]/d[2])
            return np.array((x, y0+y1))
        else:
            w = where((d[1] >= emin)&(d[1] <= emax))[0]
            return np.array((d[1][w],d[3][w]/d[2][w]))
        
    def rden(self, pref, header=None):
        fn = '%s.den'%pref
        if header is None:
            return np.loadtxt(fn, unpack=1)
        nw = 0
        with open(fn, 'r') as f:
            rs = f.readlines(20000)
            for i in range(len(rs)):
                if len(rs[i]) < 2 or rs[i][:2] != '# ':
                    nw = i
                    break
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
        return rfac.read_pot(fn, cfg=cfg, header=header)

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
        df = util.UVIP3P(r, list(d[7]), rr);
        db = util.UVIP3P(r, list(d[10]), rr);
        dp = util.UVIP3P(r, list(d[15]), rr)
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
            for i in range(nr):
                f.write('%13.7E '%dp[i]);
                if (i+1)%5 == 0:
                    f.write('\n')

    def rvg(self):
        da = np.zeros((6,self.nm,self.nr))
        for i in range(self.nm):
            for j in range(self.nr):
                pf = '%s/vg%02d_%s%s'%(self.dd,j,self.pref,self.asym[i])
                h = self.rden(pf, header='')
                da[0,i,j] = h['d0']
                da[1,i,j] = h['zf']*h['dn']
                da[2,i,j] = h['ze']*h['dn']
                da[3,i,j] = h['ub']
                da[4,i,j] = h['rps']
                da[5,i,j] = (4*np.pi/3)*(h['rps']*const.RBohr*1e-8)**3
        return da

    def ida(self, k, nd):
        r = self.rvg()
        n0 = np.max(np.min(r[k],1))*0.9
        n1 = np.min(np.max(r[k],1))*1.1
        if n0 >= n1:
            n1,n0 = n0,n1
            ro = 1
        else:
            ro = 0
            
        if n1/n0 < 1.05:
            n1 = n1*1.1
            n0 = n0*0.9
        
        if k==1:
            n0 = np.log10(n0)
            n1 = np.log10(n1)
        xa = np.linspace(n0, n1, nd)        
        ya = np.zeros((self.nm,nd))
        for i in range(self.nm):
            w = np.argsort(r[k,i])
            x0 = r[k,i][w]
            y0 = np.log10(r[-1,i][w])
            if k==1:
                u = np.where(x0 > 0)[0]
                x0 = np.log10(x0[u])
                y0 = y0[u]
            yy = np.array(util.UVIP3P(list(x0), list(y0), list(xa)))
            ya[i] = self.wm[i]*(10**yy)

        vt = sum(ya,0)
        xt = xa
        va = np.log10(vt)
        w = np.argsort(va)
        #xi = np.interp(np.log10(self.vt), va[w], xa[w])
        xi = util.UVIP3P(list(va[w]),
                         list(xa[w]),
                         float(np.log10(self.vt)))
        da = np.zeros(self.nm)
        for i in range(self.nm):
            w = np.argsort(r[k,i])
            x0 = r[k,i][w]
            y0 = np.log10(r[0,i][w])
            if k == 1:
                u = np.where(x0 > 0)[0]
                x0 = np.log10(x0[u])
                y0 = y0[u]
            #yy = np.interp(xi, x0, y0)
            yy = util.UVIP3P(list(x0), list(y0), xi)
            da[i] = 10**yy
        if k==1:
            xa = 10**xa
            xi = 10**xi

        return xi, da, xa, vt, ya, r, ro
        
    def runvg(self, v0, v1):
        self.xs = []
        for i in range(self.nm):
            iva = 10**np.linspace(np.log10(v0[i]), np.log10(v1[i]), self.nr)
            ida = self.ms[i]*1.67e-24/iva
            for j in range(self.nr):
                pf = '%s/vg%02d_%s%s'%(self.dd,j,self.pref,self.asym[i])
                self.xs.append((self.asym[i],ida[j],pf,self.cc[i]))
        self.nmr = self.nm*self.nr
        if (self.nc <= 1):
            self.ploop(0)
        else:
            nc = min(self.nc, self.nmr)
            p = Pool(processes=nc)
            p.map(self.ploop, range(nc))

    def cleanvg(self):
        for i in range(self.nm):
            for j in range(self.nr):
                pf = '%s/vg%02d_%s%s'%(self.dd,j,self.pref,self.asym[i])
                c = 'rm -rf %s.*'%pf
                os.system(c)
                
    def run(self, tol=0.1, imd=1, cvg=1, p=0.0, ptol=0.01):
        if p > 0:
            pref = '%s/%s%s'%(self.dd,self.pref,self.asym)
            self.run()
            pn = self.rden(pref, header='tp')
            biter = 0
            if pn < p:
                p0 = pn
                d0 = self.d
                while (True):
                    biter = biter+1
                    xp = 2*p/pn
                    if xp > 10:
                        xp = 10.0
                    self.d = self.d * xp
                    self.run()
                    pn = self.rden(pref, header='tp')
                    print('iter0: %3d %12.5E %12.5E %12.5E %12.5E %12.5E'%(biter, d0, p0, self.d, pn, xp))
                    if pn >= p:
                        break
                    else:
                        p0 = pn
                        d0 = self.d
                p1 = pn
                d1 = self.d
            elif pn > p:
                p1 = pn
                d1 = self.d
                while (True):
                    biter = biter+1
                    xp = 0.5*p/pn
                    if xp < 0.1:
                        xp = 0.1
                    self.d = self.d * xp
                    self.run()
                    pn = self.rden(pref, header='tp')
                    print('iter1: %3d %12.5E %12.5E %12.5E %12.5E %12.5E'%(biter, d1, p1, self.d, pn, xp))
                    if pn <= p:
                        break
                    else:
                        p1 = pn
                        d1 = self.d
                p0 = pn
                d0 = self.d
            print('biter: %3d %12.5E %12.5E %12.5E %12.5E %12.5E'%(biter, d0, p0, d1, p1, p))
            if abs(1-pn/p) < ptol:
                return
            piter = 0
            xp = np.log(p)
            while ((p1-p0)/p > ptol):
                piter = piter+1
                if p1/p0 > 1e5:
                    dn = 0.5*(d0+d1)
                else:
                    xp0 = np.log(p0)
                    xp1 = np.log(p1)
                    xd0 = np.log(d0)
                    xd1 = np.log(d1)
                    xf = (xp-xp0)/(xp1-xp0)
                    dn = np.exp(xd0*(1-xf) + xd1*xf)
                self.d = dn
                self.run()
                pn = self.rden(pref, header='tp')
                print('piter: %3d %12.5E %12.5E %12.5E %12.5E %12.5E %12.5E'%(piter, d0, d1, dn, pn, p, 1-pn/p))
                if abs(1-pn/p) < ptol:
                    return
                if pn < p:
                    d0 = dn
                    p0 = pn
                else:
                    d1 = dn
                    p1 = pn
                if piter > 100:
                    print('pressure bisection failed')
                    break
            return
        if imd <= 0:
            self.irun(tol=tol)
            return
        t0 = time.time()    
        if self.wm is None:
            self.run1z(self.asym)
            return
        v0 = np.zeros(self.nm)
        v1 = v0.copy()
        for i in range(self.nm):
            v0[i] = self.vmin*self.vt
            v1[i] = (1.25/self.wm[i])*self.vt
        niter = 0
        while (True):
            niter += 1
            miter = 0
            self.niter = niter
            for i in range(self.nm):
                print('aavm: %s %2d %2d %2s %12.5E %12.5E'%(self.ids, niter, miter, self.asym[i], v0[i], v1[i]))
            self.runvg(v0, v1)
            x,d,xa,va,ya,ra,ro = self.ida(imd, self.nr*2)
            va0 = np.min(va)
            va1 = np.max(va)
            while (True):
                miter += 1
                self.miter = miter
                if miter >= self.mmiter:                    
                    print('aavf: %s %2d %2d %10.3E %10.3E %10.3E %10.3E'%(self.ids,niter,miter,self.vt,va0,va1,time.time()-t0))
                    return
                print('aavl: %s %2d %2d %10.3E %10.3E %10.3E %10.3E'%(self.ids, niter,miter,self.vt, va0, va1, time.time()-t0))
                if ro > 0 or self.vt < va0 or self.vt > va1:
                    x0 = np.log(va)
                    y0 = np.log(self.vt)
                    w = np.argsort(x0)
                    for i in range(self.nm):
                        x1 = np.log(ya[i])
                        y1 = np.exp(util.UVIP3P(list(x0[w]), list(x1[w]), y0))
                        y1 = y1/self.wm[i]
                        v1[i] = min(1.25*self.vt/self.wm[i],
                                    y1/np.sqrt(self.vmin))
                        v0[i] = v1[i]*self.vmin
                else:
                    break        
                for i in range(self.nm):
                    print('aavm: %s %2d %2d %2s %12.5E %12.5E'%(self.ids, niter, miter, self.asym[i], v0[i], v1[i]))
                self.runvg(v0, v1)
                x,d,xa,va,ya,ra,ro = self.ida(imd, self.nr*2)
                va0 = np.min(va)
                va1 = np.max(va)
            self.xs = []
            for i in range(self.nm):
                pf = '%s/%s%s'%(self.dd,self.pref,self.asym[i])
                self.xs.append((self.asym[i],d[i],pf,self.cc[i]))
            self.nmr = self.nm
            if (self.nc <= 1):
                self.ploop(0)
            else:
                nc = min(self.nc,self.nm)
                p = Pool(processes=nc)
                p.map(self.ploop, range(nc))
            ys = np.zeros(self.nm)
            vi = ys.copy()
            for i in range(self.nm):
                pf = '%s/%s%s'%(self.dd,self.pref,self.asym[i])
                h = self.rden(pf, header='')
                if imd == 1:
                    ys[i] = h['zf']*h['dn']
                elif imd == 2:
                    ys[i] = h['ze']*h['dn']
                else:
                    ys[i] = h['ub']
                vi[i] = (4*np.pi/3)*(h['rps']*const.RBohr*1e-8)**3
            ym = np.mean(ys)
            dy = (np.max(ys)-np.min(ys))/2.0
            if (dy <= tol*abs(ym)):
                print('aacv: %s %2d %2d %10.3E %10.3E %10.3E %10.3E %s'%(self.ids,niter,miter,dy,ym,x,time.time()-t0,self.dd))
                break
            w = np.argsort(va)
            x0 = np.log10(va)
            for i in range(self.nm):
                y0 = np.log10(ya[i])
                x1 = np.log10(np.array([self.vt/2.0, self.vt*2.0]))
                yy = np.array(util.UVIP3P(list(x0[w]), list(y0[w]), list(x1)))
                v0[i] = 10**min(yy)/self.wm[i]
                v1[i] = 10**max(yy)/self.wm[i]
            if niter > self.mniter:
                print('aacf: %s %2d %2d %10.3E %10.3E %10.3E'%(self.ids,niter,miter,dy,ym,x))
                break
            print('aait: %s %2d %2d %10.3E %10.3E %10.3E %10.3E'%(self.ids,niter,miter,dy,ym,x,time.time()-t0))
        if cvg:
            self.cleanvg()
        
    def irun(self, dtol=0.1, init=True):
        if self.wm is None:
            self.run1z(self.asym)
            return
        t0 = time.time()
        nm = len(self.wm)
        dtol1 = dtol*1.25
        dtol2 = dtol/1.25
        SetOption('orbital:sc_rsf', 1.0)
        SetOption('orbital:sc_rbf', 1.0)
        if init:
            for i in range(nm):
                self.z = self.zm[i]
                self.d = self.ds[i]
                print('init run: %s %3d %12.5E %12.5E %12.5E %10.3E'%(self.ids,self.z, self.dm, self.d, self.t, time.time()-t0))
                self.run1z(self.asym[i])
        self.eden = 0.0
        niter = 0
        wst0 = 0.9
        wst1 = 0.1
        while (True):
            x = np.exp(-niter/5.)
            wst = wst0*x + wst1*(1-x)
            niter += 1
            if niter > self.maxiter:
                print('maxiter reached in outer loop: %d %s'%(niter,self.ids))
                return
            eden = 0.0
            vt = 0.0
            for i in range(nm):
                self.z = self.zm[i]
                asym = self.asym[i]
                r = self.rden('%s/%s%s'%(self.dd,self.pref,asym), header='')
                eden += (self.dm/(1.67*self.mm))*self.wm[i]*max(0.0,r['zf'])
                vt += self.wm[i]*(r['rps']*const.RBohr*1e-8)**3
            vt *= 4*np.pi/3.0
            vr = vt/self.vt
            done = 0
            if (abs(eden-self.eden) < dtol*max(1E-24,eden) and
                abs(vr-1.0) < dtol1):
                done = 1
            eden0 = self.eden
            self.eden = eden
            print('eden beg: %3d %s %12.5E %12.5E %12.5E %12.5E %10.3E %10.3E %10.3E'%(niter, self.ids, self.dm, self.t, eden0, eden, vr, wst, time.time()-t0))
            if done:
                break
            if niter > 1:
                eden = eden0*(1-wst) + eden*wst
            b = vr**wst
            self.ds *= b
            eden *= b
            for i in range(nm):
                self.z = self.zm[i]
                asym = self.asym[i]
                self.d = self.ds[i]
                ni = 0
                while (True):
                    x = np.exp(-ni/5.)
                    wst = wst0*x + wst1*(1-x)
                    ni += 1
                    if ni > self.maxiter:
                        print('maxiter reached in inner loop: %s %d %d'%(self.ids,i,ni))
                        return
                    r = self.rden('%s/%s%s'%(self.dd,self.pref, asym),
                                  header='')
                    db = abs(r['zf'])*r['dn']
                    z0 = eden/r['dn']
                    zb = db/r['dn']
                    print('eden itr: %3d %3d %s %3d %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E %10.3E'%(niter, ni, self.ids, self.z, eden, db, z0, zb, self.d, self.t, wst, time.time()-t0))
                    if abs(z0-zb)/max(1e-3,z0) < dtol2:
                        break
                    db = max(1e-3*eden,db)
                    if db > 0:
                        dn = self.d*(eden/db)
                        self.d = self.d*(1-wst) + dn*wst
                        self.run1z(asym)                    
                self.ds[i] = self.d

# AA model at fixed pressure or mass density. t in eV, p in Mbar, d in g/cc
def AAPD(z, t, d=1.0, p=0.0, dd=None, pref='', ptol=0.01, bqp=1e30):
    a = AA(z=z, t=t, d=d, dd=dd, pref=pref, bqp=bqp)
    a.run(p=p, ptol=ptol)
    
