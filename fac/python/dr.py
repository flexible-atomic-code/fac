import string
import array
import const
from math import sqrt, exp, pi, log10, log
import fac

def read_lev(fn):
    energy = array.array('f')
    f = open(fn, "r")
    while (1):
        s = f.readline()
        if not s:
            break
        s = string.split(s)
        try:
            i = long(s[0])
            energy.append(float(s[1]))
        except:
            continue
    f.close()
    return energy

def read_ai(fn, free):
    channels = []
    ai = array.array('f')
    tai = array.array('f')
    e = array.array('f')
    bound = array.array('i')
    f = open(fn, 'r')
    valid_bound = 0
    while (1):
        s = f.readline()
        if not s:
            break
        if string.find(s, 'DR Channel:') >= 0:
            s = string.split(s)
            channels.append([long(s[2]), 0])
            continue
        if string.find(s, 'Total') >= 0:
            if not valid_bound:
                continue
            s = string.split(s)
            tai.append(float(s[2]))
            channels[-1][1] = channels[-1][1] + 1
            valid_bound = 0
            continue
        s = string.split(s)
        try:
            if long(s[2]) == free:
                bound.append(long(s[0]))
                e.append(float(s[4]))
                ai.append(float(s[6]))
                valid_bound = 1
        except:
            continue
        
    f.close()

    return channels, bound, e, ai, tai

def read_tr(fn, up):
    rd = array.array('f')
    trd = array.array('f')
    trdp = array.array('f')
    te = array.array('f')
    nlow = array.array('i')
    low = array.array('i')
    for i in range(len(up)):
        nlow.append(0)

    f = open(fn, 'r')
    valid_bound = 0
    i = 0
    while (1):
        s = f.readline()
        if not s:
            break
        if string.find(s, 'Total') >= 0:
            if not valid_bound:
                continue
            s = string.split(s)
            trd.append(float(s[2]))
            trdp.append(float(s[3]))
            i = i+1
            valid_bound = 0
            continue
        s = string.split(s)
        try:
            k = long(s[0])
            while up[i] < k:
                i = i+1
                trd.append(0.0)
                trdp.append(0.0)
                
            if k == up[i]:
                valid_bound = 1
                nlow[i] = nlow[i]+1
                low.append(long(s[2]))
                te.append(float(s[4]))
                rd.append(float(s[7]))
        except:
            continue

    f.close()
    return nlow, low, te, rd, trd, trdp

def branching(tai, trdp, trd, *cas):
    br = array.array('f')
    nb = len(tai)

    i = 0

    while i < nb:
        br.append(trdp[i]/(tai[i]+trdp[i]))
        i = i+1

    n = len(cas)
    if (n == 0):
        return br
    
    if (n != 4):
        print "Including cascades needs 4 additonal args"
        print "skipped"
        return br
        
    #bound = cas[0]
    #nlow = cas[1]
    #low = cas[2]
    #rd = cas[3]

    j0 = 0
    i = 0
    while i < nb:
        if (trd[i] < 1E-30):
            i = i+1
            continue
        if (trd[i]-trdp[i]) < (1E-2)*trd[i]:
            i = i+1
            continue
        t = 0.0
        j = j0
        j1 = cas[1][i]+j0
        while j < j1:
            if cas[2][j] < 0:
                j = j+1
                continue
            try:
                k = (cas[0]).index(cas[2][j])
            except:
                j = j+1
                continue
            t = t+cas[3][j]*br[k]
            j = j+1
        j0 = j1
        if (trdp[i] < 1E-10*trd[i]):
            br[i] = t/(tai[i]+trd[i])
        else:
            br[i] = br[i]*((trdp[i] + t)/(trdp[i] + br[i]*(trd[i]-trdp[i])))
            
        i = i+1
        
    return br

def drstrength(e, ai, br):
    s = 0.0
    for i in range(len(ai)):
        s = s + e[i]*ai[i]*br[i]
    return s

def drcross(outfile, e0, ai, br, tai, trd, de, vm = 0):
    emin = min(e0)
    emax = max(e0)
    dei = emax-emin

    sig8 = 10.0*de
    emin = emin - sig8
    emax = emax + sig8
    
    dei = 0.1*de

    ene = array.array('f')
    cross = array.array('f')

    f = open(outfile, 'w')
    f.write('# Energy (eV) DR_Cross (10^-20 cm2)\n')

    if vm > 0:
        vf = const.c10 * sqrt(2.0/const.Me_eV)
        
    e = emin
    while e < emax:
        c = 0.0;
        for i in range(len(e0)):
            y = e - e0[i]
            if (abs(y) > sig8):
                continue
            w = (tai[i] + trd[i])*const.Ryd_eV*2.0
            y = voigt(y, de, w)
            y = y*ai[i]*br[i]
            if (vm > 0):
                y = y*sqrt(e0[i])*vf
            c = c + y
            
        if (c > 1E-10):
            ene.append(e)
            cross.append(c)
            f.write('%11.4E %11.4E\n'%(e, c))
        e = e + dei
 
    f.close()
    return ene, cross

def voigt(x, wg, wl):
    if wg > 10*wl:
        delta = 1.0
        wv = wg
    elif wl > 10*wg:
        delta = 0.0
        wv = wl
    else:
        wv = wl + sqrt(wl*wl+4.0*wg*wg)
        wv = wv*0.5
        delta = (wg/wv)**2

    ln2 = log(2.0)
    norm = (1.0-delta)*pi + delta*sqrt(pi/ln2)
    norm = norm*wv*0.5

    t = 2.0*x/wv
    t = t*t
    v = delta*exp(-t*ln2) + (1.0-delta)/(1.0+t)
    v = v/norm
    return v
    
def drrate(outfile, e, ai, br, tmin, tmax, ntemp):
    tmin = log10(tmin)
    tmax = log10(tmax)
    deltat = (tmax-tmin)/(ntemp-1.0)
    temp = array.array('f')
    temp.append(tmin)
    for i in range(ntemp-1):
        temp.append(temp[i]+deltat)
    for i in range(ntemp):
        temp[i] = 10**(temp[i])
        
    rate = array.array('f')
    for i in range(ntemp):
        rate.append(0.0)

    p = 4.787307/sqrt(const.Me_eV)
    for i in range(ntemp):
        t = temp[i]
        y = p / pow(t, 1.5)
        for j in range(len(e)):
            x = e[j]/t
            if x > 60.0:
                continue
            x = exp(-x)
            x = x*e[j]*y
            x = x*ai[j]*br[j]
            rate[i] = rate[i]+x

    f = open(outfile, 'w')
    f.write('#Temp (eV) DR_RateCoeff(10^-10 cm3/s)\n')
    for i in range(ntemp):
        f.write('%10.3E %10.3E\n'%(temp[i], rate[i]))
    f.close()
    return temp, rate

def sumrate(outfile, prefix, nrec, nopen, dn, nmax = 0, f = 0):
    rates = []
    if (len(nrec) == 0):
        return
    
    for n in nrec:
        rt_file = constructfn(prefix, 'rt', n, dn, f)
        (t, r) = readcol(rt_file)
        rates.append(r)
    nt = len(t)
    range_nrec = range(len(nrec))
    
    if nmax > 0 and nmax < nrec[-1]:
        nlimit = nmax
    else:
        nlimit = nrec[-1]
    
    trate = array.array('f',[0.0])
    trate = trate*nt

    n = range(nrec[0], nlimit+1)
    nmissing = []
    for i in n:
        if not i in nrec:
            nmissing.append(i)
    log_nrec = []
    for i in nrec:
        log_nrec.append(log10(i))
    for s in range(nt):
        log_rt = []
        for i in range_nrec:
            if (nmax == 0 or (nmax > 0 and i <= nmax)):
                trate[s] = trate[s] + rates[i][s]
            if rates[i][s] > 0.0:
                log_rt.append(log10(rates[i][s]))
            else:
                log_rt.append(-50.0)

        k0 = 0
        for k in nopen:
            try:
                k1 = nrec.index(k)
            except:
                continue
            if k1 == k0:
                continue
            y2 = fac.Spline(log_nrec[k0:k1], log_rt[k0:k1])
            for nm in nmissing:
                if nrec[k0] < nm < nrec[k1]:
                    r = fac.Splint(log_nrec[k0:k1], log_rt[k0:k1],
                                   y2, log10(nm))
                    r = pow(10, r)
                    trate[s] = trate[s]+r
                    
            k0 = k1

        if k0 < len(nrec)-1:
            y2 = fac.Spline(log_nrec[k0:], log_rt[k0:])
            for nm in nmissing:
                if nm > nrec[k0]:
                    r = fac.Splint(log_nrec[k0:], log_rt[k0:],
                                   y2, log10(nm))
                    r = pow(10, r)
                    trate[s] = trate[s]+r
                
        if (nmax > nrec[-1]):
            f = extrapolate(nrec[-1], nmax)
            r = rates[-1][s]*f
            trate[s] = trate[s]+r 
                
    f = open(outfile, 'w')
    for i in range(nt):
        f.write('%10.3E %10.3E\n'%(t[i], trate[i]))
    f.close()

def readcol(fn):
    a = array.array('f')
    b = array.array('f')
    f = open(fn, 'r')
    while (1):
        s = f.readline()
        if not s:
            break
        s = string.split(s)
        try:
            a.append(float(s[0]))
            b.append(float(s[1]))
        except:
            continue
    return a, b
        
def drs(outfile, nlow, ai, tai, te, trd, rd, wmin, wmax, dw):
    sigma = dw / const.FWHM
    sig8 = 8.0*sigma
    wavelength = array.array('f')
    spec = array.array('f')
    dwi = 0.1*sigma

    nb = len(ai)
    p = 1.0/(sqrt(2*pi)*sigma)

    si = array.array('f')
    for i in range(nb):
        si.append(p*ai[i]/(tai[i]+trd[i]))

    f = open(outfile, 'w')
    x = wmin
    while x < wmax:
        s = 0.0
        j0 = 0
        for i in range(nb):
            for j in range(j0, nlow[i]+j0):
                w = const.hc / te[j]
                d = x - w
                if (abs(d) > sig8):
                    continue
                
                d = d/sigma
                d = d*d*0.5
                s = s + si[i]*rd[j]*exp(-d)
        
            j0 = j0+nlow[i]

        if (s > 1E-10):
            wavelength.append(x)
            spec.append(s)
            f.write('%12.5E %10.3E\n'%(x, s))

        x = x + dwi
        
    f.close()

    return wavelength, spec

def constructfn(prefix, ext, n = -1, c = -1, i = 0):
    if n > 0:
        fn = '%s%02d'%(prefix, n)
    else:
        fn = '%s'%prefix
        
    if c >= 0:
        fn = '%s_%1d.%s'%(fn, c, ext)
    else:
        fn = '%s.%s'%(fn, ext)
    
    if i > 0:
        fn = '%s%d'%(fn, i)
        
    return fn

def drall(prefix, nrec, dn, f = 0, **p):
    if len(nrec) == 0:
        return
    
    if p.has_key('strength'):
        str = [0.0]*len(nrec)
    for k in range(len(nrec)):
        n = nrec[k]
        ai_file = constructfn(prefix, 'ai', n, dn)
        tr_file = constructfn(prefix, 'tr', n, dn)

        try:
            (channels, bound, ee, ai, tai) = read_ai(ai_file, f)
            (nlow, low, te, rd, trd, trdp) = read_tr(tr_file, bound)
        except IOError:
            continue

        br = branching(tai, trdp, trd) 

                
        if p.has_key('drs'):
            drs_file = constructfn(prefix, 'drs', n, dn, f)
            drs_p = p['drs']
            (w, s) = drs(drs_file, nlow, ai, tai, te, trd,
                         rd, drs_p[0], drs_p[1], drs_p[2])
        if p.has_key('rate'):
            rt_file = constructfn(prefix, 'rt', n, dn, f)
            rate_p = p['rate']
            (t, r) = drrate(rt_file, ee, ai, br,
                            rate_p[0], rate_p[1], rate_p[2])
        if p.has_key('cross'):
            cr_file = constructfn(prefix, 'cr', n, dn, f)
            cross_p = p['cross']
            (e, c) = drcross(cr_file, ee, ai, br, tai, trd,
                             cross_p[0], cross_p[1])
        if p.has_key('strength'):
            str[k] = drstrength(ee, ai, br)
            
    if p.has_key('strength'):
        str_file = constructfn(prefix, 'str', -1, dn)
        f = open(str_file, 'w')
        for i in range(len(nrec)):
            f.write('%-2d %15.8E\n'%(nrec[i], str[i]))
        f.close()


def dr(prefix, gr, rec_level, ground, ai_group, decay_group,
       channels, ngrid, nmin, nj0, max_kl, nspec, correlations,
       only_levels = 0):
    
    n_dn = range(len(channels))
    nrec = range(min(nmin), max(max(ngrid))+1)

    fac.SetRecSpectator(nspec[0], nspec[1])
    if len(nspec) == 3:
        ncorr = nspec[2]
    else:
        ncorr = nspec[0]
        
    ngr = len(gr)
    for n in nrec:
        
        print 'nrec = %d'%n
        lv_file = constructfn(prefix, 'lv', n)

        has_rec_states = 0

        if n < ncorr:
            for c in correlations:
                m = []
                for i in range(len(c)):
                    for j in n_dn:
                        if c[i] in decay_group[j]:
                            if nj0[j] == 0 or n <= nj0[j]:
                                m.append(max_kl[j])
                        if c[i] in ai_group[j]:
                            if n in ngrid[j]:
                                m.append(max_kl[j])
                if len(m) == 0:
                    continue
                mmax = max(m)
                fac.SetRecPWOptions(mmax, mmax)
                s = string.join(c, ', ')
                print "  recombined states for groups %s"%s
                fac.RecStates(n, c)
                has_rec_states = 1
        else:
            for i in range(ngr):
                m = []
                for j in n_dn:
                    if gr[i] in decay_group[j]:
                        if nj0[j] == 0 or n <= nj0[j]:
                            m.append(max_kl[j])
                    if gr[i] in ai_group[j]:
                        if n in ngrid[j]:
                            m.append(max_kl[j])
                        
                if len(m) == 0:
                    continue
                mmax = max(m)
                fac.SetRecPWOptions(mmax)
                print "  recombined states for group %s"%gr[i]
                fac.RecStates(n, [gr[i]])
                has_rec_states = 1
                
        if not has_rec_states:
            print '  No recombined states for this n'
            continue
    
        fac.LevelTable(lv_file, n)

        if only_levels:
            continue
    
        for i in n_dn:
            if not n in ngrid[i]:
                continue
        
            print "  DR channel %d"%channels[i]
        
            if nj0[i] > 0:
                nj = range(nmin[i], min([nj0[i], n])+1)
                if n > nj0[i]:
                    nj.append(n)
            else:
                nj = range(nmin[i], n+1)

            fac.SetRecPWLimits(0, max_kl[i])
            fac.SetPEGrid(0)
            ai_file = constructfn(prefix, 'ai', n, channels[i])
            tr_file = constructfn(prefix, 'tr', n, channels[i])
            fac.DRTable(ground[i], (ai_group[i], n), (decay_group[i], nj),
                        rec_level, ai_file, tr_file, channels[i])
            fac.FreeOrbitals(0)

        fac.FreeRecAngZ()
        if min(nj0) > 0 and n > max(nj0):
            fac.FreeOrbitals(n)
        
def getngrid(nopen, nmax):
    ngrid = []
    for i in range(len(nopen)):
        ngrid.append(nopen[i])
        if (i == len(nopen)-1):
            nlimit = nmax
        else:
            nlimit = nopen[i+1]

        dn = 1
        n = nopen[i]+1
        while n < nlimit:
            dn = dn*2
            ngrid.append(n)
            n = n+dn
            if n >= nmax:
                ngrid.append(nmax)
                return ngrid
            
        if nlimit < nmax:
            if ngrid[-1] != nlimit-1:
                ngrid.append(nlimit-1)
        else:
            ngrid.append(nlimit)
            break
        
    return ngrid
    
def extrapolate(n0, n1):
    f = 0.0
    a = float(n0)**3
    n0 = n0+1
    while (n0 <= n1):
        f = f + (1.0/n0)**3
        n0 = n0 + 1
    return f*a






