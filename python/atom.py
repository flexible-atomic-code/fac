from fac import *
from types import *
import time

# reinitialize FAC
def reinit(m = 0):
    Reinit(m)

    return

# the index of n-th shell in a complex
def get_index(n, complex):
    for i in range(len(complex)):
        if (complex[i][0] == n):
            return i
        elif (complex[i][0] > n):
            return -i-1
        
    return -i-2

def get_terms(n, nq):
    if (n == 1):
        return ['1s%d'%nq]
    elif (n == 2):
        if (nq == 1):
            return ['2s1', '2p1']
        elif (nq == 2):
            return ['2s2', '2s1 2p1', '2p2']
        else:
            a = '2s2 2p%d'%(nq-2)
            b = '2s1 2p%d'%(nq-1)
            return [a, b]
    elif (n == 3):
        if (nq == 1):
            return ['3s1', '3p1', '3d1']
        elif (nq == 2):
            return ['3s2', '3s1 3p1', '3s1 3d1', '3p2', '3p1 3d1', '3d2']
        elif (nq < 8):
            a = ['3s2 3p%d'%(nq-2)]
	    a.append('3s1 3p%d'%(nq-1))
            a.append('3s1 3p%d 3d1'%(nq-2))
            if (nq == 3):
	        a.append('3s2 3d1')
	    else:
                a.append('3s2 3p%d 3d1'%(nq-3))
            return a
        else:
	    if (nq != 8):
                a = ['3s2 3p6 3d%d'%(nq-8)]
	    else:
	        a = ['3s2 2p6']
            a.append('3s2 3p5 3d%d'%(nq-7))
            a.append('3s1 3p6 3d%d'%(nq-7))
            return a
    else:
        if (nq > 2):
            raise 'only n <= 3 supports nq > 2'
        else:
            return ['%d*%d'%(n, nq)]
    
# a configuration complex
class COMPLEX:
    def __init__(self, nele = 0):
        self.complex = []
        self.terms = []
        self.nrec = 0
        self.nrec_ext = 0
        self.name = ' '

        if (nele > 0):
            self.set_ground(nele)
            
        return

    def set_name(self, c = 0):
        if (type(c) == StringType):
            self.name = c
            return

        cname = 0
        for s in self.complex:
            if (type(cname) == IntType):
                cname = '%d*%d'%(s[0], s[1])
            else:
                cname = cname + ' %d*%d'%(s[0], s[1])
        if (type(cname) == IntType):
            cname = 'bare'
            
        self.name = cname
        
        return

    def set_terms(self):
        t = []
        for c in self.complex:
            t.append(get_terms(c[0], c[1]))
        if (len(t) == 0):
            self.terms = ['']
            return
        a = t[0]
        for b in t[1:]:
            c = []
            for i in a:
                for j in b:
                    c.append(i+' '+j)
            a = c
        self.terms = a

        return

    def set_ground(self, nelectrons):
        n = 1
        nele = nelectrons
        g = []
        while (nele > 0):
            nqm = 2*n*n
            if (nele < nqm):
                g.append((n, nele))
            else:
                g.append((n, nqm))
            nele = nele - nqm
            n = n+1

        self.complex = g
        self.set_name()
        self.set_terms()
        
        return 

    def set_excited(self, i0, i1, cbase):
        base = cbase.complex
        if (i0 >= i1):
            return -1
        i = get_index(i0, base)
        if (i < 0):
            return -1
        j = get_index(i1, base)
        ex = base[:]
        if (j >= 0):
            nq1 = ex[j][1]
            if (nq1 == 2*i1*i1):
                return 0
            else:
                ex[j] = (i1, nq1+1)
        else:
            j = -j - 1
            ex.insert(j, (i1, 1))

        nq0 = ex[i][1]
        if (nq0 == 1):
            ex.remove(ex[i])
        else:
            ex[i] = (i0, nq0-1)
    
        self.complex = ex
        self.set_name()
        self.set_terms()

        return 0

    def set_ionized(self, i0, cbase):
        base = cbase.complex
        i = get_index(i0, base)
        if (i < 0):
            return -1
        ion = base[:]
        nq0 = ion[i][1]
        if (nq0 == 1):
            ion.remove(ion[i])
        else:
            ion[i] = (i0, nq0-1)

        self.complex = ion
        self.set_name()
        self.set_terms()
        
        return 0

    def set_recombined(self, i0, cbase):
        self.set_name(cbase.name)
        self.nrec = i0
        self.complex = ([cbase.name], i0)
        
        return 0
    

# a group of complexes
class CGROUP:
    def __init__(self, type):
        self.cgroup = []
        self.type = type

        return

    def add_complex(self, complex):
        self.cgroup.append(complex)

        return


# a generic atomic ion
class ATOM:    
    def __init__(self, nele, asym = 0):
	self.process = {'ce': 1, 
		        'tr': 1, 
			'rr': 1, 
			'ci': 1,
			'ai': 1}

        self.nele = nele
        self.nele_max = [0, 2, 10, 28]
        
        self.grd_complex = COMPLEX(nele)
        self.exc_complex = []
        self.ion_complex = CGROUP('ion')

        self.angz_cut0 = 1E-3
        self.tr_cut0 = 1E-3
        self.ai_cut0 = 1E-3
        if (self.nele <= self.nele_max[1]):
            self.n_shells = 1
            self.nexc_max = [7, 7, 5]
            self.nexc_rec = [0, 0, 0]
            self.nrec_max = [[10,0], [8,10], [6,6]]
            self.nrec_ext = [25, 25, 6]
            self.rec_pw_max = [9, 6, 6]
            self.n_decay = [[10,10], [-1,-1], [-1,-1]]
            self.angz_cut1 = 1E-3
            self.tr_cut1 = 1E-3
            self.ai_cut1 = 1E-3
            self.angz_cut2 = 1E-3
            self.tr_cut2 = 1E-3
            self.ai_cut2 = 1E-3
        elif (self.nele <= self.nele_max[2]):
            self.n_shells = 2
            self.nexc_max = [4, 4, 4, 4]
            self.nexc_rec = [7, 0, 0, 0]
            self.nrec_max = [[10,0], [6,6], [6,10], [5,5]]
            self.nrec_ext = [25, 10, 25, 5]
            self.rec_pw_max = [9, 6, 6, 6]
            self.n_decay = [[10,10], [-1,-1], [-1,-1], [-1,-1]]
            self.angz_cut1 = 1E-2
            self.tr_cut1 = 1E-2
            self.ai_cut1 = 1E-2
            self.angz_cut2 = 1E-1
            self.tr_cut2 = 5E-2
            self.ai_cut2 = 5E-2
        elif (self.nele <= self.nele_max[3]):
            self.n_shells = 3
            self.nexc_max = [4, 4, 4, 4, 0]
            self.nexc_rec = [0, 0, 0, 0, 0]
            self.nrec_max = [[10,0], [6,10], [6,6], [5,5], [5,5]]
            self.nrec_ext = [25, 25, 10, 10, 5]
            self.rec_pw_max = [9, 6, 6, 6, 6]
            self.n_decay = [[10,9], [-1,-1], [-1,-1], [-1,-1], [-1,-1]]
            self.angz_cut1 = 1E-2
            self.tr_cut1 = 1E-2
            self.ai_cut1 = 1E-2
            self.angz_cut2 = 1E-1
            self.tr_cut2 = 5E-2
            self.ai_cut2 = 5E-2
        else:
            raise 'ion with NELE >= %d not supported'%(self.nele_max[3])

        if (type(asym) == StringType):
            self.set_atom(asym)
            
        return

    
    def set_atom(self, asym):
        self.asym = asym
        self.bfiles = {'en': '%s%02db.en'%(asym, self.nele),
                       'tr': '%s%02db.tr'%(asym, self.nele),
                       'ce': '%s%02db.ce'%(asym, self.nele),
                       'rr': '%s%02db.rr'%(asym, self.nele),
                       'ai': '%s%02db.ai'%(asym, self.nele),
                       'ci': '%s%02db.ci'%(asym, self.nele)}
        
        self.afiles = {'en': '%s%02da.en'%(asym, self.nele),
                       'tr': '%s%02da.tr'%(asym, self.nele),
                       'ce': '%s%02da.ce'%(asym, self.nele),
                       'rr': '%s%02da.rr'%(asym, self.nele),
                       'ai': '%s%02da.ai'%(asym, self.nele),
                       'ci': '%s%02da.ci'%(asym, self.nele)}

        return


    def add_excited(self, n0, n1, ibase):
        if (ibase == -1):
            base = self.grd_complex
        else:
            base = self.exc_complex[0].cgroup[ibase]
            
        n1.sort()
        
        cg = CGROUP('exc')
        for i1 in n1:
            ex = COMPLEX()
            if (ex.set_excited(n0, i1, base) == 0):
                cg.add_complex(ex)

        self.exc_complex.append(cg)
        
        return
    

    def add_ionized(self, n0, ibase):
        if (ibase == -1):
            base = self.grd_complex
        else:
            base = self.exc_complex[0].cgroup[ibase]
        
        ion = COMPLEX()
        ion.set_ionized(n0, base)
        self.ion_complex.add_complex(ion)

        return


    def add_recombined(self, n0, ibase):
        base = self.ion_complex.cgroup[ibase]
        n0.sort()
        cg = self.exc_complex[ibase]
        for i0 in n0:
            rec = COMPLEX()
            if (rec.set_recombined(i0, base) == 0):
                cg.add_complex(rec)
                
        return


    def run_tr_ce(self, a, b, c, d, tr = [-1], ce = 1):
        if (type(a) == StringType or type(a) == IntType):
            low = [a]
        else:
            low = a
        if (type(b) == StringType or type(b) == IntType):
            up = [b]
        else:
            up = b
            
	if (self.process['tr'] != 0):
            for m in tr:
                s = 'TR: %s -> %s %d'%(str(a), str(b), m)
                Print(s)
                TransitionTable(self.bfiles['tr'], low, up, m)

        if (ce != 0 and self.process['ce'] != 0):
            if (type(c) == StringType or type(c) == IntType):
                low = [c]
            else:
                low = c
            if (type(d) == StringType or type(d) == IntType):
                up = [d]
            else:
                up = d
            s = 'CE: %s -> %s'%(str(c), str(d))
            Print(s)
            CETable(self.bfiles['ce'], low, up)

        Reinit(radial = 1, excitation = 1, structure = 1)
        
        return


    def run_ci_rr(self, a, b, c, d, ci = 1, rr = 1):
        if (rr != 0 and self.process['rr'] != 0):
            if (type(c) == StringType or type(c) == IntType):
                low = [c]
            else:
                low = c
            if (type(d) == StringType or type(d) == IntType):
                up = [d]
            else:
                up = d
            s = 'RR: %s -> %s'%(c, d)
            Print(s)
            
            RRTable(self.bfiles['rr'], low, up)
            
            SetPEGrid(0)
            SetUsrPEGrid(0)
            SetRRTEGrid(0)
            
        if (ci != 0 and self.process['ci'] != 0):
            if (type(a) == StringType or type(a) == IntType):
                low = [a]
            else:
                low = a
            if (type(b) == StringType or type(b) == IntType):
                up = [b]
            else:
                up = b
            s = 'CI: %s -> %s'%(a, b)
            Print(s)
            CITable(self.bfiles['ci'], low, up)


        if (ci != 0 or rr != 0):
            Reinit(radial = 1, recombination = 1,
                   ionization = 1, structure = 1)
        
        return


    def run_ai(self, a, b, k):
	if (self.process['ai'] == 0):
	    return

        s = 'AI: %s -> %s %d'%(a, b, k)
        Print(s)
        if (type(a) == StringType or type(a) == IntType):
            low = [a]
        else:
            low = a
        if (type(b) == StringType or type(b) == IntType):
            up = [b]
        else:
            up = b

        AITable(self.bfiles['ai'], low, up, k)

        Reinit(radial = 1, recombination = 1, structure = 1)
        
        return
    

    def run_en(self):        
        SetAtom(self.asym)
        SetAngZCut(self.angz_cut0)

        g = [self.grd_complex.name]
        OptimizeRadial(g)

        # ground and the first excited complex are interacting
        Print('Structure: ground complex')
        c = self.exc_complex[0].cgroup[0].name
        g.append(c)
        Structure(self.bfiles['en'], g)
        
        Print('Structure: ionized complexes')
        c = self.ion_complex.cgroup
        for a in c:
            Structure(self.bfiles['en'], [a.name])

        Print('Structure: excited complexes')
        for i in range(len(self.exc_complex)):
            c = self.exc_complex[i].cgroup
            SetRecPWOptions(self.rec_pw_max[i])
            SetRecSpectator(self.nexc_max[i]+1)
            for j in range(len(c)):
                a = c[j]
                if (i == 0 and j == 0):
                    continue
                if (a.nrec > 0):
                    RecStates(self.bfiles['en'], [a.name], a.nrec)
                else:
                    if (a.name != ' '):
                        Structure(self.bfiles['en'], [a.name])

        return


    def run_tr_ce_all(self):
        g = [self.grd_complex.name]
        SetAngZCut(self.angz_cut0)
        SetTransitionCut(self.tr_cut0)
        SetAICut(self.ai_cut0)
        self.run_tr_ce(g, g, g, g, tr=[-1,1,-2,2], ce=1)
        for i in range(self.n_shells):
            c = self.exc_complex[i].cgroup
            for j in range(len(c)):
                if (c[j].name == ' '):
                    continue
                a1 = g
                if (c[j].nrec == 0):
                    b = [c[j].name]
                else:
                    b = ([c[j].name], c[j].nrec)
                if (i == 0 and j == 0):
                    a2 = g
                else:
                    a2 = [0]
                if (i == 0 and j == 0):
                    tr = [-1, 1, -2, 2]
                    SetAngZCut(self.angz_cut0)
                    SetTransitionCut(self.tr_cut0)
                else:
                    tr = [-1]
                    SetAngZCut(self.angz_cut1)
                    SetTransitionCut(self.tr_cut1)
                if (c[j].nrec <= self.nexc_rec[i]):
                    ce = 1
                else:
                    ce = 0
                self.run_tr_ce(a1, b, a2, b, tr=tr, ce=ce)

        if (self.nele < 3):
            ce0 = 1
        else:
            ce0 = 0
        for i in range(len(self.exc_complex)):
            c = self.exc_complex[i].cgroup
            for j in range(len(c)):
                if (c[j].name == ' '):
                    continue
                if (c[j].nrec == 0):
                    b = [c[j].name]
                    SetAngZCut(self.angz_cut1)
                    SetTransitionCut(self.tr_cut1)
                else:
                    b = ([c[j].name], c[j].nrec)
                    SetAngZCut(self.angz_cut2)
                    SetTransitionCut(self.tr_cut2)
                if (i == 0 and j == 0):
                    ce = ce0
                else:
                    ce = 0
                for m in range(j+1):
                    if (c[m].name == ' '):
                        continue
                    if (c[j].nrec <= self.n_decay[i][0] and
                        c[m].nrec <= self.n_decay[i][1]):
                        if (i == 0 and j == 0):
                            tr = [-1, 1]
                        else:
                            tr = [-1]
                    else:
                        tr = []
                    if (c[m].nrec == 0):
                        a = [c[m].name]
                    else:
                        a = ([c[m].name], c[m].nrec)
                    if (ce != 0 or tr != []):
                        self.run_tr_ce(a, b, a, b, tr=tr, ce=ce)

                if (i == 0):
                    continue
                d = self.exc_complex[0].cgroup
                for m in range(len(d)):
                    if (d[m].name == ' '):
                        continue
                    n1 = d[m].nrec
                    if (n1 == 0):
                        a = [d[m].name]
                        n1 = d[m].complex[-1][0]
                        SetAngZCut(self.angz_cut1)
                        SetTransitionCut(self.tr_cut1)
                    else:
                        a = ([d[m].name], n1)
                        SetAngZCut(self.angz_cut2)
                        SetTransitionCut(self.tr_cut2)
                    n2 = c[j].nrec
                    if (n2 == 0):
                        n2 = c[j].complex[-1][0]
                    if (n1 != n2):
                        continue
                    self.run_tr_ce(a, b, a, b, tr=[-1], ce=0)
                        
        return

    
    def run_ci_rr_all(self):
        SetAngZCut(self.angz_cut0)
        g = [self.grd_complex.name]
        for i in range(self.n_shells):
            b = [self.ion_complex.cgroup[i].name]
            self.run_ci_rr(g, b, g, b, ci=1, rr=1)

        b = [self.ion_complex.cgroup[0].name]
        for i in range(len(self.exc_complex[0].cgroup)):
            a = self.exc_complex[0].cgroup[i].name
            if (a == ' '):
                continue
            a = [a]
            n = self.exc_complex[0].cgroup[i].nrec
            if (n > 0):
                a = (a, n)
            self.run_ci_rr(a, b, a, b, ci=1, rr=1)
             
        return

    
    def run_ai_all(self):
        for i in range(len(self.exc_complex)):
            c = self.exc_complex[i].cgroup
            if (i == 0):
                if (self.nele == self.nele_max[self.n_shells-1]+1):
                    continue
            for j in range(len(c)):
                a = c[j].name
                if (a == ' '):
                    continue
                a = [a]
                n = c[j].nrec
                if (n > 0):
                    a = (a, n)
                    SetAngZCut(self.angz_cut1)
                    SetAICut(self.ai_cut1)
                else:
                    SetAngZCut(self.angz_cut0)
                    SetAICut(self.ai_cut0)
                if (i == 0):
                    p = 1
                else:
                    p = 2
                for k in range(p):
                    if (n > self.nrec_max[i][k] and
                        c[j].nrec_ext == 0):
                        continue
                    if (self.nele == self.nele_max[self.n_shells-1]+1):
                        t = k
                    else:
                        if (k == 1 and i >= self.n_shells):
                            t = self.n_shells
                        else:
                            t = k
                    b = [self.ion_complex.cgroup[t].name]
                    m = k*10+i
                    if (c[j].nrec_ext > 0):
                        m = n + (self.nrec_max[i][k]+1)*1000
                        m = -m
                    self.run_ai(a, b, m)
                    
        return
    
    
    def print_table(self, type = [], v = 0):
        if (v != 0):
            MemENTable(self.bfiles['en'])

        for t in type:
            PrintTable(self.bfiles[t], self.afiles[t], v)
            
        return
  

    def set_configs(self):
        n = range(1, self.n_shells+1)
        n.reverse()
        nexc = self.nexc_max
        nrec = []
        for i in self.nrec_max:
            nrec.append(max(i))
        for i in n:
            j = self.n_shells - i
            self.add_excited(i, range(nexc[j]+1), -1)
            self.add_ionized(i, -1)
            r = range(nexc[j]+1, nrec[j]+1)
            if (self.nrec_ext[j] > nrec[j]):
                r.append(self.nrec_ext[j])
            self.add_recombined(r, j)
            self.exc_complex[j].cgroup[-1].nrec_ext = nrec[j]+1

        j = self.n_shells
        if (self.nele == self.nele_max[self.n_shells-1]+1):
            i = self.n_shells-1
            mrec = nrec[j+1]
            r = range(nexc[j]+1, mrec+1)
            if (self.nrec_ext[j+1] > mrec):
                r.append(self.nrec_ext[j+1])
        else:
            i = self.n_shells
            mrec = nrec[j]
            r = range(nexc[j]+1, mrec+1)
            if (self.nrec_ext[j] > mrec):
                r.append(self.nrec_ext[j])
        if (i > 0):
            self.add_excited(i, range(j+1, nexc[j]+1), 0)
            self.add_ionized(i, 0)
            self.add_recombined(r, j)
            self.exc_complex[j].cgroup[-1].nrec_ext = mrec+1

        if (self.nele > self.nele_max[self.n_shells-1]+1):
            i = self.n_shells
            j = self.n_shells+1
            self.add_excited(i, range(j+1, nexc[j]+1), 1)
            self.add_ionized(i, 1)
            r = range(nexc[j]+1, nrec[j]+1)
            if (self.nrec_ext[j] > nrec[j]):
                r.append(self.nrec_ext[j])
            self.add_recombined(r, j)
            self.exc_complex[j].cgroup[-1].nrec_ext = nrec[j]+1
            
        # set up configurations
        c = self.grd_complex.name
        s = 'adding complex: %s'%(str(c))
        Print(s)
        for t in self.grd_complex.terms:
            Config(t, group = c)
            
        for cg in self.exc_complex:
            for i in range(len(cg.cgroup)):
                if (cg.cgroup[i].nrec > 0):
                    continue
                c = cg.cgroup[i].name
                if (c == ' '):
                    continue
                s = 'adding complex: %s'%(str(c))
                Print(s)
                for t in cg.cgroup[i].terms:
                    Config(t, group = c)

        cg = self.ion_complex.cgroup
        for i in range(len(cg)):
            c = cg[i].name
            s = 'adding complex: %s'%(str(c))
            Print(s)
            if (c == 'bare'):
                Config('', group = c)
            else:
                for t in cg[i].terms:
                    Config(t, group = c)
 
        return


    def run_fac(self):
        start = time.time()
        start0 = start

        # solve the structure
        s = 'EN...'
        Print(s)
        self.run_en()        
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print s

        # transition rates and collisional excitation
        start = stop  
        s = 'TR and/or CE...'
        Print(s)  
        self.run_tr_ce_all() 
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print s

        # CI and/or RR
        start = stop
        s = 'CI and/or RR...'
        Print(s)
        self.run_ci_rr_all()
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print s 

        start = stop 
        s = 'AI...'
        Print(s)
        self.run_ai_all()
        stop = time.time()
        s = 'Done %10.3E s'%(stop-start)
        print s

        s = 'Total Time: %10.3E s'%(stop-start0)
        print s

        return


def atomic_data(nele, asym, iprint = -1, **process):
    if (type(asym) == StringType):
        a = [asym]
    else:
        a = asym

    if (type(nele) == IntType):
        n = [nele]
    else:
        n = nele

    for m in n:
        atom = ATOM(m)
	if (process.has_key('no_ce')):
	    atom.process['ce'] = not process['no_ce']
	if (process.has_key('no_tr')):
	    atom.process['tr'] = not process['no_tr']
	if (process.has_key('no_rr')):
	    atom.process['rr'] = not process['no_rr']
	if (process.has_key('no_ci')):
	    atom.process['ci'] = not process['no_ci']
	if (process.has_key('no_ai')):
	    atom.process['ai'] = not process['no_ai']

        atom.set_configs()
        for b in a:
            s = 'NELE = %d'%m
            Print(s)
            s = 'ATOM = %s'%b
            z = ATOMICSYMBOL.index(b.lower())
            s = '%s, Z = %d'%(s, z)
            if (z < 1.5*m):
                Print('%s, skipping'%s)
                continue
            Print(s)
            atom.set_atom(b)
            atom.run_fac()
            if (iprint >= 0):
                atom.print_table(['en', 'tr', 'ce', 'rr', 'ai', 'ci'], iprint)
                
            Reinit(1)

        Reinit(0)

    return
